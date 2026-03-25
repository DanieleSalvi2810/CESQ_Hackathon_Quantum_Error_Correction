#!/usr/bin/env python3
import json
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx
    import numpy as np
    from pymatching import Matching
except ModuleNotFoundError as exc:
    print(
        f"Missing dependency: {exc.name}. Install with: python3 -m pip install --user pymatching",
        file=sys.stderr,
    )
    raise SystemExit(3)


def wrap_index(index: int, size: int) -> int:
    return index % size


def flatten_matrix(matrix: list[list[int]]) -> np.ndarray:
    return np.array([value & 1 for row in matrix for value in row], dtype=np.uint8)


def validate_syndromes(syndrome_plaquette: list[list[int]], syndrome_cross: list[list[int]]) -> int:
    if not syndrome_plaquette or not syndrome_plaquette[0]:
        raise ValueError("syndrome_plaquette is empty")
    if not syndrome_cross or not syndrome_cross[0]:
        raise ValueError("syndrome_cross is empty")

    n_rows = len(syndrome_plaquette)
    n_cols = len(syndrome_plaquette[0])

    if len(syndrome_cross) != n_rows:
        raise ValueError("syndrome_cross row count mismatch")
    if any(len(row) != n_cols for row in syndrome_plaquette):
        raise ValueError("syndrome_plaquette is ragged")
    if any(len(row) != n_cols for row in syndrome_cross):
        raise ValueError("syndrome_cross is ragged")
    if n_rows != n_cols:
        raise ValueError("expected square syndrome matrices (D x D)")

    return n_rows


def build_x_check_matrix(d: int) -> np.ndarray:
    n_checks = d * d
    n_qubit_rows = 2 * d
    n_qubit_cols = d
    n_qubits = n_qubit_rows * n_qubit_cols
    h_x = np.zeros((n_checks, n_qubits), dtype=np.uint8)

    for i in range(n_qubit_rows):
        for j in range(n_qubit_cols):
            qubit_index = i * n_qubit_cols + j
            if i % 2 == 0:
                toggles = [(i // 2, j), (i // 2, j + 1)]
            else:
                toggles = [((i - 1) // 2, j), ((i + 1) // 2, j)]

            for row, col in toggles:
                row_wrapped = wrap_index(row, d)
                col_wrapped = wrap_index(col, d)
                check_index = row_wrapped * d + col_wrapped
                h_x[check_index, qubit_index] ^= 1

    return h_x


def build_z_check_matrix(d: int) -> np.ndarray:
    n_checks = d * d
    n_qubit_rows = 2 * d
    n_qubit_cols = d
    n_qubits = n_qubit_rows * n_qubit_cols
    h_z = np.zeros((n_checks, n_qubits), dtype=np.uint8)

    for i in range(n_qubit_rows):
        for j in range(n_qubit_cols):
            qubit_index = i * n_qubit_cols + j
            if i % 2 == 0:
                toggles = [(i // 2 - 1, j), (i // 2, j)]
            else:
                toggles = [((i - 1) // 2, j - 1), ((i - 1) // 2, j)]

            for row, col in toggles:
                row_wrapped = wrap_index(row, d)
                col_wrapped = wrap_index(col, d)
                check_index = row_wrapped * d + col_wrapped
                h_z[check_index, qubit_index] ^= 1

    return h_z


def make_orientation_weights(d: int, horizontal_weight: float, vertical_weight: float) -> np.ndarray:
    n_qubit_rows = 2 * d
    n_qubit_cols = d
    n_qubits = n_qubit_rows * n_qubit_cols
    weights = np.empty(n_qubits, dtype=float)

    for i in range(n_qubit_rows):
        row_weight = horizontal_weight if i % 2 == 0 else vertical_weight
        for j in range(n_qubit_cols):
            qubit_index = i * n_qubit_cols + j
            weights[qubit_index] = row_weight

    return weights


def print_channel_pairing(label: str, matching: Matching, syndrome: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    matched_pairs = matching.decode_to_matched_dets_array(syndrome)
    matched_edges = matching.decode_to_edges_array(syndrome)
    print(f"{label} matched detector pairs: {matched_pairs.tolist()}")
    print(f"{label} path edges: {matched_edges.tolist()}")
    return matched_pairs, matched_edges


def make_positions(graph: nx.Graph, d: int) -> dict[int, tuple[float, float]]:
    pos: dict[int, tuple[float, float]] = {}
    overflow_x = float(d + 1)
    overflow_y = 0.0
    for node in graph.nodes():
        node_int = int(node)
        if 0 <= node_int < d * d:
            row = node_int // d
            col = node_int % d
            pos[node_int] = (float(col), -float(row))
        else:
            pos[node_int] = (overflow_x, -overflow_y)
            overflow_y += 1.0
    return pos


def edges_from_fault_ids(graph: nx.Graph, fault_ids: set[int]) -> list[tuple[int, int]]:
    marked_edges: list[tuple[int, int]] = []
    for u, v, attrs in graph.edges(data=True):
        edge_fault_ids = attrs.get("fault_ids", set())
        edge_fault_ids_int = {int(fid) for fid in edge_fault_ids}
        if edge_fault_ids_int & fault_ids:
            marked_edges.append((int(u), int(v)))
    return marked_edges


def normalize_edge(u: int, v: int) -> tuple[int, int]:
    return (u, v) if u <= v else (v, u)


def edge_set_from_array(edges: np.ndarray) -> set[tuple[int, int]]:
    out: set[tuple[int, int]] = set()
    for edge in edges.tolist():
        if len(edge) != 2:
            continue
        out.add(normalize_edge(int(edge[0]), int(edge[1])))
    return out


def split_paths_from_matched_edges(matched_edges: np.ndarray) -> list[list[tuple[int, int]]]:
    path_graph = nx.Graph()
    for edge in matched_edges.tolist():
        if len(edge) != 2:
            continue
        path_graph.add_edge(int(edge[0]), int(edge[1]))

    path_components: list[list[tuple[int, int]]] = []
    for component_nodes in nx.connected_components(path_graph):
        sub = path_graph.subgraph(component_nodes)
        path_components.append([(int(u), int(v)) for u, v in sub.edges()])
    return path_components


def draw_nodes_with_syndrome(
    ax: plt.Axes,
    graph: nx.Graph,
    pos: dict[int, tuple[float, float]],
    syndrome: np.ndarray,
) -> None:
    active_nodes = {int(i) for i, bit in enumerate(syndrome.tolist()) if int(bit) == 1}
    node_colors = []
    for node, attrs in graph.nodes(data=True):
        node_int = int(node)
        is_boundary = bool(attrs.get("is_boundary", False))
        if node_int in active_nodes:
            node_colors.append("#D1495B")
        elif is_boundary:
            node_colors.append("#FFFFFF")
        else:
            node_colors.append("#8AB6D6")

    nx.draw_networkx_nodes(
        graph,
        pos=pos,
        ax=ax,
        node_color=node_colors,
        edgecolors="#2F2F2F",
        linewidths=0.6,
        node_size=220,
    )
    nx.draw_networkx_labels(graph, pos=pos, ax=ax, font_size=7)


def draw_channel(
    ax: plt.Axes,
    label: str,
    matching: Matching,
    syndrome: np.ndarray,
    matched_pairs: np.ndarray,
    matched_edges: np.ndarray,
    error_fault_ids: set[int],
    show_paths: bool,
    show_error_edges: bool,
    d: int,
) -> None:
    graph = matching.to_networkx()
    pos = make_positions(graph, d)

    nx.draw_networkx_edges(graph, pos=pos, ax=ax, edge_color="#C7CBD1", width=1.0)

    if show_paths:
        path_components = split_paths_from_matched_edges(matched_edges)
        palette = ["#D1495B", "#2A9D8F", "#4361EE", "#F4A261", "#6A4C93", "#3A86FF", "#E63946", "#4CC9F0"]
        for idx, path_edges in enumerate(path_components):
            if not path_edges:
                continue
            nx.draw_networkx_edges(
                graph,
                pos=pos,
                ax=ax,
                edgelist=path_edges,
                edge_color=palette[idx % len(palette)],
                width=3.0,
            )

    error_edges = edges_from_fault_ids(graph, error_fault_ids)
    if show_error_edges and error_edges:
        nx.draw_networkx_edges(
            graph,
            pos=pos,
            ax=ax,
            edgelist=error_edges,
            edge_color="#111111",
            width=5.0,
        )

    draw_nodes_with_syndrome(ax, graph, pos, syndrome)

    error_edges_count = len(error_edges)
    if show_paths:
        ax.set_title(f"{label}: pairs={len(matched_pairs)}, path_edges={len(matched_edges)}", fontsize=10)
    else:
        ax.set_title(f"{label}: flipped_qubit_edges={error_edges_count}", fontsize=10)

    legend_items = []
    if show_paths:
        legend_items.append(plt.Line2D([0], [0], color="#D1495B", lw=3, label="Matching path (colored)"))
    if show_error_edges:
        legend_items.append(plt.Line2D([0], [0], color="#111111", lw=5, label="Physical error edge"))
    if legend_items:
        ax.legend(handles=legend_items, loc="upper right", fontsize=8, frameon=True)
    ax.set_axis_off()


def draw_xor_channel(
    ax: plt.Axes,
    label: str,
    matching: Matching,
    syndrome: np.ndarray,
    matched_edges: np.ndarray,
    error_fault_ids: set[int],
    d: int,
) -> None:
    graph = matching.to_networkx()
    pos = make_positions(graph, d)
    nx.draw_networkx_edges(graph, pos=pos, ax=ax, edge_color="#D5D8DC", width=1.0)

    path_edges = edge_set_from_array(matched_edges)
    flip_edges = {normalize_edge(u, v) for u, v in edges_from_fault_ids(graph, error_fault_ids)}
    xor_edges = path_edges.symmetric_difference(flip_edges)

    # Keep only closed-loop components: every node in the component must have degree 2.
    xor_graph = nx.Graph()
    xor_graph.add_edges_from(xor_edges)
    loop_components: list[list[tuple[int, int]]] = []
    for component_nodes in nx.connected_components(xor_graph):
        sub = xor_graph.subgraph(component_nodes)
        if all(deg == 2 for _, deg in sub.degree()):
            loop_components.append([(int(u), int(v)) for u, v in sub.edges()])

    palette = ["#D1495B", "#2A9D8F", "#4361EE", "#F4A261", "#6A4C93", "#3A86FF", "#E63946", "#4CC9F0"]
    n_path_edges = 0
    n_flip_edges = 0
    for idx, component_edges in enumerate(loop_components):
        path_part = []
        flip_part = []
        for u, v in component_edges:
            key = normalize_edge(u, v)
            if key in path_edges:
                path_part.append((u, v))
            elif key in flip_edges:
                flip_part.append((u, v))
        if path_part:
            n_path_edges += len(path_part)
            nx.draw_networkx_edges(
                graph,
                pos=pos,
                ax=ax,
                edgelist=path_part,
                edge_color=palette[idx % len(palette)],
                width=3.4,
            )
        if flip_part:
            n_flip_edges += len(flip_part)
            nx.draw_networkx_edges(
                graph,
                pos=pos,
                ax=ax,
                edgelist=flip_part,
                edge_color="#111111",
                width=5.0,
            )

    draw_nodes_with_syndrome(ax, graph, pos, syndrome)
    ax.set_title(
        f"{label}: loops={len(loop_components)}, path_edges={n_path_edges}, flip_edges={n_flip_edges}",
        fontsize=10,
    )
    legend_items = [
        plt.Line2D([0], [0], color="#D1495B", lw=3, label="Loop edge from path"),
        plt.Line2D([0], [0], color="#111111", lw=5, label="Loop edge from physical error"),
    ]
    ax.legend(handles=legend_items, loc="upper right", fontsize=8, frameon=True)
    ax.set_axis_off()


def save_pairing_overview_figure(
    d: int,
    matching_x: Matching,
    syndrome_x: np.ndarray,
    pairs_x: np.ndarray,
    edges_x: np.ndarray,
    error_fault_ids_x: set[int],
    matching_z: Matching,
    syndrome_z: np.ndarray,
    pairs_z: np.ndarray,
    edges_z: np.ndarray,
    error_fault_ids_z: set[int],
    use_weights: bool,
) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(20, 12), constrained_layout=True)
    flat_axes = axes.ravel()
    draw_channel(
        flat_axes[0], "X flips", matching_x, syndrome_x, pairs_x, edges_x, error_fault_ids_x, False, True, d
    )
    draw_channel(
        flat_axes[1], "X paths", matching_x, syndrome_x, pairs_x, edges_x, error_fault_ids_x, True, False, d
    )
    draw_xor_channel(flat_axes[2], "X XOR edge-wise", matching_x, syndrome_x, edges_x, error_fault_ids_x, d)
    draw_channel(flat_axes[3], "Z flips", matching_z, syndrome_z, pairs_z, edges_z, error_fault_ids_z, False, True, d)
    draw_channel(flat_axes[4], "Z paths", matching_z, syndrome_z, pairs_z, edges_z, error_fault_ids_z, True, False, d)
    draw_xor_channel(flat_axes[5], "Z XOR edge-wise", matching_z, syndrome_z, edges_z, error_fault_ids_z, d)

    mode_tag = "weighted" if use_weights else "unweighted"
    fig.suptitle(f"Pairing Overview ({mode_tag})", fontsize=12)
    output_path = Path("pairing_overview.png")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def main() -> int:
    if len(sys.argv) not in (3, 5):
        print(
            "Usage: python3 print_pairings_with_pymatching.py <syndrome.json> <error_matrix.json> "
            "[<horizontal_weight> <vertical_weight>]",
            file=sys.stderr,
        )
        return 1

    syndrome_path = Path(sys.argv[1])
    error_matrix_path = Path(sys.argv[2])
    use_weights = len(sys.argv) == 5

    horizontal_weight = 1.0
    vertical_weight = 1.0
    if use_weights:
        try:
            horizontal_weight = float(sys.argv[3])
            vertical_weight = float(sys.argv[4])
        except ValueError:
            print("Weights must be valid floats.", file=sys.stderr)
            return 1
        if horizontal_weight <= 0 or vertical_weight <= 0:
            print("Weights must be > 0.", file=sys.stderr)
            return 1

    try:
        with syndrome_path.open("r", encoding="utf-8") as f:
            syndrome_data = json.load(f)
        with error_matrix_path.open("r", encoding="utf-8") as f:
            error_matrix = json.load(f)

        syndrome_plaquette = syndrome_data["syndrome_plaquette"]
        syndrome_cross = syndrome_data["syndrome_cross"]
        d = validate_syndromes(syndrome_plaquette, syndrome_cross)

        n_qubit_rows = 2 * d
        n_qubit_cols = d
        if (
            len(error_matrix) != n_qubit_rows
            or any(not isinstance(row, list) for row in error_matrix)
            or any(len(row) != n_qubit_cols for row in error_matrix)
        ):
            raise ValueError("error_matrix shape mismatch (expected 2D x D)")

        syndrome_x = flatten_matrix(syndrome_plaquette)
        syndrome_z = flatten_matrix(syndrome_cross)

        h_x = build_x_check_matrix(d)
        h_z = build_z_check_matrix(d)

        if use_weights:
            weights = make_orientation_weights(d, horizontal_weight, vertical_weight)
            matching_x = Matching(h_x, weights=weights)
            matching_z = Matching(h_z, weights=weights)
            print(
                f"Pairing view (weighted): horizontal_weight={horizontal_weight}, "
                f"vertical_weight={vertical_weight}"
            )
        else:
            matching_x = Matching(h_x)
            matching_z = Matching(h_z)
            print("Pairing view (unweighted)")

        error_fault_ids_x: set[int] = set()
        error_fault_ids_z: set[int] = set()
        for i in range(n_qubit_rows):
            for j in range(n_qubit_cols):
                value = int(error_matrix[i][j])
                qid = i * n_qubit_cols + j
                if value in (1, 3):
                    error_fault_ids_x.add(qid)
                if value in (2, 3):
                    error_fault_ids_z.add(qid)

        pairs_x, edges_x = print_channel_pairing("X", matching_x, syndrome_x)
        pairs_z, edges_z = print_channel_pairing("Z", matching_z, syndrome_z)
        overview_path = save_pairing_overview_figure(
            d,
            matching_x,
            syndrome_x,
            pairs_x,
            edges_x,
            error_fault_ids_x,
            matching_z,
            syndrome_z,
            pairs_z,
            edges_z,
            error_fault_ids_z,
            use_weights,
        )
        print(f"Saved pairing overview image to {overview_path}")
    except Exception as exc:
        print(f"Pairing extraction failed: {exc}", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
