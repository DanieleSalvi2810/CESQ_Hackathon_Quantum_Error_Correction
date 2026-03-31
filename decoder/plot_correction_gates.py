#!/usr/bin/env python3
import json
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    import networkx as nx
    import numpy as np
    from pymatching import Matching
except ModuleNotFoundError as exc:
    print(
        f"Missing dependency: {exc.name}. Install with: python3 -m pip install --user matplotlib networkx numpy pymatching",
        file=sys.stderr,
    )
    raise SystemExit(3)


def wrap_index(index: int, size: int) -> int:
    return index % size


def validate_correction_matrix(correction_matrix: list[list[int]]) -> int:
    if not correction_matrix or not correction_matrix[0]:
        raise ValueError("correction_matrix is empty")

    n_rows = len(correction_matrix)
    n_cols = len(correction_matrix[0])

    if any(len(row) != n_cols for row in correction_matrix):
        raise ValueError("correction_matrix is ragged")
    if n_rows != 2 * n_cols:
        raise ValueError("expected correction_matrix shape (2D x D)")

    for row in correction_matrix:
        for value in row:
            if int(value) not in (0, 1, 2, 3):
                raise ValueError("correction_matrix values must be in {0,1,2,3}")

    return n_cols


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


def lattice_coords(node: int, d: int) -> tuple[int, int] | None:
    if 0 <= node < d * d:
        return (node // d, node % d)
    return None


def torus_edge_segments(
    u: int,
    v: int,
    pos: dict[int, tuple[float, float]],
    d: int,
    stub: float = 0.45,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    if u not in pos or v not in pos:
        return []

    rc_u = lattice_coords(u, d)
    rc_v = lattice_coords(v, d)
    if rc_u is None or rc_v is None:
        return [(pos[u], pos[v])]

    ru, cu = rc_u
    rv, cv = rc_v

    if ru == rv and abs(cu - cv) == d - 1:
        y = -float(ru)
        return [((-stub, y), (0.0, y)), ((float(d - 1), y), (float(d - 1) + stub, y))]

    if cu == cv and abs(ru - rv) == d - 1:
        x = float(cu)
        top_y = 0.0
        bottom_y = -float(d - 1)
        return [((x, top_y), (x, top_y + stub)), ((x, bottom_y), (x, bottom_y - stub))]

    return [(pos[u], pos[v])]


def draw_edges_torus(
    ax: plt.Axes,
    pos: dict[int, tuple[float, float]],
    edges: list[tuple[int, int]],
    d: int,
    color: str,
    width: float,
    alpha: float = 1.0,
    zorder: int = 1,
) -> None:
    segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
    for u, v in edges:
        segments.extend(torus_edge_segments(int(u), int(v), pos, d))
    if not segments:
        return
    collection = LineCollection(segments, colors=color, linewidths=width, alpha=alpha, zorder=zorder)
    ax.add_collection(collection)


def edges_from_fault_ids(graph: nx.Graph, fault_ids: set[int]) -> list[tuple[int, int]]:
    marked_edges: list[tuple[int, int]] = []
    for u, v, attrs in graph.edges(data=True):
        edge_fault_ids = attrs.get("fault_ids", set())
        edge_fault_ids_int = {int(fid) for fid in edge_fault_ids}
        if edge_fault_ids_int & fault_ids:
            marked_edges.append((int(u), int(v)))
    return marked_edges


def draw_gate_channel(
    ax: plt.Axes,
    matching: Matching,
    d: int,
    fault_ids: set[int],
    title: str,
    color: str,
) -> None:
    graph = matching.to_networkx()
    pos = make_positions(graph, d)

    draw_edges_torus(ax, pos, list(graph.edges()), d, color="#C7CBD1", width=1.0, zorder=1)

    marked_edges = edges_from_fault_ids(graph, fault_ids)
    if marked_edges:
        draw_edges_torus(ax, pos, marked_edges, d, color=color, width=5.0, zorder=3)

    nx.draw_networkx_nodes(
        graph,
        pos=pos,
        ax=ax,
        node_color="#8AB6D6",
        edgecolors="#2F2F2F",
        linewidths=0.6,
        node_size=220,
    )
    nx.draw_networkx_labels(graph, pos=pos, ax=ax, font_size=7)

    ax.set_title(f"{title}: gate_edges={len(marked_edges)}", fontsize=11)
    ax.legend(
        handles=[plt.Line2D([0], [0], color=color, lw=5, label=title)],
        loc="upper right",
        fontsize=9,
        frameon=True,
    )
    ax.set_axis_off()


def save_correction_gates_figure(correction_matrix: list[list[int]], output_path: Path) -> None:
    d = validate_correction_matrix(correction_matrix)
    h_x = build_x_check_matrix(d)
    h_z = build_z_check_matrix(d)
    matching_x = Matching(h_x)
    matching_z = Matching(h_z)

    n_qubit_rows = 2 * d
    n_qubit_cols = d
    x_fault_ids: set[int] = set()
    z_fault_ids: set[int] = set()
    for i in range(n_qubit_rows):
        for j in range(n_qubit_cols):
            value = int(correction_matrix[i][j])
            qid = i * n_qubit_cols + j
            if value in (1, 3):
                x_fault_ids.add(qid)
            if value in (2, 3):
                z_fault_ids.add(qid)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
    draw_gate_channel(axes[0], matching_x, d, x_fault_ids, "Apply X gates", "#E11D48")
    draw_gate_channel(axes[1], matching_z, d, z_fault_ids, "Apply Z gates", "#2563EB")
    fig.suptitle("Correction gates on decoder graphs", fontsize=13)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> int:
    if len(sys.argv) not in (2, 3):
        print("Usage: python3 plot_correction_gates.py <correction.json> [<output.png>]", file=sys.stderr)
        return 1

    correction_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2]) if len(sys.argv) == 3 else Path("correction_gates.png")

    try:
        with correction_path.open("r", encoding="utf-8") as f:
            correction_matrix = json.load(f)
        if not isinstance(correction_matrix, list):
            raise ValueError("correction JSON must contain a matrix (list of rows)")

        save_correction_gates_figure(correction_matrix, output_path)
    except Exception as exc:
        print(f"Failed to create correction-gates image: {exc}", file=sys.stderr)
        return 2

    print(f"Saved correction gates image to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
