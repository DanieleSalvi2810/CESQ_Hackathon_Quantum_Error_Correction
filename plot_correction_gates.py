#!/usr/bin/env python3
import json
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx
except ModuleNotFoundError as exc:
    print(
        f"Missing dependency: {exc.name}. Install with: python3 -m pip install --user matplotlib networkx",
        file=sys.stderr,
    )
    raise SystemExit(3)


def validate_correction_matrix(correction_matrix: list[list[int]]) -> tuple[int, int]:
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

    return n_rows, n_cols


def normalize_edge(u: tuple[int, int], v: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    return (u, v) if u <= v else (v, u)


def qubit_to_lattice_edge(i: int, j: int, d: int) -> tuple[tuple[int, int], tuple[int, int]]:
    if i % 2 == 0:
        r = i // 2
        u = (r % d, j % d)
        v = ((r + 1) % d, j % d)
    else:
        r = (i - 1) // 2
        u = (r % d, j % d)
        v = (r % d, (j + 1) % d)
    return normalize_edge(u, v)


def build_toric_graph(d: int) -> nx.Graph:
    graph = nx.Graph()
    for r in range(d):
        for c in range(d):
            graph.add_node((r, c))
    for i in range(2 * d):
        for j in range(d):
            u, v = qubit_to_lattice_edge(i, j, d)
            graph.add_edge(u, v)
    return graph


def draw_correction_graph(correction_matrix: list[list[int]], output_path: Path) -> None:
    _, d = validate_correction_matrix(correction_matrix)
    graph = build_toric_graph(d)

    pos = {(r, c): (float(c), -float(r)) for r in range(d) for c in range(d)}

    x_edges: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    z_edges: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for i in range(2 * d):
        for j in range(d):
            value = int(correction_matrix[i][j])
            edge = qubit_to_lattice_edge(i, j, d)
            if value in (1, 3):
                x_edges.add(edge)
            if value in (2, 3):
                z_edges.add(edge)

    y_edges = x_edges & z_edges
    x_only_edges = x_edges - y_edges
    z_only_edges = z_edges - y_edges

    fig, ax = plt.subplots(1, 1, figsize=(max(7.5, d * 1.4), max(6.5, d * 1.2)), constrained_layout=True)

    nx.draw_networkx_edges(graph, pos=pos, ax=ax, edge_color="#D1D5DB", width=1.2)

    if x_only_edges:
        nx.draw_networkx_edges(
            graph,
            pos=pos,
            ax=ax,
            edgelist=list(x_only_edges),
            edge_color="#E11D48",
            width=4.0,
        )
    if z_only_edges:
        nx.draw_networkx_edges(
            graph,
            pos=pos,
            ax=ax,
            edgelist=list(z_only_edges),
            edge_color="#2563EB",
            width=4.0,
        )
    if y_edges:
        nx.draw_networkx_edges(
            graph,
            pos=pos,
            ax=ax,
            edgelist=list(y_edges),
            edge_color="#7C3AED",
            width=5.0,
        )

    nx.draw_networkx_nodes(
        graph,
        pos=pos,
        ax=ax,
        node_color="#111827",
        node_size=34,
    )

    ax.set_title(
        f"Correction on toric graph (X={len(x_only_edges)} red, Z={len(z_only_edges)} blue, Y={len(y_edges)} purple)",
        fontsize=11,
    )
    legend_items = [
        plt.Line2D([0], [0], color="#E11D48", lw=4, label="Apply X gate"),
        plt.Line2D([0], [0], color="#2563EB", lw=4, label="Apply Z gate"),
    ]
    if y_edges:
        legend_items.append(plt.Line2D([0], [0], color="#7C3AED", lw=5, label="Apply Y (X+Z)"))
    ax.legend(handles=legend_items, loc="upper right", frameon=True, fontsize=9)

    ax.set_axis_off()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> int:
    if len(sys.argv) not in (2, 3):
        print(
            "Usage: python3 plot_correction_gates.py <correction.json> [<output.png>]",
            file=sys.stderr,
        )
        return 1

    correction_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2]) if len(sys.argv) == 3 else Path("correction_gates.png")

    try:
        with correction_path.open("r", encoding="utf-8") as f:
            correction_matrix = json.load(f)
        if not isinstance(correction_matrix, list):
            raise ValueError("correction JSON must contain a matrix (list of rows)")

        draw_correction_graph(correction_matrix, output_path)
    except Exception as exc:
        print(f"Failed to create correction-gates image: {exc}", file=sys.stderr)
        return 2

    print(f"Saved correction gates image to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
