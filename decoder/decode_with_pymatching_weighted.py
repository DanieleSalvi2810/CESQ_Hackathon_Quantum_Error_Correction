#!/usr/bin/env python3
import json
import math
import sys
from pathlib import Path

try:
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


def decode_with_pymatching_weighted(
    syndrome_plaquette: list[list[int]],
    syndrome_cross: list[list[int]],
    px: float,
    pz: float,
) -> list[list[int]]:
    d = validate_syndromes(syndrome_plaquette, syndrome_cross)
    syndrome_x = flatten_matrix(syndrome_plaquette)
    syndrome_z = flatten_matrix(syndrome_cross)

    h_x = build_x_check_matrix(d)
    h_z = build_z_check_matrix(d)
    # Asymmetric+weighted policy requested:
    # - plaquette channel: horizontal=-log10(px), vertical=-log10(pz)
    # - cross channel:     horizontal=-log10(pz), vertical=-log10(px)
    plaquette_horizontal_weight = -math.log10(px)
    plaquette_vertical_weight = -math.log10(pz)
    cross_horizontal_weight = -math.log10(pz)
    cross_vertical_weight = -math.log10(px)

    weights_plaquette = make_orientation_weights(d, plaquette_horizontal_weight, plaquette_vertical_weight)
    weights_cross = make_orientation_weights(d, cross_horizontal_weight, cross_vertical_weight)

    matching_x = Matching(h_x, weights=weights_plaquette)
    matching_z = Matching(h_z, weights=weights_cross)

    correction_x = np.array(matching_x.decode(syndrome_x), dtype=np.uint8) & 1
    correction_z = np.array(matching_z.decode(syndrome_z), dtype=np.uint8) & 1

    correction_matrix: list[list[int]] = []
    n_qubit_rows = 2 * d
    n_qubit_cols = d
    for i in range(n_qubit_rows):
        row: list[int] = []
        for j in range(n_qubit_cols):
            idx = i * n_qubit_cols + j
            if(i%2==0):
                row.append(int(correction_x[idx]) + 2 * int(correction_z[idx]))
            else:
                row.append(int(correction_z[idx]) + 2 * int(correction_x[idx]))
        correction_matrix.append(row)

    return correction_matrix


def main() -> int:
    if len(sys.argv) != 5:
        print(
            "Usage: python3 decode_with_pymatching_weighted.py "
            "<syndrome.json> <correction.json> <px> <pz>",
            file=sys.stderr,
        )
        return 1

    syndrome_path = Path(sys.argv[1])
    correction_path = Path(sys.argv[2])

    try:
        px = float(sys.argv[3])
        pz = float(sys.argv[4])
    except ValueError:
        print("px and pz must be valid floats.", file=sys.stderr)
        return 1

    if px <= 0 or pz <= 0 or px > 1 or pz > 1:
        print("px and pz must be in (0, 1].", file=sys.stderr)
        return 1

    try:
        with syndrome_path.open("r", encoding="utf-8") as f:
            syndrome_data = json.load(f)

        syndrome_plaquette = syndrome_data["syndrome_plaquette"]
        syndrome_cross = syndrome_data["syndrome_cross"]
        correction_matrix = decode_with_pymatching_weighted(
            syndrome_plaquette, syndrome_cross, px, pz
        )

        with correction_path.open("w", encoding="utf-8") as f:
            json.dump(correction_matrix, f, indent=2)
            f.write("\n")
    except Exception as exc:
        print(f"Decoding failed: {exc}", file=sys.stderr)
        return 2

    print(f"Wrote correction matrix to {correction_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
