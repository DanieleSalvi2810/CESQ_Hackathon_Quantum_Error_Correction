#!/usr/bin/env python3
import json
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


def validate_syndromes(syndrome_plaquette: list[list[int]], syndrome_cross: list[list[int]]) -> tuple[int, int]:
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
    if n_cols % 2 != 0 or n_rows != n_cols // 2:
        raise ValueError("expected N_SYNDROME_ROWS = N / 2")

    return n_rows, n_cols


def build_x_check_matrix(n: int, n_syndrome_rows: int) -> np.ndarray:
    n_checks = n_syndrome_rows * n
    n_qubits = n * n
    h_x = np.zeros((n_checks, n_qubits), dtype=np.uint8)

    for i in range(n):
        for j in range(n):
            qubit_index = i * n + j
            if i % 2 == 0:
                toggles = [(i // 2, j), (i // 2, j + 1)]
            else:
                toggles = [((i - 1) // 2, j), ((i + 1) // 2, j)]

            for row, col in toggles:
                row_wrapped = wrap_index(row, n_syndrome_rows)
                col_wrapped = wrap_index(col, n)
                check_index = row_wrapped * n + col_wrapped
                h_x[check_index, qubit_index] ^= 1

    return h_x


def build_z_check_matrix(n: int, n_syndrome_rows: int) -> np.ndarray:
    n_checks = n_syndrome_rows * n
    n_qubits = n * n
    h_z = np.zeros((n_checks, n_qubits), dtype=np.uint8)

    for i in range(n):
        for j in range(n):
            qubit_index = i * n + j
            if i % 2 == 0:
                toggles = [(i // 2 - 1, j), (i // 2, j)]
            else:
                toggles = [((i - 1) // 2, j - 1), ((i - 1) // 2, j)]

            for row, col in toggles:
                row_wrapped = wrap_index(row, n_syndrome_rows)
                col_wrapped = wrap_index(col, n)
                check_index = row_wrapped * n + col_wrapped
                h_z[check_index, qubit_index] ^= 1

    return h_z


def decode_with_pymatching(syndrome_plaquette: list[list[int]], syndrome_cross: list[list[int]]) -> list[list[int]]:
    n_syndrome_rows, n = validate_syndromes(syndrome_plaquette, syndrome_cross)
    syndrome_x = flatten_matrix(syndrome_plaquette)
    syndrome_z = flatten_matrix(syndrome_cross)

    h_x = build_x_check_matrix(n, n_syndrome_rows)
    h_z = build_z_check_matrix(n, n_syndrome_rows)

    matching_x = Matching(h_x)
    matching_z = Matching(h_z)

    correction_x = np.array(matching_x.decode(syndrome_x), dtype=np.uint8) & 1
    correction_z = np.array(matching_z.decode(syndrome_z), dtype=np.uint8) & 1

    correction_matrix: list[list[int]] = []
    for i in range(n):
        row: list[int] = []
        for j in range(n):
            idx = i * n + j
            row.append(int(correction_x[idx]) + 2 * int(correction_z[idx]))
        correction_matrix.append(row)

    return correction_matrix


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: python3 decode_with_pymatching.py <syndrome.json> <correction.json>", file=sys.stderr)
        return 1

    syndrome_path = Path(sys.argv[1])
    correction_path = Path(sys.argv[2])

    try:
        with syndrome_path.open("r", encoding="utf-8") as f:
            syndrome_data = json.load(f)

        syndrome_plaquette = syndrome_data["syndrome_plaquette"]
        syndrome_cross = syndrome_data["syndrome_cross"]
        correction_matrix = decode_with_pymatching(syndrome_plaquette, syndrome_cross)

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
