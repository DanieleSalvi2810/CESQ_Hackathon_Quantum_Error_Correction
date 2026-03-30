#!/usr/bin/env python3
"""
Standalone web demo for toric-code-like decoding with PyMatching.

This file is intentionally self-contained and does not modify any existing
project files.

Features:
- Mobile-friendly web UI
- Click/tap qubits to cycle I -> X -> Z -> Y
- Compute syndromes
- Show decoder pairing paths with animation
- Show correction gates on the original qubit graph
- Final message with correction status
- Optional QR image endpoint (if qrcode[pil] is installed)
"""

from __future__ import annotations

import argparse
import io
import json
import os
import socket
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, quote, urlparse

os.environ.setdefault("MPLCONFIGDIR", "/tmp")

try:
    import numpy as np
    from pymatching import Matching
except ModuleNotFoundError as exc:
    missing = exc.name or "dependency"
    raise SystemExit(
        f"Missing dependency: {missing}. "
        "Install with: python3 -m pip install --user numpy pymatching"
    )

try:
    import qrcode
except ModuleNotFoundError:
    qrcode = None  # Optional dependency for /qr.png endpoint.


UI_VERSION = "1.2.10"


def wrap_index(index: int, size: int) -> int:
    return index % size


def flatten_binary_matrix(matrix: list[list[int]]) -> np.ndarray:
    return np.array([int(value) & 1 for row in matrix for value in row], dtype=np.uint8)


def validate_error_matrix(raw_matrix: object) -> tuple[list[list[int]], int]:
    if not isinstance(raw_matrix, list) or not raw_matrix:
        raise ValueError("error_matrix must be a non-empty 2D list")

    matrix: list[list[int]] = []
    for row in raw_matrix:
        if not isinstance(row, list):
            raise ValueError("error_matrix must be a 2D list")
        matrix.append(row)

    n_rows = len(matrix)
    n_cols = len(matrix[0]) if matrix[0] else 0
    if n_cols == 0:
        raise ValueError("error_matrix has empty rows")
    if any(len(row) != n_cols for row in matrix):
        raise ValueError("error_matrix is ragged")

    if n_rows != 2 * n_cols:
        raise ValueError("expected error_matrix shape (2D x D)")

    normalized: list[list[int]] = []
    for i, row in enumerate(matrix):
        out_row: list[int] = []
        for j, value in enumerate(row):
            if not isinstance(value, int):
                raise ValueError(f"error_matrix[{i}][{j}] must be an integer in [0, 3]")
            if value < 0 or value > 3:
                raise ValueError(f"error_matrix[{i}][{j}] must be in [0, 3]")
            out_row.append(int(value))
        normalized.append(out_row)

    d = n_cols
    return normalized, d


def zeros(rows: int, cols: int) -> list[list[int]]:
    return [[0 for _ in range(cols)] for _ in range(rows)]


def toggle(matrix: list[list[int]], row: int, col: int) -> None:
    n_rows = len(matrix)
    n_cols = len(matrix[0])
    r = wrap_index(row, n_rows)
    c = wrap_index(col, n_cols)
    matrix[r][c] ^= 1


def generate_syndromes(error_matrix: list[list[int]]) -> tuple[list[list[int]], list[list[int]]]:
    d = len(error_matrix[0])
    n_qubit_rows = len(error_matrix)
    n_qubit_cols = len(error_matrix[0])

    syndrome_plaquette = zeros(d, d)
    syndrome_cross = zeros(d, d)

    for i in range(n_qubit_rows):
        for j in range(n_qubit_cols):
            value = error_matrix[i][j]
            even_row = (i % 2 == 0)

            if even_row:
                if value in (1, 3):
                    toggle(syndrome_plaquette, i // 2, j)
                    toggle(syndrome_plaquette, i // 2, j + 1)
                if value in (2, 3):
                    toggle(syndrome_cross, i // 2 - 1, j)
                    toggle(syndrome_cross, i // 2, j)
            else:
                if value in (1, 3):
                    toggle(syndrome_plaquette, (i - 1) // 2, j)
                    toggle(syndrome_plaquette, (i + 1) // 2, j)
                if value in (2, 3):
                    toggle(syndrome_cross, (i - 1) // 2, j - 1)
                    toggle(syndrome_cross, (i - 1) // 2, j)

    return syndrome_plaquette, syndrome_cross


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


def channel_bits_from_error_matrix(error_matrix: list[list[int]]) -> tuple[np.ndarray, np.ndarray]:
    x_bits: list[int] = []
    z_bits: list[int] = []
    for row in error_matrix:
        for value in row:
            x_bits.append(1 if value in (1, 3) else 0)
            z_bits.append(1 if value in (2, 3) else 0)
    return np.array(x_bits, dtype=np.uint8), np.array(z_bits, dtype=np.uint8)


def vector_to_matrix(vector: np.ndarray, d: int) -> list[list[int]]:
    n_qubit_rows = 2 * d
    n_qubit_cols = d
    out = zeros(n_qubit_rows, n_qubit_cols)
    idx = 0
    for i in range(n_qubit_rows):
        for j in range(n_qubit_cols):
            out[i][j] = int(vector[idx]) & 1
            idx += 1
    return out


def combine_corrections(correction_x: np.ndarray, correction_z: np.ndarray, d: int) -> list[list[int]]:
    n_qubit_rows = 2 * d
    n_qubit_cols = d
    out = zeros(n_qubit_rows, n_qubit_cols)
    idx = 0
    for i in range(n_qubit_rows):
        for j in range(n_qubit_cols):
            out[i][j] = int(correction_x[idx]) + 2 * int(correction_z[idx])
            idx += 1
    return out


def xor_matrices(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    rows = len(a)
    cols = len(a[0])
    out = zeros(rows, cols)
    for i in range(rows):
        for j in range(cols):
            out[i][j] = int(a[i][j]) ^ int(b[i][j])
    return out


def flatten_matrix_as_int_list(matrix: list[list[int]]) -> list[int]:
    return [int(v) for row in matrix for v in row]


def numpy_pairs_to_list(arr: np.ndarray) -> list[list[int]]:
    if arr.size == 0:
        return []
    out: list[list[int]] = []
    for item in arr.tolist():
        if not isinstance(item, list):
            continue
        if len(item) != 2:
            continue
        out.append([int(item[0]), int(item[1])])
    return out


def detector_positions(d: int, edges: list[list[int]], pairs: list[list[int]]) -> dict[str, list[float]]:
    nodes: set[int] = set(range(d * d))
    for u, v in edges:
        nodes.add(int(u))
        nodes.add(int(v))
    for u, v in pairs:
        nodes.add(int(u))
        nodes.add(int(v))

    pos: dict[str, list[float]] = {}
    overflow_y = 0.0
    for node in sorted(nodes):
        if 0 <= node < d * d:
            row = node // d
            col = node % d
            pos[str(node)] = [float(col), float(row)]
        else:
            pos[str(node)] = [float(d + 1), float(overflow_y)]
            overflow_y += 1.0
    return pos


def nonzero_count(matrix: list[list[int]]) -> int:
    count = 0
    for row in matrix:
        for value in row:
            if int(value) != 0:
                count += 1
    return count


def decode_request(error_matrix: list[list[int]], horizontal_weight: float, vertical_weight: float) -> dict[str, object]:
    if horizontal_weight <= 0 or vertical_weight <= 0:
        raise ValueError("horizontal_weight and vertical_weight must be > 0")

    validated_error_matrix, d = validate_error_matrix(error_matrix)
    syndrome_plaquette, syndrome_cross = generate_syndromes(validated_error_matrix)

    syndrome_x = flatten_binary_matrix(syndrome_plaquette)
    syndrome_z = flatten_binary_matrix(syndrome_cross)

    h_x = build_x_check_matrix(d)
    h_z = build_z_check_matrix(d)
    weights = make_orientation_weights(d, horizontal_weight, vertical_weight)

    matching_x = Matching(h_x, weights=weights)
    matching_z = Matching(h_z, weights=weights)

    correction_x = np.array(matching_x.decode(syndrome_x), dtype=np.uint8) & 1
    correction_z = np.array(matching_z.decode(syndrome_z), dtype=np.uint8) & 1
    correction_matrix = combine_corrections(correction_x, correction_z, d)

    matched_pairs_x = numpy_pairs_to_list(matching_x.decode_to_matched_dets_array(syndrome_x))
    matched_edges_x = numpy_pairs_to_list(matching_x.decode_to_edges_array(syndrome_x))
    matched_pairs_z = numpy_pairs_to_list(matching_z.decode_to_matched_dets_array(syndrome_z))
    matched_edges_z = numpy_pairs_to_list(matching_z.decode_to_edges_array(syndrome_z))

    residual_matrix = xor_matrices(validated_error_matrix, correction_matrix)
    residual_syndrome_plaquette, residual_syndrome_cross = generate_syndromes(residual_matrix)

    stabilizers_cleared = all(v == 0 for v in flatten_matrix_as_int_list(residual_syndrome_plaquette)) and all(
        v == 0 for v in flatten_matrix_as_int_list(residual_syndrome_cross)
    )
    all_qubits_identity_after_correction = nonzero_count(residual_matrix) == 0

    if stabilizers_cleared and all_qubits_identity_after_correction:
        status_message = "All physical errors are canceled out."
        status_level = "success"
    elif stabilizers_cleared:
        status_message = (
            "Syndrome cleared: stabilizer checks are satisfied, "
            "but a non-trivial lattice configuration remains (possible logical/stabilizer cycle)."
        )
        status_level = "warning"
    else:
        status_message = "Incomplete correction: active syndromes remain."
        status_level = "error"

    return {
        "d": d,
        "error_matrix": validated_error_matrix,
        "weights": {
            "horizontal_weight": horizontal_weight,
            "vertical_weight": vertical_weight,
        },
        "syndromes": {
            "plaquette": syndrome_plaquette,
            "cross": syndrome_cross,
            "flat_x": syndrome_x.astype(int).tolist(),
            "flat_z": syndrome_z.astype(int).tolist(),
        },
        "pairings": {
            "x": {
                "pairs": matched_pairs_x,
                "edges": matched_edges_x,
                "positions": detector_positions(d, matched_edges_x, matched_pairs_x),
            },
            "z": {
                "pairs": matched_pairs_z,
                "edges": matched_edges_z,
                "positions": detector_positions(d, matched_edges_z, matched_pairs_z),
            },
        },
        "correction_matrix": correction_matrix,
        "residual_matrix": residual_matrix,
        "residual_syndromes": {
            "plaquette": residual_syndrome_plaquette,
            "cross": residual_syndrome_cross,
        },
        "stabilizers_cleared": stabilizers_cleared,
        "all_qubits_identity_after_correction": all_qubits_identity_after_correction,
        "status": {
            "level": status_level,
            "message": status_message,
        },
    }


def make_zero_error_matrix(d: int) -> list[list[int]]:
    return zeros(2 * d, d)


def build_page(
  default_d: int,
  default_h_weight: float,
  default_v_weight: float,
  default_share_url: str,
  ui_version: str,
) -> str:
    html = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Quantum Error Correction Live Demo</title>
  <style>
    :root {
      --bg-1: #0e1a2b;
      --bg-2: #162b42;
      --card: rgba(255, 255, 255, 0.94);
      --ink: #0f1723;
      --muted: #5f6b7a;
      --accent-1: #e4572e;
      --accent-2: #2a9d8f;
      --accent-3: #3a86ff;
      --accent-4: #ffb703;
      --error: #d62839;
      --warning: #f4a261;
      --success: #1b9e77;
      --grid: #d7dde4;
    }

    * {
      box-sizing: border-box;
    }

    body {
      margin: 0;
      color: var(--ink);
      background:
        radial-gradient(circle at 10% 10%, rgba(255, 183, 3, 0.22), transparent 33%),
        radial-gradient(circle at 90% 80%, rgba(42, 157, 143, 0.25), transparent 36%),
        linear-gradient(140deg, var(--bg-1), var(--bg-2) 60%, #1d3557);
      font-family: "Space Grotesk", "Avenir Next", "Segoe UI", sans-serif;
      min-height: 100vh;
    }

    .wrap {
      max-width: 1280px;
      margin: 0 auto;
      padding: 18px 14px 28px;
      animation: rise .5s ease-out both;
    }

    @keyframes rise {
      from {
        opacity: 0;
        transform: translateY(10px);
      }
      to {
        opacity: 1;
        transform: translateY(0);
      }
    }

    .hero {
      background: linear-gradient(120deg, rgba(255,255,255,.16), rgba(255,255,255,.07));
      border: 1px solid rgba(255,255,255,.25);
      border-radius: 18px;
      color: #f9fcff;
      padding: 16px 18px;
      margin-bottom: 14px;
      backdrop-filter: blur(3px);
    }

    .hero h1 {
      margin: 0;
      letter-spacing: .2px;
      font-size: clamp(1.3rem, 2.6vw, 2rem);
    }

    .hero p {
      margin: 8px 0 0;
      opacity: 0.95;
    }

    .grid {
      display: grid;
      gap: 12px;
      grid-template-columns: 1.4fr 1fr;
    }

    .full {
      grid-column: 1 / -1;
    }

    .card {
      background: var(--card);
      border: 1px solid rgba(15, 23, 35, .08);
      border-radius: 16px;
      padding: 12px;
      box-shadow: 0 8px 30px rgba(0, 0, 0, .16);
      overflow: hidden;
    }

    h2 {
      margin: 0 0 8px;
      font-size: 1.03rem;
      letter-spacing: .2px;
    }

    .title-with-version {
      display: inline-flex;
      align-items: center;
      gap: 8px;
      flex-wrap: wrap;
    }

    .version-tag {
      display: inline-block;
      padding: 2px 8px;
      border-radius: 999px;
      font-size: .72rem;
      font-weight: 700;
      letter-spacing: .2px;
      background: #eaf2ff;
      color: #1d4ed8;
      border: 1px solid #bcd3ff;
      line-height: 1.4;
    }

    .controls {
      display: grid;
      gap: 8px;
      grid-template-columns: repeat(6, minmax(0, 1fr));
      margin-bottom: 10px;
    }

    .controls label {
      display: block;
      font-size: .82rem;
      color: var(--muted);
      margin-bottom: 4px;
    }

    .controls input {
      width: 100%;
      border: 1px solid #cfd8e3;
      border-radius: 10px;
      padding: 8px 10px;
      font-size: .95rem;
      background: #fbfdff;
    }

    .btn-row {
      grid-column: 1 / -1;
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
    }

    button {
      border: 0;
      border-radius: 10px;
      padding: 9px 12px;
      font-weight: 600;
      color: white;
      background: linear-gradient(135deg, #2563eb, #1d4ed8);
      cursor: pointer;
      transition: transform .15s ease, filter .15s ease;
    }

    button.secondary {
      background: linear-gradient(135deg, #0f766e, #115e59);
    }

    button.ghost {
      color: #123;
      background: #edf2f7;
    }

    button:active {
      transform: translateY(1px);
    }

    button:hover {
      filter: brightness(1.06);
    }

    .hint {
      margin: 6px 0 0;
      font-size: .82rem;
      color: var(--muted);
    }

    .status {
      margin-top: 9px;
      border-radius: 10px;
      padding: 10px 11px;
      font-weight: 600;
      border: 1px solid transparent;
      transition: all .2s ease;
      min-height: 42px;
      display: flex;
      align-items: center;
    }

    .status.success {
      color: #0d5a45;
      background: #e6f8f2;
      border-color: #9fe3cc;
    }

    .status.warning {
      color: #7a4308;
      background: #fff5ea;
      border-color: #f8cf9c;
    }

    .status.error {
      color: #7e0f1f;
      background: #ffecef;
      border-color: #f5a8b4;
    }

    .status.neutral {
      color: #2a3442;
      background: #f4f7fb;
      border-color: #d7e2ef;
    }

    .legend {
      display: flex;
      flex-wrap: wrap;
      gap: 7px 10px;
      margin-top: 8px;
      font-size: .8rem;
      color: #354254;
    }

    .chip {
      border-radius: 999px;
      padding: 4px 9px;
      background: #eef4fa;
      border: 1px solid #d6e2ef;
      display: inline-flex;
      align-items: center;
      gap: 6px;
    }

    .dot {
      width: 10px;
      height: 10px;
      border-radius: 50%;
      display: inline-block;
    }

    .dot.x { background: #d1495b; }
    .dot.z { background: #2f74d0; }
    .dot.y { background: #f39c12; }
    .dot.c { background: #2a9d8f; }

    svg {
      width: 100%;
      border-radius: 12px;
      background: #f9fbff;
      border: 1px solid #dbe5f1;
    }

    .pair-grid {
      display: grid;
      gap: 10px;
      grid-template-columns: 1fr 1fr;
    }

    .pair-controls {
      display: flex;
      align-items: center;
      gap: 8px;
      margin: 2px 0 8px;
      color: #334155;
      font-size: .86rem;
      font-weight: 600;
    }

    .pair-controls input {
      width: 16px;
      height: 16px;
      margin: 0;
    }

    .pair-caption {
      margin: 0 0 4px;
      color: #334155;
      font-weight: 600;
      font-size: .9rem;
    }

    .syndrome-panels {
      display: grid;
      gap: 10px;
      grid-template-columns: 1fr 1fr;
    }

    .mini-grid {
      display: grid;
      gap: 3px;
      margin-top: 6px;
      justify-content: start;
    }

    .mini-cell {
      width: 15px;
      height: 15px;
      border-radius: 4px;
      border: 1px solid #c8d6e6;
      background: #f5f8fc;
    }

    .mini-cell.active {
      background: #e4572e;
      border-color: #c0421f;
    }

    .mini-cell.active.cross {
      background: #2a9d8f;
      border-color: #1f7f74;
    }

    .share-wrap {
      display: grid;
      gap: 8px;
      grid-template-columns: 1fr auto;
      margin-top: 8px;
      align-items: center;
    }

    .share-wrap input {
      border: 1px solid #cfd8e3;
      border-radius: 10px;
      padding: 8px 10px;
      width: 100%;
      background: #fbfdff;
    }

    .qr-box {
      margin-top: 8px;
      display: flex;
      gap: 10px;
      align-items: center;
      flex-wrap: wrap;
    }

    .qr-box img {
      width: 160px;
      height: 160px;
      border-radius: 10px;
      border: 1px solid #d5dfeb;
      background: white;
      object-fit: contain;
    }

    .mono {
      font-family: "IBM Plex Mono", "Cascadia Mono", "Consolas", monospace;
      font-size: .8rem;
      color: #334155;
      white-space: pre-wrap;
      margin: 0;
      max-width: 100%;
      overflow: auto;
    }

    @media (max-width: 980px) {
      .grid {
        grid-template-columns: 1fr;
      }
      .controls {
        grid-template-columns: repeat(2, minmax(0, 1fr));
      }
      .pair-grid, .syndrome-panels {
        grid-template-columns: 1fr;
      }
    }
  </style>
</head>
<body>
  <div class="wrap">
    <section class="hero">
      <h1>Quantum Error Correction Live Demo</h1>
      <p>Set qubit flips (X, Z, Y), compute syndromes, animate pairing paths, and show corrections.</p>
    </section>

    <div class="grid">
      <section class="card">
        <h2 class="title-with-version">Input on the Qubit Graph <span class="version-tag">v__UI_VERSION__</span></h2>
        <div class="controls">
          <div>
            <label for="distanceInput">Code distance (D)</label>
            <input id="distanceInput" type="number" min="2" max="14" step="1" value="__DEFAULT_D__" />
          </div>
          <div>
            <label for="hWeightInput">Horizontal weight</label>
            <input id="hWeightInput" type="number" min="0.0001" step="0.1" value="__DEFAULT_HW__" />
          </div>
          <div>
            <label for="vWeightInput">Vertical weight</label>
            <input id="vWeightInput" type="number" min="0.0001" step="0.1" value="__DEFAULT_VW__" />
          </div>
          <div>
            <label for="pxInput">Random pX</label>
            <input id="pxInput" type="number" min="0" max="1" step="0.01" value="0.10" />
          </div>
          <div>
            <label for="pzInput">Random pZ</label>
            <input id="pzInput" type="number" min="0" max="1" step="0.01" value="0.10" />
          </div>
          <div>
            <label for="seedInput">Random seed</label>
            <input id="seedInput" type="number" step="1" value="" placeholder="empty = random" />
          </div>
          <div class="btn-row">
            <button id="randomBtn" class="ghost">Randomize Errors</button>
            <button id="clearBtn" class="ghost">Clear</button>
            <button id="decodeBtn">Decode and Animate</button>
            <button id="viewResidualBtn" class="ghost">View Corrected State</button>
          </div>
        </div>
        <svg id="qubitSvg"></svg>
        <p class="hint">Tap/click each edge qubit: I -> X -> Z -> Y (implicit toric wrap-around, flat view). Dashed lines show suggested corrections.</p>
        <div class="legend">
          <span class="chip"><span class="dot x"></span>X</span>
          <span class="chip"><span class="dot z"></span>Z</span>
          <span class="chip"><span class="dot y"></span>Y (X+Z)</span>
          <span class="chip"><span class="dot c"></span>Proposed correction X</span>
          <span class="chip"><span class="dot" style="background:#2563eb"></span>Proposed correction Z</span>
        </div>
        <div id="statusBox" class="status neutral">Waiting for decoding.</div>
      </section>

      <section class="card">
        <h2>QR and Sharing</h2>
        <div class="share-wrap">
          <input id="shareUrlInput" />
          <button id="copyBtn" class="secondary">Copy URL</button>
        </div>
        <div class="qr-box">
          <img id="qrImage" alt="QR demo" />
          <div>
            <button id="qrBtn" class="ghost">Generate QR</button>
            <p class="hint" id="qrHint">If /qr.png is not available, install qrcode[pil].</p>
          </div>
        </div>
        <p class="hint">For phone demos, use a reachable URL (LAN or public HTTPS tunnel).</p>
        <pre id="summary" class="mono"></pre>
      </section>

      <section class="card full">
        <h2>Syndromes</h2>
        <div class="syndrome-panels">
          <div>
            <p class="pair-caption">Plaquette syndrome (X-channel)</p>
            <div id="plaquetteGrid" class="mini-grid"></div>
          </div>
          <div>
            <p class="pair-caption">Cross syndrome (Z-channel)</p>
            <div id="crossGrid" class="mini-grid"></div>
          </div>
        </div>
      </section>

      <section class="card full">
        <h2>Pairing Paths with Animation</h2>
        <label class="pair-controls" for="showOriginalErrorsToggle">
          <input id="showOriginalErrorsToggle" type="checkbox" checked />
          Show original errors
        </label>
        <div class="pair-grid">
          <div>
            <p class="pair-caption">X channel</p>
            <svg id="xPairSvg"></svg>
          </div>
          <div>
            <p class="pair-caption">Z channel</p>
            <svg id="zPairSvg"></svg>
          </div>
        </div>
      </section>
    </div>
  </div>

  <script>
    const NS = "http://www.w3.org/2000/svg";
    const LABELS = ["I", "X", "Z", "Y"];
    const ERROR_COLORS = {
      0: "#94a3b8",
      1: "#d1495b",
      2: "#2f74d0",
      3: "#f39c12",
    };
    const CORR_X_COLOR = "#2a9d8f";
    const CORR_Z_COLOR = "#60a5fa";
    const DEFAULT_SEED_PLACEHOLDER = "empty = random";
    const DEFAULT_SHARE_URL = __DEFAULT_SHARE_URL__;
    const state = {
      d: __DEFAULT_D__,
      matrix: [],
      lastResult: null,
      showResidual: false,
      pairAnimationToken: 0,
      pairAnimationFrameBySvg: {},
      pairAnimationTokenBySvg: {},
    };

    function createZeroMatrix(d) {
      return Array.from({ length: 2 * d }, () => Array.from({ length: d }, () => 0));
    }

    function seededRng(seed) {
      let x = Number(seed) >>> 0;
      return () => {
        x ^= x << 13;
        x ^= x >>> 17;
        x ^= x << 5;
        return ((x >>> 0) / 4294967296);
      };
    }

    function randomizeErrors(d, px, pz, seed) {
      const rng = seededRng(seed);
      const matrix = createZeroMatrix(d);
      for (let i = 0; i < 2 * d; i++) {
        for (let j = 0; j < d; j++) {
          const rx = rng();
          const rz = rng();
          let value = 0;
          if (rx < px && rz < pz) value = 3;
          else if (rx < px) value = 1;
          else if (rz < pz) value = 2;
          matrix[i][j] = value;
        }
      }
      return matrix;
    }

    function updateStatus(level, message) {
      const box = document.getElementById("statusBox");
      box.className = "status " + (level || "neutral");
      box.textContent = message;
    }

    function el(name, attrs = {}) {
      const node = document.createElementNS(NS, name);
      Object.entries(attrs).forEach(([k, v]) => node.setAttribute(k, String(v)));
      return node;
    }

    function drawQubitGraph() {
      const svg = document.getElementById("qubitSvg");
      const d = state.d;
      const inputMatrix = state.matrix;
      const matrix = state.showResidual && state.lastResult ? state.lastResult.residual_matrix : inputMatrix;
      const correction = state.lastResult && !state.showResidual ? state.lastResult.correction_matrix : null;
      const margin = 28;
      const cell = Math.max(32, Math.min(56, 430 / d));
      const width = margin * 2 + d * cell;
      const height = margin * 2 + d * cell;
      
      // Extend viewBox to show wrap-around segments on all sides.
      const segmentExtraSpace = cell * 0.8;
      const viewBoxMinX = Math.floor(-segmentExtraSpace);
      const viewBoxMinY = Math.floor(-segmentExtraSpace);
      const viewBoxMaxX = Math.ceil(width + segmentExtraSpace);
      const viewBoxMaxY = Math.ceil(height + segmentExtraSpace);
      const viewBoxWidth = viewBoxMaxX - viewBoxMinX;
      const viewBoxHeight = viewBoxMaxY - viewBoxMinY;
      svg.setAttribute("viewBox", `${viewBoxMinX} ${viewBoxMinY} ${viewBoxWidth} ${viewBoxHeight}`);
      svg.innerHTML = "";

      const bg = el("rect", { x: viewBoxMinX, y: viewBoxMinY, width: viewBoxWidth, height: viewBoxHeight, fill: "#f9fbff" });
      svg.appendChild(bg);

      // Draw only torus vertices (d x d). Using d+1 created an isolated corner dot.
      for (let r = 0; r < d; r++) {
        for (let c = 0; c < d; c++) {
          svg.appendChild(el("circle", {
            cx: margin + c * cell,
            cy: margin + r * cell,
            r: 2.5,
            fill: "#64748b",
          }));
        }
      }

      // Collect all edges to render in two passes
      const edges = [];

      const placeEdge = (i, j, x1, y1, x2, y2) => {
        edges.push({ i, j, x1, y1, x2, y2 });
      };

      // Collect all edges
      for (let i = 0; i < 2 * d; i++) {
        for (let j = 0; j < d; j++) {
          if (i % 2 === 0) {
            const row = i / 2;
            const x1 = margin + j * cell;
            const y = margin + row * cell;
            const x2 = margin + (j + 1) * cell;
            placeEdge(i, j, x1, y, x2, y);
            // Horizontal wrap: draw ghost edge on left side
            if (j === d - 1) {
              const ghostX1 = x1 - d * cell;
              const ghostX2 = x2 - d * cell;
              placeEdge(i, j, ghostX1, y, ghostX2, y);
            }
          } else {
            const row = (i - 1) / 2;
            const x = margin + j * cell;
            const y1 = margin + row * cell;
            const y2 = margin + (row + 1) * cell;
            placeEdge(i, j, x, y1, x, y2);
            // Vertical wrap: draw ghost edge on top side
            if (row === d - 1) {
              const ghostY1 = y1 - d * cell;
              const ghostY2 = y2 - d * cell;
              placeEdge(i, j, x, ghostY1, x, ghostY2);
            }
          }
        }
      }

      // First pass: draw all grey background edges and hit areas
      for (const { i, j, x1, y1, x2, y2 } of edges) {
        svg.appendChild(el("line", {
          x1, y1, x2, y2,
          stroke: "#c4ceda",
          "stroke-width": 5,
          "stroke-linecap": "round",
        }));

        const hit = el("line", {
          x1, y1, x2, y2,
          stroke: "transparent",
          "stroke-width": 18,
          "stroke-linecap": "round",
          cursor: "pointer",
        });
        hit.addEventListener("click", () => {
          inputMatrix[i][j] = (inputMatrix[i][j] + 1) % 4;
          state.lastResult = null;
          updateStatus("neutral", "Input changed. Press 'Decode and Animate'.");
          drawAll();
        });
        svg.appendChild(hit);
      }

      // Second pass: draw all error overlays on top
      for (const { i, j, x1, y1, x2, y2 } of edges) {
        const value = matrix[i][j];
        if (value !== 0) {
          svg.appendChild(el("line", {
            x1, y1, x2, y2,
            stroke: ERROR_COLORS[value],
            "stroke-width": 7.5,
            "stroke-linecap": "round",
          }));
          const tx = (x1 + x2) / 2;
          const ty = (y1 + y2) / 2;
          const isVerticalEdge = Math.abs(x2 - x1) < Math.abs(y2 - y1);
          const text = el("text", {
            x: isVerticalEdge ? tx + 8 : tx,
            y: isVerticalEdge ? ty + 3 : ty - 7,
            "text-anchor": isVerticalEdge ? "start" : "middle",
            "font-size": 10,
            "font-weight": 700,
            fill: "#0f172a",
          });
          text.textContent = LABELS[value];
          svg.appendChild(text);
        }

        if (correction && correction[i]) {
          const cVal = Number(correction[i][j] || 0);
          const hasXCorrection = (cVal & 1) !== 0;
          const hasZCorrection = (cVal & 2) !== 0;

          if (hasXCorrection || hasZCorrection) {
            let xShift = 0;
            let yShift = 0;
            if (hasXCorrection && hasZCorrection) {
              // Split overlays slightly so both dashed corrections remain visible on the same edge.
              if (Math.abs(x2 - x1) >= Math.abs(y2 - y1)) {
                yShift = 1.6;
              } else {
                xShift = 1.6;
              }
            }

            if (hasXCorrection) {
              svg.appendChild(el("line", {
                x1: x1 - xShift,
                y1: y1 - yShift,
                x2: x2 - xShift,
                y2: y2 - yShift,
                stroke: CORR_X_COLOR,
                "stroke-width": 3.1,
                "stroke-dasharray": "5 4",
                "stroke-linecap": "round",
              }));
            }

            if (hasZCorrection) {
              svg.appendChild(el("line", {
                x1: x1 + xShift,
                y1: y1 + yShift,
                x2: x2 + xShift,
                y2: y2 + yShift,
                stroke: CORR_Z_COLOR,
                "stroke-width": 3.1,
                "stroke-dasharray": "5 4",
                "stroke-linecap": "round",
              }));
            }
          }
        }
      }

      // Third pass: draw wider hitboxes for error edges on top
      for (const { i, j, x1, y1, x2, y2 } of edges) {
        const value = matrix[i][j];
        if (value !== 0) {
          const hit = el("line", {
            x1, y1, x2, y2,
            stroke: "transparent",
            "stroke-width": 26,
            "stroke-linecap": "round",
            cursor: "pointer",
          });
          hit.addEventListener("click", () => {
            inputMatrix[i][j] = (inputMatrix[i][j] + 1) % 4;
            state.lastResult = null;
            state.showResidual = false;
            updateStatus("neutral", "Input changed. Press 'Decode and Animate'.");
            drawAll();
          });
          svg.appendChild(hit);
        }
      }
    }

    function drawBinaryGrid(targetId, matrix, activeClass) {
      const target = document.getElementById(targetId);
      target.innerHTML = "";
      if (!matrix || matrix.length === 0) {
        return;
      }
      const rows = matrix.length;
      const cols = matrix[0].length;
      target.style.gridTemplateColumns = `repeat(${cols}, 15px)`;
      target.style.gridTemplateRows = `repeat(${rows}, 15px)`;

      for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
          const cell = document.createElement("div");
          cell.className = "mini-cell" + (matrix[i][j] ? ` active ${activeClass || ""}` : "");
          target.appendChild(cell);
        }
      }
    }

    function toPosMap(d, positions) {
      const out = new Map();
      Object.entries(positions || {}).forEach(([k, val]) => {
        out.set(Number(k), { x: Number(val[0]), y: Number(val[1]) });
      });
      for (let node = 0; node < d * d; node++) {
        if (!out.has(node)) {
          const row = Math.floor(node / d);
          const col = node % d;
          out.set(node, { x: col, y: row });
        }
      }
      return out;
    }

    function drawPairingSvg(svgId, pairingData, flatSyndrome, color, d, sourceMatrix, channel, showOriginalErrors) {
      const svg = document.getElementById(svgId);
      const margin = 26;
      const cell = Math.max(24, Math.min(46, 380 / d));
      const usePlaquetteCenters = channel === "z";

      if (state.pairAnimationFrameBySvg[svgId]) {
        cancelAnimationFrame(state.pairAnimationFrameBySvg[svgId]);
        state.pairAnimationFrameBySvg[svgId] = 0;
      }
      const thisToken = ++state.pairAnimationToken;
      state.pairAnimationTokenBySvg[svgId] = thisToken;
      
      // Extend viewBox to show wrap-around segments poking out on all sides.
      const segmentExtraSpace = cell * 0.8;
      const viewBoxMinX = Math.floor(margin - segmentExtraSpace);
      const viewBoxMinY = Math.floor(margin - segmentExtraSpace);
      const viewBoxMaxX = Math.ceil(margin + d * cell + segmentExtraSpace);
      const viewBoxMaxY = Math.ceil(margin + d * cell + segmentExtraSpace);
      const viewBoxWidth = viewBoxMaxX - viewBoxMinX;
      const viewBoxHeight = viewBoxMaxY - viewBoxMinY;
      svg.setAttribute("viewBox", `${viewBoxMinX} ${viewBoxMinY} ${viewBoxWidth} ${viewBoxHeight}`);
      svg.innerHTML = "";

      const posMap = toPosMap(d, pairingData ? pairingData.positions : {});
      const toVisualPoint = (p) => {
        if (!usePlaquetteCenters) {
          return p;
        }
        return { x: p.x + 0.5, y: p.y + 0.5 };
      };
      const toPx = (p) => ({
        x: margin + p.x * cell,
        y: margin + p.y * cell,
      });
      const toPairPx = (p) => toPx(toVisualPoint(p));

      const drawReferenceQubitGrid = () => {
        // Match the top qubit-grid style in the pairing panes.
        const gridEdges = [];
        const placeGridEdge = (x1, y1, x2, y2) => {
          gridEdges.push({ x1, y1, x2, y2 });
        };

        for (let i = 0; i < 2 * d; i++) {
          for (let j = 0; j < d; j++) {
            if (i % 2 === 0) {
              const row = i / 2;
              const x1 = margin + j * cell;
              const y = margin + row * cell;
              const x2 = margin + (j + 1) * cell;
              placeGridEdge(x1, y, x2, y);
              if (j === d - 1) {
                placeGridEdge(x1 - d * cell, y, x2 - d * cell, y);
              }
            } else {
              const row = (i - 1) / 2;
              const x = margin + j * cell;
              const y1 = margin + row * cell;
              const y2 = margin + (row + 1) * cell;
              placeGridEdge(x, y1, x, y2);
              if (row === d - 1) {
                placeGridEdge(x, y1 - d * cell, x, y2 - d * cell);
              }
            }
          }
        }

        for (const seg of gridEdges) {
            svg.appendChild(el("line", {
              x1: seg.x1,
              y1: seg.y1,
              x2: seg.x2,
              y2: seg.y2,
              stroke: "#c4ceda",
              "stroke-width": 5,
              "stroke-linecap": "round",
            }));
        }
      };

      const addToricSegments = (pu, pv, uNode, vNode, pushSegment) => {
        const isGridNode = (n) => Number.isInteger(n) && n >= 0 && n < d * d;
        const dx = pv.x - pu.x;
        const dy = pv.y - pu.y;

        // Draw toric wrap-around edges as two short segments that exit and re-enter the grid.
        if (isGridNode(uNode) && isGridNode(vNode) && Math.abs(dx) === d - 1 && dy === 0) {
          const leftGhostX = usePlaquetteCenters ? -1 : -1;
          const rightGhostX = usePlaquetteCenters ? d-0.2 : d;
          if (dx > 0) {
            pushSegment(pu, { x: leftGhostX, y: pu.y });
            pushSegment(pv, { x: rightGhostX, y: pv.y });
          } else {
            pushSegment(pu, { x: rightGhostX, y: pu.y });
            pushSegment(pv, { x: leftGhostX, y: pv.y });
          }
          return true;
        }

        if (isGridNode(uNode) && isGridNode(vNode) && Math.abs(dy) === d - 1 && dx === 0) {
          const topGhostY = usePlaquetteCenters ? -1 : -1;
          const bottomGhostY = usePlaquetteCenters ? d-0.2: d;
          if (dy > 0) {
            pushSegment(pu, { x: pu.x, y: topGhostY });
            pushSegment(pv, { x: pv.x, y: bottomGhostY });
          } else {
            pushSegment(pu, { x: pu.x, y: bottomGhostY });
            pushSegment(pv, { x: pv.x, y: topGhostY });
          }
          return true;
        }

        pushSegment(pu, pv);
        return true;
      };

      drawReferenceQubitGrid();

      for (let r = 0; r < d; r++) {
        for (let c = 0; c < d; c++) {
          const p = toPairPx({ x: c, y: r });
          svg.appendChild(el("circle", {
            cx: p.x,
            cy: p.y,
            r: 2.4,
            fill: "#b8c5d6",
          }));
        }
      }

      if (showOriginalErrors && sourceMatrix && sourceMatrix.length === 2 * d) {
        const includeError = (value) => {
          if (channel === "x") return value === 1 || value === 3;
          return value === 2 || value === 3;
        };
        const overlayColor = channel === "x" ? "#d1495b" : "#2f74d0";

        const addOverlaySegment = (from, to) => {
          const projectPoint = channel === "z" ? toPairPx : toPx;
          const a = projectPoint(from);
          const b = projectPoint(to);
          svg.appendChild(el("line", {
            x1: a.x,
            y1: a.y,
            x2: b.x,
            y2: b.y,
            stroke: overlayColor,
            "stroke-width": 4.4,
            "stroke-linecap": "round",
            "stroke-dasharray": "6 4",
            opacity: 0.36,
          }));
        };

        for (let i = 0; i < 2 * d; i++) {
          for (let j = 0; j < d; j++) {
            const value = Number(sourceMatrix[i][j] || 0);
            if (!includeError(value)) continue;

            if (channel === "z") {
              if (i % 2 === 0) {
                const row = i / 2;
                // Dual-lattice orientation: horizontal physical edge -> vertical center link.
                const pu = { x: j, y: (row - 1 + d) % d };
                const pv = { x: j, y: row };
                const uNode = pu.y * d + pu.x;
                const vNode = pv.y * d + pv.x;
                addToricSegments(pu, pv, uNode, vNode, addOverlaySegment);
              } else {
                const row = (i - 1) / 2;
                // Dual-lattice orientation: vertical physical edge -> horizontal center link.
                const pu = { x: (j - 1 + d) % d, y: row };
                const pv = { x: j, y: row };
                const uNode = pu.y * d + pu.x;
                const vNode = pv.y * d + pv.x;
                addToricSegments(pu, pv, uNode, vNode, addOverlaySegment);
              }
              continue;
            }

            if (i % 2 === 0) {
              const row = i / 2;
              if (j === d - 1) {
                addOverlaySegment({ x: j, y: row }, { x: d, y: row });
                addOverlaySegment({ x: 0, y: row }, { x: -1, y: row });
              } else {
                addOverlaySegment({ x: j, y: row }, { x: j + 1, y: row });
              }
            } else {
              const row = (i - 1) / 2;
              if (row === d - 1) {
                addOverlaySegment({ x: j, y: row }, { x: j, y: d });
                addOverlaySegment({ x: j, y: 0 }, { x: j, y: -1 });
              } else {
                addOverlaySegment({ x: j, y: row }, { x: j, y: row + 1 });
              }
            }
          }
        }
      }

      const edges = (pairingData && pairingData.edges) ? pairingData.edges : [];
      const lines = [];
      edges.forEach((edge, idx) => {
        const [u, v] = edge;
        const pu = posMap.get(Number(u));
        const pv = posMap.get(Number(v));
        if (!pu || !pv) return;
        const uNode = Number(u);
        const vNode = Number(v);

        const addAnimatedLine = (from, to, orderIdx) => {
          const a = toPairPx(from);
          const b = toPairPx(to);
          const len = Math.hypot(b.x - a.x, b.y - a.y);
          const line = el("line", {
            x1: a.x,
            y1: a.y,
            x2: b.x,
            y2: b.y,
            stroke: color,
            "stroke-width": 3.2,
            "stroke-linecap": "round",
            opacity: 1,
            "stroke-dasharray": `${len}`,
            "stroke-dashoffset": `${len}`,
          });
          svg.appendChild(line);
          lines.push({ line, idx: orderIdx, len });
        };

        addToricSegments(pu, pv, uNode, vNode, (from, to) => addAnimatedLine(from, to, idx));
      });

      const activeNodes = new Set();
      (flatSyndrome || []).forEach((bit, idx) => {
        if (bit) activeNodes.add(idx);
      });

      posMap.forEach((p, node) => {
        const q = toPairPx(p);
        const isActive = activeNodes.has(node);
        svg.appendChild(el("circle", {
          cx: q.x,
          cy: q.y,
          r: isActive ? 6.5 : 4,
          fill: isActive ? "#ef476f" : "#ffffff",
          stroke: isActive ? "#a3163f" : "#7b8ca2",
          "stroke-width": isActive ? 1.7 : 1.1,
        }));
      });

      if (lines.length === 0) {
        return;
      }

      const perEdgeDelayMs = 90;
      const drawDurationMs = 520;
      const forwardHoldMs = 1520;
      const maxOrderIdx = lines.reduce((acc, item) => Math.max(acc, item.idx), 0);
      const oneWayMs = maxOrderIdx * perEdgeDelayMs + drawDurationMs;
      const cycleMs = 2 * oneWayMs + forwardHoldMs;
      let loopStartTime = performance.now();

      const animateFrame = (now) => {
        if (state.pairAnimationTokenBySvg[svgId] !== thisToken) {
          return;
        }

        const elapsed = now - loopStartTime;
        const phase = elapsed % cycleMs;
        let mirroredElapsed;
        if (phase <= oneWayMs) {
          mirroredElapsed = phase;
        } else if (phase <= oneWayMs + forwardHoldMs) {
          mirroredElapsed = oneWayMs;
        } else {
          mirroredElapsed = cycleMs - phase;
        }

        for (const item of lines) {
          const t = (mirroredElapsed - item.idx * perEdgeDelayMs) / drawDurationMs;
          if (t <= 0) {
            item.line.setAttribute("stroke-dashoffset", `${item.len}`);
            continue;
          }
          if (t >= 1) {
            item.line.setAttribute("stroke-dashoffset", "0");
            continue;
          }
          item.line.setAttribute("stroke-dashoffset", `${item.len * (1 - t)}`);
        }

        state.pairAnimationFrameBySvg[svgId] = requestAnimationFrame(animateFrame);
      };

      state.pairAnimationFrameBySvg[svgId] = requestAnimationFrame(animateFrame);
    }

    function summarize(result) {
      const pairsX = result?.pairings?.x?.pairs?.length || 0;
      const edgesX = result?.pairings?.x?.edges?.length || 0;
      const pairsZ = result?.pairings?.z?.pairs?.length || 0;
      const edgesZ = result?.pairings?.z?.edges?.length || 0;
      const residualPlaquetteOnes = (result?.residual_syndromes?.plaquette || [])
        .flat()
        .reduce((acc, v) => acc + (Number(v) ? 1 : 0), 0);
      const residualCrossOnes = (result?.residual_syndromes?.cross || [])
        .flat()
        .reduce((acc, v) => acc + (Number(v) ? 1 : 0), 0);
      const s = [
        `D=${result.d}`,
        `pairs_x=${pairsX} edges_x=${edgesX}`,
        `pairs_z=${pairsZ} edges_z=${edgesZ}`,
        `residual_syndrome_plaquette_ones=${residualPlaquetteOnes}`,
        `residual_syndrome_cross_ones=${residualCrossOnes}`,
        `stabilizers_cleared=${result.stabilizers_cleared}`,
        `all_qubits_identity_after_correction=${result.all_qubits_identity_after_correction}`,
      ];
      return s.join("\n");
    }

    function drawAll() {
      drawQubitGraph();
      const showOriginalErrors = document.getElementById("showOriginalErrorsToggle")?.checked === true;

      if (!state.lastResult) {
        drawBinaryGrid("plaquetteGrid", createZeroMatrix(state.d).slice(0, state.d), "");
        drawBinaryGrid("crossGrid", createZeroMatrix(state.d).slice(0, state.d), "cross");
        drawPairingSvg("xPairSvg", { edges: [], positions: {} }, [], "#d1495b", state.d, state.matrix, "x", showOriginalErrors);
        drawPairingSvg("zPairSvg", { edges: [], positions: {} }, [], "#2f74d0", state.d, state.matrix, "z", showOriginalErrors);
        document.getElementById("summary").textContent = "";
        return;
      }

      const res = state.lastResult;
      drawBinaryGrid("plaquetteGrid", res.syndromes.plaquette, "");
      drawBinaryGrid("crossGrid", res.syndromes.cross, "cross");
      drawPairingSvg("xPairSvg", res.pairings.x, res.syndromes.flat_x, "#d1495b", state.d, res.error_matrix, "x", showOriginalErrors);
      drawPairingSvg("zPairSvg", res.pairings.z, res.syndromes.flat_z, "#2f74d0", state.d, res.error_matrix, "z", showOriginalErrors);
      document.getElementById("summary").textContent = summarize(res);
      updateStatus(res.status.level, res.status.message);
    }

    async function runDecode() {
      const hWeight = Number(document.getElementById("hWeightInput").value);
      const vWeight = Number(document.getElementById("vWeightInput").value);
      if (!(hWeight > 0) || !(vWeight > 0)) {
        updateStatus("error", "Weights must be > 0.");
        return;
      }

      updateStatus("neutral", "Decoding in progress...");
      try {
        const response = await fetch("/api/decode", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            error_matrix: state.matrix,
            horizontal_weight: hWeight,
            vertical_weight: vWeight,
          }),
        });
        const data = await response.json();
        if (!response.ok) {
          throw new Error(data.error || "API error");
        }
        state.lastResult = data;
        drawAll();
      } catch (err) {
        updateStatus("error", `Error: ${err.message}`);
      }
    }

    function applyDistance(normalizeInput = false) {
      const distanceInput = document.getElementById("distanceInput");
      const raw = Number(distanceInput.value);
      const minBound = Number(distanceInput.min);
      const maxBound = Number(distanceInput.max);
      const minD = Number.isFinite(minBound) ? minBound : 2;
      const maxD = Number.isFinite(maxBound) ? maxBound : 14;
      const parsed = Number.isFinite(raw) ? Math.floor(raw) : state.d;
      const d = Math.max(minD, Math.min(maxD, parsed));

      if (normalizeInput || (Number.isFinite(raw) && (raw < minD || raw > maxD))) {
        distanceInput.value = String(d);
      }

      state.d = d;
      state.matrix = createZeroMatrix(d);
      state.lastResult = null;
      drawAll();
      updateStatus("neutral", `D set to ${d}.`);
    }

    function clearMatrix() {
      const raw = Number(document.getElementById("distanceInput").value);
      const d = Math.max(2, Math.min(14, Number.isFinite(raw) ? Math.floor(raw) : state.d));
      state.d = d;
      state.matrix = createZeroMatrix(d);
      state.lastResult = null;
      state.showResidual = false;
      drawAll();
      updateStatus("neutral", `Cleared and set D to ${d}.`);
    }

    function randomMatrix() {
      const px = Number(document.getElementById("pxInput").value);
      const pz = Number(document.getElementById("pzInput").value);
      const seedInput = document.getElementById("seedInput");
      const seedRaw = seedInput.value.trim();
      const usedAutoSeed = seedRaw === "";
      let seed;
      if (!(px >= 0 && px <= 1 && pz >= 0 && pz <= 1)) {
        updateStatus("error", "pX and pZ must be in [0,1].");
        return;
      }

      if (seedRaw === "") {
        seed = Math.floor(Math.random() * 4294967296) >>> 0;
      } else {
        const parsedSeed = Number(seedRaw);
        if (!Number.isFinite(parsedSeed)) {
          updateStatus("error", "Seed must be a finite number or left empty.");
          return;
        }
        seed = Math.floor(parsedSeed) >>> 0;
      }

      state.matrix = randomizeErrors(state.d, px, pz, seed);
      if (usedAutoSeed) {
        // Keep the field empty but show the sampled seed as suggestion text.
        seedInput.value = "";
        seedInput.placeholder = String(seed);
        seedInput.dataset.suggestedSeed = String(seed);
      } else {
        seedInput.dataset.suggestedSeed = "";
        seedInput.placeholder = DEFAULT_SEED_PLACEHOLDER;
      }
      state.lastResult = null;
      state.showResidual = false;
      drawAll();
      updateStatus("neutral", `Random errors generated (seed=${seed}).`);
    }

    async function copyUrl() {
      const url = document.getElementById("shareUrlInput").value.trim();
      try {
        await navigator.clipboard.writeText(url);
        updateStatus("neutral", "URL copied to clipboard.");
      } catch {
        updateStatus("warning", "Copy failed: copy the URL field manually.");
      }
    }

    async function generateQr() {
      const url = document.getElementById("shareUrlInput").value.trim();
      if (!url) {
        updateStatus("warning", "Enter a valid URL before generating the QR code.");
        return;
      }
      const qrImage = document.getElementById("qrImage");
      qrImage.src = `/qr.png?url=${encodeURIComponent(url)}&t=${Date.now()}`;
      qrImage.onerror = () => {
        document.getElementById("qrHint").textContent =
          "QR unavailable: install qrcode[pil] or use an external QR generator.";
      };
      qrImage.onload = () => {
        document.getElementById("qrHint").textContent =
          "QR ready: scan it with your phone camera.";
      };
    }

    function init() {
      state.matrix = createZeroMatrix(state.d);
      const shareInput = document.getElementById("shareUrlInput");
      const host = window.location.hostname || "";
      const isLocalhost = host === "localhost" || host === "127.0.0.1" || host === "::1";
      shareInput.value = isLocalhost ? DEFAULT_SHARE_URL : window.location.href;

      const seedInput = document.getElementById("seedInput");
      seedInput.dataset.suggestedSeed = "";
      seedInput.addEventListener("keydown", (ev) => {
        if (ev.key !== "Tab") return;
        if (seedInput.value.trim() !== "") return;
        const suggestedSeed = seedInput.dataset.suggestedSeed || "";
        if (suggestedSeed === "") return;
        // Accept suggested seed with Tab, then continue normal focus navigation.
        seedInput.value = suggestedSeed;
        seedInput.dataset.suggestedSeed = "";
        seedInput.placeholder = DEFAULT_SEED_PLACEHOLDER;
      });

      document.getElementById("clearBtn").addEventListener("click", clearMatrix);
      document.getElementById("randomBtn").addEventListener("click", randomMatrix);
      const distanceInput = document.getElementById("distanceInput");
      distanceInput.addEventListener("input", () => applyDistance(false));
      distanceInput.addEventListener("change", () => applyDistance(true));
      distanceInput.addEventListener("blur", () => applyDistance(true));
      document.getElementById("decodeBtn").addEventListener("click", runDecode);
      document.getElementById("viewResidualBtn").addEventListener("click", () => {
        if (!state.lastResult) {
          updateStatus("warning", "Decode first to view corrected state.");
          return;
        }
        state.showResidual = !state.showResidual;
        const btn = document.getElementById("viewResidualBtn");
        btn.textContent = state.showResidual ? "View with Corrections" : "View Corrected State";
        drawAll();
      });
      document.getElementById("copyBtn").addEventListener("click", copyUrl);
      document.getElementById("qrBtn").addEventListener("click", generateQr);
      document.getElementById("showOriginalErrorsToggle").addEventListener("change", drawAll);

      drawAll();
      generateQr();
    }

    init();
  </script>
</body>
</html>
"""
    html = html.replace("__DEFAULT_D__", str(default_d))
    html = html.replace("__DEFAULT_HW__", str(default_h_weight))
    html = html.replace("__DEFAULT_VW__", str(default_v_weight))
    html = html.replace("__DEFAULT_SHARE_URL__", json.dumps(default_share_url))
    html = html.replace("__UI_VERSION__", ui_version)
    return html


def make_qr_png(url: str) -> bytes:
    if qrcode is None:
        raise RuntimeError("qrcode[pil] is not installed")
    if not url or len(url) > 2048:
        raise ValueError("invalid URL for QR")
    img = qrcode.make(url)
    out = io.BytesIO()
    img.save(out, format="PNG")
    return out.getvalue()


class DemoHandler(BaseHTTPRequestHandler):
    page_html: str = ""
    default_d: int = 5
    default_horizontal_weight: float = 1.0
    default_vertical_weight: float = 1.0
    lan_url: str = "http://127.0.0.1:8000/"

    server_version = "HackathonLiveDemo/1.0"

    def log_message(self, format: str, *args: object) -> None:
        # Keep normal access logs, concise and useful during demo.
        super().log_message(format, *args)

    def _send_json(self, payload: object, status: int = 200) -> None:
        raw = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def _send_text(self, text: str, status: int = 200, content_type: str = "text/plain; charset=utf-8") -> None:
        raw = text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)

        if parsed.path == "/":
            self._send_text(self.page_html, content_type="text/html; charset=utf-8")
            return

        if parsed.path == "/health":
            self._send_json({"ok": True})
            return

        if parsed.path == "/api/template":
            query = parse_qs(parsed.query)
            d_raw = query.get("d", [str(self.default_d)])[0]
            try:
                d = int(d_raw)
            except ValueError:
                self._send_json({"error": "d must be an integer"}, status=400)
                return
            if d < 2 or d > 30:
                self._send_json({"error": "d must be in [2,30]"}, status=400)
                return
            self._send_json({"d": d, "error_matrix": make_zero_error_matrix(d)})
            return

        if parsed.path == "/qr.png":
            query = parse_qs(parsed.query)
            url = query.get("url", [""])[0].strip()
            if not url:
                url = self.lan_url
            try:
                png = make_qr_png(url)
            except RuntimeError as exc:
                self._send_text(str(exc), status=HTTPStatus.NOT_IMPLEMENTED)
                return
            except Exception as exc:
                self._send_text(f"QR generation failed: {exc}", status=400)
                return

            self.send_response(200)
            self.send_header("Content-Type", "image/png")
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(png)))
            self.end_headers()
            self.wfile.write(png)
            return

        self._send_json({"error": "Not found"}, status=404)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path != "/api/decode":
            self._send_json({"error": "Not found"}, status=404)
            return

        content_length = self.headers.get("Content-Length")
        if content_length is None:
            self._send_json({"error": "Missing Content-Length"}, status=411)
            return

        try:
            length = int(content_length)
        except ValueError:
            self._send_json({"error": "Invalid Content-Length"}, status=400)
            return

        if length <= 0 or length > 2_000_000:
            self._send_json({"error": "Request body too large or empty"}, status=413)
            return

        raw = self.rfile.read(length)
        try:
            payload = json.loads(raw.decode("utf-8"))
        except Exception:
            self._send_json({"error": "Invalid JSON body"}, status=400)
            return

        if not isinstance(payload, dict):
            self._send_json({"error": "Body must be a JSON object"}, status=400)
            return

        error_matrix = payload.get("error_matrix")
        horizontal_weight = payload.get("horizontal_weight", self.default_horizontal_weight)
        vertical_weight = payload.get("vertical_weight", self.default_vertical_weight)

        try:
            hw = float(horizontal_weight)
            vw = float(vertical_weight)
        except (TypeError, ValueError):
            self._send_json({"error": "horizontal_weight and vertical_weight must be numeric"}, status=400)
            return

        try:
            decoded = decode_request(error_matrix=error_matrix, horizontal_weight=hw, vertical_weight=vw)
        except ValueError as exc:
            self._send_json({"error": str(exc)}, status=400)
            return
        except Exception as exc:
            self._send_json({"error": f"Decoding failed: {exc}"}, status=500)
            return

        self._send_json(decoded, status=200)


def detect_local_ip() -> str:
    # Best effort: discover the LAN IP for sharing with phones in the same Wi-Fi.
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        sock.connect(("8.8.8.8", 80))
        ip = sock.getsockname()[0]
    except Exception:
        ip = "127.0.0.1"
    finally:
        sock.close()
    return ip


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Standalone hackathon web demo for quantum error correction.")
    parser.add_argument("--host", default="0.0.0.0", help="Host bind address (default: 0.0.0.0)")
    parser.add_argument("--port", type=int, default=8000, help="TCP port (default: 8000)")
    parser.add_argument("--d", type=int, default=5, help="Default code distance for UI (default: 5)")
    parser.add_argument(
        "--horizontal-weight",
        type=float,
        default=1.0,
        help="Default horizontal weight in UI (default: 1.0)",
    )
    parser.add_argument(
        "--vertical-weight",
        type=float,
        default=1.0,
        help="Default vertical weight in UI (default: 1.0)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.port < 1 or args.port > 65535:
        raise SystemExit("Port must be in [1, 65535]")
    if args.d < 2 or args.d > 14:
        raise SystemExit("Default D must be in [2, 14]")
    if args.horizontal_weight <= 0 or args.vertical_weight <= 0:
        raise SystemExit("Default weights must be > 0")

    local_ip = detect_local_ip()
    lan_url = f"http://{local_ip}:{args.port}/"

    page = build_page(
        default_d=args.d,
        default_h_weight=args.horizontal_weight,
        default_v_weight=args.vertical_weight,
        default_share_url=lan_url,
      ui_version=UI_VERSION,
    )

    DemoHandler.page_html = page
    DemoHandler.default_d = args.d
    DemoHandler.default_horizontal_weight = args.horizontal_weight
    DemoHandler.default_vertical_weight = args.vertical_weight
    DemoHandler.lan_url = lan_url

    server = ThreadingHTTPServer((args.host, args.port), DemoHandler)
    qr_default_url = f"http://127.0.0.1:{args.port}/qr.png?url={quote(lan_url, safe='')}"

    print("\nLive demo server ready.")
    print(f"Local:  http://127.0.0.1:{args.port}/")
    print(f"LAN:    {lan_url}")
    print("Use the LAN URL for phones in the same Wi-Fi.")
    if qrcode is None:
        print("QR endpoint disabled (optional dependency missing): python3 -m pip install --user qrcode[pil]")
    else:
        print("QR endpoint enabled.")
        print(f"Open this in your browser to get a LAN QR instantly: {qr_default_url}")
    print("Press Ctrl+C to stop.\n")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
