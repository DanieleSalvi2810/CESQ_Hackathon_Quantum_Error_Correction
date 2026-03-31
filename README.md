# CESQ Hackathon - Quantum Error Correction Live Demo

Interactive toric-code-style decoding demo built for live presentations.

The main deliverable of this repository is a single-file web app:
- [hackathon_live_demo.py](hackathon_live_demo.py)

It lets you inject Pauli errors on a lattice, decode with PyMatching, animate pairing paths, and inspect residual syndromes in real time.

## Website: What You Can Show Live

The web UI is designed for both desktop and phone demos:
- Click/tap each qubit edge to cycle errors: `I -> X -> Z -> Y`
- Choose code distance `D` (2 to 14)
- Set anisotropic decoder weights (horizontal vs vertical)
- Generate random error patterns with `pX`, `pZ`, and optional seed
- Run decode and animate X/Z matching paths
- Overlay proposed correction gates on the original graph
- Toggle corrected-state view (residual)
- Generate/copy share URL and QR code for audience phones

The backend serves the page and decode API directly from Python `http.server`, so no Node/web framework is needed.

## Quick Start (Recommended)

### 1. Install dependencies

```bash
python3 -m pip install --user -r requirements.txt
```

### 2. Run the live demo server

```bash
python3 hackathon_live_demo.py
```

### 3. Open in browser

- Local: `http://127.0.0.1:8000/`
- LAN (for phones): printed in terminal at startup

If `qrcode[pil]` is installed, a QR endpoint is automatically enabled.

## Launch Options

`hackathon_live_demo.py` supports:

```bash
python3 hackathon_live_demo.py \
	--host 0.0.0.0 \
	--port 8000 \
	--d 5 \
	--horizontal-weight 1.0 \
	--vertical-weight 1.0
```

Arguments:
- `--host`: bind address (default `0.0.0.0`)
- `--port`: TCP port (default `8000`)
- `--d`: default UI distance, range `[2,14]` (default `5`)
- `--horizontal-weight`: default horizontal matching weight, `> 0`
- `--vertical-weight`: default vertical matching weight, `> 0`

## API Endpoints

Served by [hackathon_live_demo.py](hackathon_live_demo.py):

- `GET /` -> HTML web app
- `GET /health` -> health check JSON
- `GET /api/template?d=<int>` -> zero-initialized error matrix template
- `POST /api/decode` -> runs decoding and returns:
	- computed syndromes
	- pairing edges/pairs for X and Z channels
	- correction matrix
	- residual matrix and residual syndromes
	- status message (`success`, `warning`, `error`)
- `GET /qr.png?url=<encoded-url>` -> QR image (requires `qrcode[pil]`)

## Suggested Demo Flow (3-5 minutes)

1. Start at `D=5`, all qubits identity.
2. Inject manual errors (mix `X`, `Z`, `Y`) with clicks.
3. Press Decode and show animated pairing in both channels.
4. Point out correction overlays on the qubit graph.
5. Toggle corrected-state view and explain residual interpretation:
	 - Stabilizers cleared + identity residual -> full cancellation
	 - Stabilizers cleared + non-identity residual -> possible logical/stabilizer cycle
	 - Remaining syndrome -> incomplete correction
6. Randomize with fixed seed for reproducibility.
7. Share QR with audience for phone interaction.

## Visual Asset

The repo includes a toric loop animation:
- [toric_loops.gif](toric_loops.gif)

You can regenerate it with:

```bash
python3 make_gif.py --dpi 300
```

## Other Scripts (Non-Website)

- [toric_codespace.py](toric_codespace.py)
	- Generates logical-sector basis information for symmetric toric code distance `d`.

- [decoder/decode_with_pymatching.py](decoder/decode_with_pymatching.py)
	- CLI decoder from syndrome JSON to correction JSON.

- [decoder/decode_with_pymatching_weighted.py](decoder/decode_with_pymatching_weighted.py)
	- Weighted/asymmetric CLI decoder variant.

These are useful for offline experiments, but the website experience is centered on [hackathon_live_demo.py](hackathon_live_demo.py).

## Project Layout

```text
.
├── hackathon_live_demo.py        # standalone web UI + API server + decoder integration
├── requirements.txt              # Python dependencies
├── make_gif.py                   # toric animation generator
├── toric_loops.gif               # generated animation asset
├── toric_codespace.py            # toric codespace utility script
└── decoder/                      # decoder utilities and static JSON test data
```

## Troubleshooting

- `Missing dependency: pymatching`:
	- Run `python3 -m pip install --user -r requirements.txt`

- QR endpoint disabled:
	- Install optional package: `python3 -m pip install --user qrcode[pil]`

- Phone cannot reach the demo:
	- Run with `--host 0.0.0.0`
	- Use the printed LAN URL
	- Ensure laptop and phone are on the same Wi-Fi
	- Check local firewall rules

- Port already in use:
	- Start with another port, for example `--port 8010`