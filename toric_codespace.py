import argparse
import itertools

def generate_toric_codespace(d):
    """
    Generates the physical qubits, X-stabilizers, and the computational 
    basis strings for the logical |00> state of a distance-d Toric Code.
    """
    # In a Toric code, data qubits live on the edges of a d x d grid.
    # Total qubits = 2 * d^2 (horizontal and vertical edges)
    n_qubits = 2 * d * d
    
    # Helper to map 2D grid coordinates to a 1D array index
    def q_idx(i, j, orientation):
        idx = (i % d) * d + (j % d)
        return idx if orientation == 'h' else idx + (d * d)

    # 1. Generate X-Stabilizers (Vertices)
    # Each vertex touches 4 edges (2 horizontal, 2 vertical)
    x_stabilizers = []
    for i in range(d):
        for j in range(d):
            stab = [
                q_idx(i, j, 'h'),       # Right
                q_idx(i, j-1, 'h'),     # Left
                q_idx(i, j, 'v'),       # Down
                q_idx(i-1, j, 'v')      # Up
            ]
            x_stabilizers.append(stab)
            
    # The product of all X-stabilizers is the Identity matrix. 
    # To avoid generating duplicates, we drop one dependent stabilizer.
    independent_x_stabs = x_stabilizers[:-1]
    
    def build_sector_strings(logical_loop_indices):
        """Generate computational-basis strings for one logical sector."""
        base_state = [0] * n_qubits
        for q in logical_loop_indices:
            base_state[q] ^= 1

        codespace_strings = []
        for bits in itertools.product([0, 1], repeat=len(independent_x_stabs)):
            state = list(base_state)
            for stab_idx, apply_stab in enumerate(bits):
                if apply_stab:
                    for q in independent_x_stabs[stab_idx]:
                        state[q] ^= 1
            codespace_strings.append("".join(map(str, state)))

        codespace_strings.sort()
        return codespace_strings

    # Two independent non-contractible logical X loops (convention):
    # - qubit 1: horizontal loop on horizontal edges, row 0
    # - qubit 2: vertical loop on vertical edges, column 0
    logical_x1 = [q_idx(0, j, 'h') for j in range(d)]
    logical_x2 = [q_idx(i, 0, 'v') for i in range(d)]

    logical_loops = {
        "00": [],
        "01": logical_x2,
        "10": logical_x1,
        "11": logical_x1 + logical_x2,
    }

    print(f"--- Distance {d} Symmetric Toric Code ---")
    print(f"Physical Qubits: {n_qubits}")
    print(f"Independent X-Stabilizers: {len(independent_x_stabs)}")
    print(f"States in each logical sector superposition: 2^{len(independent_x_stabs)}")
    print("Sign convention: every listed basis term has a + sign in this X-stabilizer construction.")
    print("Convention: |10> = X1_bar|00>, |01> = X2_bar|00>, |11> = X1_bar X2_bar |00>\n")

    sector_strings = {label: build_sector_strings(loop) for label, loop in logical_loops.items()}
    sector_sets = {label: set(strings) for label, strings in sector_strings.items()}

    labels = ("00", "01", "10", "11")
    overlap_rows = []
    all_orthogonal = True
    for a, b in itertools.combinations(labels, 2):
        overlap = len(sector_sets[a].intersection(sector_sets[b]))
        overlap_rows.append((a, b, overlap))
        if overlap != 0:
            all_orthogonal = False

    all_distinct = len({frozenset(sector_sets[label]) for label in labels}) == len(labels)

    # Representative computational-basis vector from each logical sector.
    print("Representative computational-basis vectors (with sign):")
    print(f"|00> rep: +|{sector_strings['00'][0]}⟩")
    print(f"|01> rep: +|{sector_strings['01'][0]}⟩")
    print(f"|10> rep: +|{sector_strings['10'][0]}⟩")
    print(f"|11> rep: +|{sector_strings['11'][0]}⟩\n")

    # Explicit listing for small distances.
    if d <= 3:
        for label in ("00", "01", "10", "11"):
            print(f"Logical |{label}> sector basis strings ({len(sector_strings[label])} states):")
            for s in sector_strings[label]:
                print(f"+|{s}⟩")
            print()
    else:
        print("Distance too large to print full sectors explicitly.")
        print("Run with d <= 3 to print all basis strings for |00>, |01>, |10>, |11>.")

        
    print("Pairwise overlap check (must be 0 for orthogonality):")
    for a, b, overlap in overlap_rows:
        print(f"<{a}|{b}> support overlap count = {overlap}")
    print(f"All four sectors distinct: {all_distinct}")
    print(f"All four sectors orthogonal: {all_orthogonal}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate symmetric toric code information for a given distance.",
    )
    parser.add_argument(
        "-d",
        "--distance",
        type=int,
        default=2,
        help="Code distance (positive integer). Default: 2",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.distance < 1:
        raise SystemExit("Distance must be >= 1")
    generate_toric_codespace(d=args.distance)
    
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())