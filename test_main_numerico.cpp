#include <cassert>
#include <iostream>

#define N 10
#define N_SYNDROME_ROWS (N / 2)

void generate_syndrome_matrices(
    const int error_matrix[N][N],
    int syndrome_plaquette[N_SYNDROME_ROWS][N],
    int syndrome_cross[N_SYNDROME_ROWS][N]
);

int test_syndrome_row_index_from_global_row(int row) {
    return row / 2;
}

void clear_full_matrix(int matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = 0;
        }
    }
}

void clear_syndrome_matrix(int matrix[N_SYNDROME_ROWS][N]) {
    for (int i = 0; i < N_SYNDROME_ROWS; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = 0;
        }
    }
}

void assert_syndrome_matrix_equals(
    const int actual[N_SYNDROME_ROWS][N],
    const int expected[N_SYNDROME_ROWS][N],
    const char* test_name,
    const char* matrix_name
) {
    for (int i = 0; i < N_SYNDROME_ROWS; i++) {
        for (int j = 0; j < N; j++) {
            if (actual[i][j] != expected[i][j]) {
                std::cerr << "Test fallito: " << test_name
                          << " su " << matrix_name
                          << " alla cella (" << i << ", " << j << ")"
                          << ". atteso=" << expected[i][j]
                          << ", ottenuto=" << actual[i][j] << "\n";
                assert(false);
            }
        }
    }
}

void run_single_error_case(
    int i,
    int j,
    int error_value,
    const int expected_positions[][2],
    int expected_count,
    const char* test_name
) {
    int error_matrix[N][N];
    int syndrome_plaquette[N_SYNDROME_ROWS][N];
    int syndrome_cross[N_SYNDROME_ROWS][N];
    int expected_plaquette[N_SYNDROME_ROWS][N];
    int expected_cross[N_SYNDROME_ROWS][N];

    clear_full_matrix(error_matrix);
    clear_syndrome_matrix(syndrome_plaquette);
    clear_syndrome_matrix(syndrome_cross);
    clear_syndrome_matrix(expected_plaquette);
    clear_syndrome_matrix(expected_cross);

    error_matrix[i][j] = error_value;

    for (int k = 0; k < expected_count; k++) {
        const int global_row = expected_positions[k][0];
        const int col = expected_positions[k][1];
        const int row = test_syndrome_row_index_from_global_row(global_row);

        if (global_row % 2 == 0) {
            expected_plaquette[row][col] = expected_plaquette[row][col] ? 0 : 1;
        } else {
            expected_cross[row][col] = expected_cross[row][col] ? 0 : 1;
        }
    }

    generate_syndrome_matrices(error_matrix, syndrome_plaquette, syndrome_cross);
    assert_syndrome_matrix_equals(syndrome_plaquette, expected_plaquette, test_name, "plaquette");
    assert_syndrome_matrix_equals(syndrome_cross, expected_cross, test_name, "cross");
}

void test_toggle_cancellation() {
    int error_matrix[N][N];
    int syndrome_plaquette[N_SYNDROME_ROWS][N];
    int syndrome_cross[N_SYNDROME_ROWS][N];
    int expected_plaquette[N_SYNDROME_ROWS][N];
    int expected_cross[N_SYNDROME_ROWS][N];

    clear_full_matrix(error_matrix);
    clear_syndrome_matrix(syndrome_plaquette);
    clear_syndrome_matrix(syndrome_cross);
    clear_syndrome_matrix(expected_plaquette);
    clear_syndrome_matrix(expected_cross);

    // Due X su riga pari con overlap sulla cella (2,4): deve annullarsi (XOR)
    error_matrix[2][3] = 1;
    error_matrix[2][4] = 1;

    expected_plaquette[1][3] = 1;
    expected_plaquette[1][5] = 1;

    generate_syndrome_matrices(error_matrix, syndrome_plaquette, syndrome_cross);
    assert_syndrome_matrix_equals(
        syndrome_plaquette,
        expected_plaquette,
        "cancellazione_toggle_xor",
        "plaquette"
    );
    assert_syndrome_matrix_equals(
        syndrome_cross,
        expected_cross,
        "cancellazione_toggle_xor",
        "cross"
    );
}

int main() {
    {
        const int expected[][2] = {{2, 3}, {2, 4}};
        run_single_error_case(2, 3, 1, expected, 2, "x_su_riga_pari");
    }

    {
        const int expected[][2] = {{1, 3}, {3, 3}};
        run_single_error_case(2, 3, 2, expected, 2, "z_su_riga_pari");
    }

    {
        const int expected[][2] = {{2, 3}, {4, 3}};
        run_single_error_case(3, 3, 1, expected, 2, "x_su_riga_dispari");
    }

    {
        const int expected[][2] = {{3, 3}, {3, 2}};
        run_single_error_case(3, 3, 2, expected, 2, "z_su_riga_dispari");
    }

    {
        const int expected[][2] = {{2, 3}, {4, 3}, {3, 3}, {3, 2}};
        run_single_error_case(3, 3, 3, expected, 4, "y_su_riga_dispari");
    }

    {
        const int expected[][2] = {{5, 0}, {5, 9}};
        run_single_error_case(5, 0, 2, expected, 2, "wrap_j_meno_uno");
    }

    {
        const int expected[][2] = {{9, 4}, {1, 4}};
        run_single_error_case(0, 4, 2, expected, 2, "wrap_i_meno_uno");
    }

    {
        const int expected[][2] = {{0, 9}, {0, 0}};
        run_single_error_case(0, 9, 1, expected, 2, "wrap_j_piu_uno");
    }

    test_toggle_cancellation();

    std::cout << "Tutti i test sulla generazione della sindrome sono passati.\n";
    return 0;
}
