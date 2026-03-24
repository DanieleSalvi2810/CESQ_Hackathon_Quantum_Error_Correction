#include <cstdlib>
#include <stdio.h>

#define N 10
#define N_SYNDROME_ROWS (N / 2)

int wrap_index(int index) {
    if (index < 0) {
        return N - 1;
    }
    if (index >= N) {
        return 0;
    }
    return index;
}

int syndrome_row_index_from_global_row(int row) {
    return row / 2;
}

int* get_syndrome(
    int syndrome_plaquette[N_SYNDROME_ROWS][N],
    int syndrome_cross[N_SYNDROME_ROWS][N],
    int i,
    int j
) {
    const int global_row = wrap_index(i);
    const int col = wrap_index(j);
    const int row = syndrome_row_index_from_global_row(global_row);

    if (global_row % 2 == 0) {
        return &syndrome_plaquette[row][col];
    }
    return &syndrome_cross[row][col];
}

void toggle_syndrome(
    int syndrome_plaquette[N_SYNDROME_ROWS][N],
    int syndrome_cross[N_SYNDROME_ROWS][N],
    int i,
    int j
) {
    int* syndrome = get_syndrome(syndrome_plaquette, syndrome_cross, i, j);
    *syndrome = *syndrome ? 0 : 1;
}

void generate_syndrome_matrices(
    const int error_matrix[N][N],
    int syndrome_plaquette[N_SYNDROME_ROWS][N],
    int syndrome_cross[N_SYNDROME_ROWS][N]
) {
    for (int i = 0; i < N_SYNDROME_ROWS; i++) {
        for (int j = 0; j < N; j++) {
            syndrome_plaquette[i][j] = 0;
            syndrome_cross[i][j] = 0;
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            const bool even_row = (i % 2 == 0);

            if (even_row) {
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
                    // X error on even rows: self + right neighbour
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i, j);
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i, j + 1);
                }
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    // Z error on even rows: up + down neighbours
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i - 1, j);
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i + 1, j);
                }
            } else {
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
                    // X error on odd rows: up + down neighbours
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i - 1, j);
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i + 1, j);
                }
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    // Z error on odd rows: self + left neighbour
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i, j);
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i, j - 1);
                }
            }
        }
    }
}

void generate_random_error_matrix(int error_matrix[N][N], double px, double pz) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            error_matrix[i][j] = 0;

            const double rx = static_cast<double>(rand()) / RAND_MAX;
            const double rz = static_cast<double>(rand()) / RAND_MAX;

            if (rx < px && rz < pz) {
                error_matrix[i][j] = 3;
            } else if (rx < px) {
                error_matrix[i][j] = 1;
            } else if (rz < pz) {
                error_matrix[i][j] = 2;
            }
        }
    }
}

void print_full_matrix(const int matrix[N][N], const char* title) {
    printf("%s\n", title);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

void print_syndrome_matrix(const int matrix[N_SYNDROME_ROWS][N], const char* title) {
    printf("%s\n", title);
    for (int i = 0; i < N_SYNDROME_ROWS; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {
    int error_matrix[N][N];
    int syndrome_plaquette[N_SYNDROME_ROWS][N];
    int syndrome_cross[N_SYNDROME_ROWS][N];
    const double px = 0.1;
    const double pz = 0.1;

    generate_random_error_matrix(error_matrix, px, pz);
    
    for (int i = 0; i < N ; i++) {
        for (int j = 0; j < N; j++) {
            error_matrix[i][j] = 0;
        }
    }

    error_matrix[3][2] = 1;
    error_matrix[5][4] = 2;

    generate_syndrome_matrices(error_matrix, syndrome_plaquette, syndrome_cross);
    print_full_matrix(error_matrix, "Error matrix:");
    print_syndrome_matrix(syndrome_plaquette, "Syndrome plaquette:");
    print_syndrome_matrix(syndrome_cross, "Syndrome cross:");

    return 0;
}
