#include <cstdlib>
#include <stdio.h>
#include <vector>
#include <map>
#include <cmath>
#include <climits>

#define N 6
#define N_SYNDROME_ROWS (N / 2)


int error_matrix[N][N];
int syndrome_plaquette[N_SYNDROME_ROWS][N];
int syndrome_cross[N_SYNDROME_ROWS][N];
const double px = 0.1;
const double pz = 0.1;

int wrap_index(int index) {
    if (index < 0) {
        return N - 1;
    }
    if (index >= N) {
        return 0;
    }
    return index;
}


void toggle_plaquette(int i, int j) {
    int row = wrap_index(i);
    int col = wrap_index(j);
    syndrome_plaquette[row][col] = syndrome_plaquette[row][col] ? 0 : 1;
}

void toggle_cross(int i, int j) {
    int row = wrap_index(i);
    int col = wrap_index(j);
    syndrome_cross[row][col] = syndrome_cross[row][col] ? 0 : 1;
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
                //x error on horizontal qubit - toggles plaquettes left and right
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
<<<<<<< Updated upstream
                    toggle_plaquette(i/2, j);
                    toggle_plaquette(i/2, j+1);
                }
                //z error on horizontal qubit - toggles crosses above and below
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_cross(i/2 -1, j);
                    toggle_cross(i/2, j);
=======
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i/2, j);
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i/2, j+1);
                }
                //z error on horizontal qubit - toggles crosses above and below
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i/2, j);
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i/2 - 1, j);
>>>>>>> Stashed changes
                }
                //odd row
            } else {
                //x error on vertical qubit - toggles crosses above and below
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
<<<<<<< Updated upstream
                    toggle_plaquette((i-1)/2, j);
                    toggle_plaquette((i+1)/2, j);
                }
                //z error on vertical qubit - toggles plaquettes left and right
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_cross((i-1)/2, j-1);
                    toggle_cross((i-1)/2, j);
=======
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i/2, j);
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i/2 + 1, j);
                }
                //z error on vertical qubit - toggles plaquettes left and right
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i/2, j);
                    toggle_syndrome(syndrome_plaquette, syndrome_cross, i/2, j + 1);
>>>>>>> Stashed changes
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

// MWPM Decoder using greedy matching algorithm
void decode_with_mwpm(
    const int syndrome_plaquette[N_SYNDROME_ROWS][N],
    const int syndrome_cross[N_SYNDROME_ROWS][N],
    int correction_matrix[N][N]
) {
    // Step 1: Extract active syndrome positions (defects)
    std::vector<std::pair<int, int>> defects;
<<<<<<< Updated upstream

    for (int i = 0; i < N_SYNDROME_ROWS; i++) {
        for (int j = 0; j < N; j++) {
            if (syndrome_plaquette[i][j] == 1) {
                defects.push_back({2*i, j}); // Even rows for plaquettes
            }
            if (syndrome_cross[i][j] == 1) {
                defects.push_back({2*i + 1, j}); // Odd rows for crosses
            }
=======

    for (int i = 0; i < N_SYNDROME_ROWS; i++) {
        for (int j = 0; j < N; j++) {
            if (syndrome_plaquette[i][j] == 1) {
                defects.push_back({2*i, j}); // Even rows for plaquettes
            }
            if (syndrome_cross[i][j] == 1) {
                defects.push_back({2*i + 1, j}); // Odd rows for crosses
            }
        }
    }

    printf("Found %zu defects\n", defects.size());

    if (defects.size() == 0) {
        printf("No defects to match\n");
        return;
    }

    // Step 2: Greedy matching - pair defects by minimum distance
    std::vector<bool> matched(defects.size(), false);

    for (size_t i = 0; i < defects.size(); i++) {
        if (matched[i]) continue;

        int min_dist = INT_MAX;
        int best_match = -1;

        // Find closest unmatched defect
        for (size_t j = i + 1; j < defects.size(); j++) {
            if (matched[j]) continue;

            // Manhattan distance on torus
            int di = abs(defects[i].first - defects[j].first);
            int dj = abs(defects[i].second - defects[j].second);

            // Wrap around (toric boundary conditions)
            di = std::min(di, N - di);
            dj = std::min(dj, N - dj);

            int dist = di + dj;

            if (dist < min_dist) {
                min_dist = dist;
                best_match = j;
            }
        }

        if (best_match != -1) {
            // Mark both as matched
            matched[i] = true;
            matched[best_match] = true;

            // Draw correction path between defects
            int r1 = defects[i].first;
            int c1 = defects[i].second;
            int r2 = defects[best_match].first;
            int c2 = defects[best_match].second;

            // Simple path: first move in row direction, then column
            int r = r1, c = c1;

            while (r != r2 || c != c2) {
                correction_matrix[r][c] = 1;

                if (r != r2) {
                    if (abs(r2 - r) <= N/2) {
                        r = (r2 > r) ? r + 1 : r - 1;
                    } else {
                        r = (r2 > r) ? r - 1 : r + 1;
                    }
                    r = (r + N) % N;
                } else if (c != c2) {
                    if (abs(c2 - c) <= N/2) {
                        c = (c2 > c) ? c + 1 : c - 1;
                    } else {
                        c = (c2 > c) ? c - 1 : c + 1;
                    }
                    c = (c + N) % N;
                }
            }

            printf("Matched defect at (%d,%d) with (%d,%d), distance=%d\n",
                   r1, c1, r2, c2, min_dist);
        }
    }

    printf("MWPM Correction computed successfully!\n");
}

int main() {
    int error_matrix[N][N];
    int syndrome_plaquette[N_SYNDROME_ROWS][N];
    int syndrome_cross[N_SYNDROME_ROWS][N];
    int correction_matrix[N][N];
    const double px = 0.1;
    const double pz = 0.1;

    // generate_random_error_matrix(error_matrix, px, pz);

    for (int i = 0; i < N ; i++) {
        for (int j = 0; j < N; j++) {
            error_matrix[i][j] = 0;
            correction_matrix[i][j] = 0;
>>>>>>> Stashed changes
        }
    }

    printf("Found %zu defects\n", defects.size());

    if (defects.size() == 0) {
        printf("No defects to match\n");
        return;
    }

    // Step 2: Greedy matching - pair defects by minimum distance
    std::vector<bool> matched(defects.size(), false);

    for (size_t i = 0; i < defects.size(); i++) {
        if (matched[i]) continue;

        int min_dist = INT_MAX;
        int best_match = -1;

        // Find closest unmatched defect
        for (size_t j = i + 1; j < defects.size(); j++) {
            if (matched[j]) continue;

            // Manhattan distance on torus
            int di = abs(defects[i].first - defects[j].first);
            int dj = abs(defects[i].second - defects[j].second);

            // Wrap around (toric boundary conditions)
            di = std::min(di, N - di);
            dj = std::min(dj, N - dj);

            int dist = di + dj;

            if (dist < min_dist) {
                min_dist = dist;
                best_match = j;
            }
        }

        if (best_match != -1) {
            // Mark both as matched
            matched[i] = true;
            matched[best_match] = true;

            // Draw correction path between defects
            int r1 = defects[i].first;
            int c1 = defects[i].second;
            int r2 = defects[best_match].first;
            int c2 = defects[best_match].second;

            // Simple path: first move in row direction, then column
            int r = r1, c = c1;

            while (r != r2 || c != c2) {
                correction_matrix[r][c] = 1;

                if (r != r2) {
                    if (abs(r2 - r) <= N/2) {
                        r = (r2 > r) ? r + 1 : r - 1;
                    } else {
                        r = (r2 > r) ? r - 1 : r + 1;
                    }
                    r = (r + N) % N;
                } else if (c != c2) {
                    if (abs(c2 - c) <= N/2) {
                        c = (c2 > c) ? c + 1 : c - 1;
                    } else {
                        c = (c2 > c) ? c - 1 : c + 1;
                    }
                    c = (c + N) % N;
                }
            }

            printf("Matched defect at (%d,%d) with (%d,%d), distance=%d\n",
                   r1, c1, r2, c2, min_dist);
        }
    }

    printf("MWPM Correction computed successfully!\n");
}

int main() {

    // generate_random_error_matrix(error_matrix, px, pz);

    for (int i = 0; i < N ; i++) {
        for (int j = 0; j < N; j++) {
            error_matrix[i][j] = 0;
            correction_matrix[i][j] = 0;
        }
    }

    error_matrix[1][1] = 1;
    error_matrix[3][2] = 2;

    generate_syndrome_matrices(error_matrix, syndrome_plaquette, syndrome_cross);
    print_full_matrix(error_matrix, "Error matrix:");
    print_syndrome_matrix(syndrome_plaquette, "Syndrome plaquette:");
    print_syndrome_matrix(syndrome_cross, "Syndrome cross:");

    // Decode using MWPM
    decode_with_mwpm(syndrome_plaquette, syndrome_cross, correction_matrix);
    print_full_matrix(correction_matrix, "\nCorrection matrix from MWPM:");

    return 0;
}
