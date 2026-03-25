#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <climits>

#define N 6
#define N_SYNDROME_ROWS (N / 2)

int error_matrix[N][N];
int syndrome_plaquette[N_SYNDROME_ROWS][N];
int syndrome_cross[N_SYNDROME_ROWS][N];
int correction_matrix[N][N];
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

int wrap_syndrome_row(int index) {
    if (index < 0) {
        return N_SYNDROME_ROWS - 1;
    }
    if (index >= N_SYNDROME_ROWS) {
        return 0;
    }
    return index;
}

void toggle_plaquette(int i, int j) {
    int row = wrap_syndrome_row(i);
    int col = wrap_index(j);
    syndrome_plaquette[row][col] = syndrome_plaquette[row][col] ? 0 : 1;
}

void toggle_cross(int i, int j) {
    int row = wrap_syndrome_row(i);
    int col = wrap_index(j);
    syndrome_cross[row][col] = syndrome_cross[row][col] ? 0 : 1;
}

void generate_syndrome_matrices() {
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
                // x error on horizontal qubit - toggles plaquettes left and right
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
                    toggle_plaquette(i / 2, j);
                    toggle_plaquette(i / 2, j + 1);
                }
                // z error on horizontal qubit - toggles crosses above and below
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_cross(i / 2 - 1, j);
                    toggle_cross(i / 2, j);
                }
            } else {
                // x error on vertical qubit - toggles crosses above and below
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
                    toggle_plaquette((i - 1) / 2, j);
                    toggle_plaquette((i + 1) / 2, j);
                }
                // z error on vertical qubit - toggles plaquettes left and right
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_cross((i - 1) / 2, j - 1);
                    toggle_cross((i - 1) / 2, j);
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

void print_full_matrix(const int matrix[N][N], const char *title) {
    printf("%s\n", title);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

void print_syndrome_matrix(const int matrix[N_SYNDROME_ROWS][N], const char *title) {
    printf("%s\n", title);
    for (int i = 0; i < N_SYNDROME_ROWS; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

template <size_t ROWS, size_t COLS>
void write_json_matrix(std::ofstream &out, const int (&matrix)[ROWS][COLS], int indent_spaces) {
    const std::string indent(indent_spaces, ' ');
    const std::string row_indent(indent_spaces + 2, ' ');

    out << "[\n";
    for (size_t i = 0; i < ROWS; i++) {
        out << row_indent << "[";
        for (size_t j = 0; j < COLS; j++) {
            out << matrix[i][j];
            if (j + 1 < COLS) {
                out << ", ";
            }
        }
        out << "]";
        if (i + 1 < ROWS) {
            out << ",";
        }
        out << "\n";
    }
    out << indent << "]";
}

bool write_syndromes_to_json(const std::string &path) {
    std::ofstream out(path);
    if (!out) {
        return false;
    }

    out << "{\n";
    out << "  \"syndrome_plaquette\": ";
    write_json_matrix(out, syndrome_plaquette, 2);
    out << ",\n";
    out << "  \"syndrome_cross\": ";
    write_json_matrix(out, syndrome_cross, 2);
    out << "\n}\n";
    return true;
}

bool read_correction_matrix_from_json(const std::string &path) {
    std::ifstream in(path);
    if (!in) {
        return false;
    }

    const std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::vector<int> values;
    values.reserve(N * N);

    std::string token;
    for (char ch : content) {
        const unsigned char u = static_cast<unsigned char>(ch);
        if (std::isdigit(u) || ch == '-') {
            token.push_back(ch);
            continue;
        }
        if (!token.empty()) {
            values.push_back(std::stoi(token));
            token.clear();
        }
    }
    if (!token.empty()) {
        values.push_back(std::stoi(token));
    }

    if (values.size() != static_cast<size_t>(N * N)) {
        return false;
    }

    size_t idx = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            correction_matrix[i][j] = values[idx++];
        }
    }
    return true;
}

void apply_correction() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            error_matrix[i][j] ^= correction_matrix[i][j];
        }
    }
}

int run_python_decoder(const std::string &syndrome_path, const std::string &correction_path) {
    const std::string cmd = "MPLCONFIGDIR=/tmp python3 decode_with_pymatching.py " + syndrome_path + " " + correction_path;
    return std::system(cmd.c_str());
}

int main() {
    const std::string syndrome_path = "syndrome.json";
    const std::string correction_path = "correction.json";

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            error_matrix[i][j] = 0;
            correction_matrix[i][j] = 0;
        }
    }

    // Esempio deterministico; puoi sostituire con generate_random_error_matrix(error_matrix, px, pz)
    error_matrix[1][1] = 1;
    error_matrix[3][2] = 2;

    generate_syndrome_matrices();
    print_full_matrix(error_matrix, "Error matrix before correction:");
    print_syndrome_matrix(syndrome_plaquette, "Syndrome plaquette:");
    print_syndrome_matrix(syndrome_cross, "Syndrome cross:");

    if (!write_syndromes_to_json(syndrome_path)) {
        fprintf(stderr, "Failed to write syndrome JSON to %s\n", syndrome_path.c_str());
        return 1;
    }

    if (run_python_decoder(syndrome_path, correction_path) != 0) {
        fprintf(stderr, "Python decoder failed. Check pymatching installation.\n");
        return 1;
    }

    if (!read_correction_matrix_from_json(correction_path)) {
        fprintf(stderr, "Failed to read correction matrix from %s\n", correction_path.c_str());
        return 1;
    }

    print_full_matrix(correction_matrix, "Correction matrix from pymatching:");
    apply_correction();
    print_full_matrix(error_matrix, "Residual matrix after applying correction:");

    generate_syndrome_matrices();
    print_syndrome_matrix(syndrome_plaquette, "Residual syndrome plaquette:");
    print_syndrome_matrix(syndrome_cross, "Residual syndrome cross:");

    return 0;
}
