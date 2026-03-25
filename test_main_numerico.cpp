#include <cctype>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

struct ErrorEntry {
    int row;
    int col;
    int value;  // 0=I, 1=X, 2=Z, 3=Y(X+Z)
};

enum class ErrorMode {
    Manual,
    Random,
};

struct TestCase {
    std::string name;
    int d;  // code distance
    ErrorMode mode;
    double px;  // usato solo in modalita' Random
    double pz;  // usato solo in modalita' Random
    int seed;   // -1 = seed da tempo corrente
    std::vector<ErrorEntry> errors;
};

// Cambia solo questa variabile per scegliere il test da eseguire.
static constexpr size_t ACTIVE_TEST_CASE = 4;

static const std::vector<TestCase> TEST_CASES = {
    #0
    {
        "base_two_errors_d3",
        3,
        ErrorMode::Manual,
        0.0,
        0.0,
        -1,
        {
            {1, 1, 1},
            {3, 2, 2},
        },
    },
    #1
    {
        "single_y_error_d3",
        3,
        ErrorMode::Manual,
        0.0,
        0.0,
        -1,
        {
            {2, 4, 3},
        },
    },
    #2
    {
        "mixed_errors_d4",
        4,
        ErrorMode::Manual,
        0.0,
        0.0,
        -1,
        {
            {1, 1, 1},
            {2, 6, 2},
            {5, 3, 3},
        },
    },
    #3
    {
        "random_generated_d3",
        3,
        ErrorMode::Random,
        0.10,
        0.10,
        -1,
        {},
    },
    #4
    {
        "random_generated_d4",
        50,
        ErrorMode::Random,
        0.15,
        0.05,
        42,
        {},
    },
};

int wrap_index(int index, int size) {
    int wrapped = index % size;
    if (wrapped < 0) {
        wrapped += size;
    }
    return wrapped;
}

void print_matrix(const std::vector<std::vector<int>> &matrix, const char *title) {
    std::cout << title << "\n";
    for (const auto &row : matrix) {
        for (int value : row) {
            std::cout << value << " ";
        }
        std::cout << "\n";
    }
}

void toggle(std::vector<std::vector<int>> &matrix, int row, int col) {
    const int n_rows = static_cast<int>(matrix.size());
    const int n_cols = static_cast<int>(matrix[0].size());
    const int r = wrap_index(row, n_rows);
    const int c = wrap_index(col, n_cols);
    matrix[r][c] = matrix[r][c] ? 0 : 1;
}

void generate_syndrome_matrices(
    const std::vector<std::vector<int>> &error_matrix,
    std::vector<std::vector<int>> &syndrome_plaquette,
    std::vector<std::vector<int>> &syndrome_cross
) {
    const int n = static_cast<int>(error_matrix.size());

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            const bool even_row = (i % 2 == 0);
            const int value = error_matrix[i][j];

            if (even_row) {
                if (value == 1 || value == 3) {
                    toggle(syndrome_plaquette, i / 2, j);
                    toggle(syndrome_plaquette, i / 2, j + 1);
                }
                if (value == 2 || value == 3) {
                    toggle(syndrome_cross, i / 2 - 1, j);
                    toggle(syndrome_cross, i / 2, j);
                }
            } else {
                if (value == 1 || value == 3) {
                    toggle(syndrome_plaquette, (i - 1) / 2, j);
                    toggle(syndrome_plaquette, (i + 1) / 2, j);
                }
                if (value == 2 || value == 3) {
                    toggle(syndrome_cross, (i - 1) / 2, j - 1);
                    toggle(syndrome_cross, (i - 1) / 2, j);
                }
            }
        }
    }
}

void generate_random_error_matrix(std::vector<std::vector<int>> &error_matrix, double px, double pz) {
    const int n = static_cast<int>(error_matrix.size());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            error_matrix[i][j] = 0;

            const double rx = static_cast<double>(std::rand()) / RAND_MAX;
            const double rz = static_cast<double>(std::rand()) / RAND_MAX;

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

template <typename MatrixT>
void write_json_matrix(std::ofstream &out, const MatrixT &matrix, int indent_spaces) {
    const std::string indent(indent_spaces, ' ');
    const std::string row_indent(indent_spaces + 2, ' ');

    out << "[\n";
    for (size_t i = 0; i < matrix.size(); i++) {
        out << row_indent << "[";
        for (size_t j = 0; j < matrix[i].size(); j++) {
            out << matrix[i][j];
            if (j + 1 < matrix[i].size()) {
                out << ", ";
            }
        }
        out << "]";
        if (i + 1 < matrix.size()) {
            out << ",";
        }
        out << "\n";
    }
    out << indent << "]";
}

bool write_syndromes_to_json(
    const std::string &path,
    const std::vector<std::vector<int>> &syndrome_plaquette,
    const std::vector<std::vector<int>> &syndrome_cross
) {
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

bool read_correction_matrix_from_json(
    const std::string &path,
    std::vector<std::vector<int>> &correction_matrix,
    int n
) {
    std::ifstream in(path);
    if (!in) {
        return false;
    }

    const std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::vector<int> values;
    values.reserve(static_cast<size_t>(n * n));

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

    if (values.size() != static_cast<size_t>(n * n)) {
        return false;
    }

    size_t idx = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            correction_matrix[i][j] = values[idx++];
        }
    }
    return true;
}

void apply_correction(
    std::vector<std::vector<int>> &error_matrix,
    const std::vector<std::vector<int>> &correction_matrix
) {
    const int n = static_cast<int>(error_matrix.size());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            error_matrix[i][j] ^= correction_matrix[i][j];
        }
    }
}

int run_python_decoder(const std::string &syndrome_path, const std::string &correction_path) {
    const std::string cmd = "MPLCONFIGDIR=/tmp python3 decode_with_pymatching.py " + syndrome_path + " " + correction_path;
    return std::system(cmd.c_str());
}

int main() {
    if (TEST_CASES.empty()) {
        std::cerr << "No test cases defined.\n";
        return 1;
    }
    if (ACTIVE_TEST_CASE >= TEST_CASES.size()) {
        std::cerr << "ACTIVE_TEST_CASE out of range. Available: 0.." << (TEST_CASES.size() - 1) << "\n";
        return 1;
    }

    const TestCase &tc = TEST_CASES[ACTIVE_TEST_CASE];
    const int d = tc.d;
    const int n = 2 * d;
    const int n_syndrome_rows = d;

    if (d <= 0) {
        std::cerr << "Invalid test case: d must be > 0.\n";
        return 1;
    }

    std::vector<std::vector<int>> error_matrix(n, std::vector<int>(n, 0));
    std::vector<std::vector<int>> syndrome_plaquette(n_syndrome_rows, std::vector<int>(n, 0));
    std::vector<std::vector<int>> syndrome_cross(n_syndrome_rows, std::vector<int>(n, 0));
    std::vector<std::vector<int>> correction_matrix(n, std::vector<int>(n, 0));

    if (tc.mode == ErrorMode::Random) {
        if (tc.px < 0.0 || tc.px > 1.0 || tc.pz < 0.0 || tc.pz > 1.0) {
            std::cerr << "Invalid random probabilities in test '" << tc.name
                      << "': px and pz must be in [0,1].\n";
            return 1;
        }

        const unsigned int seed_to_use =
            tc.seed >= 0 ? static_cast<unsigned int>(tc.seed) : static_cast<unsigned int>(std::time(nullptr));
        std::srand(seed_to_use);
        generate_random_error_matrix(error_matrix, tc.px, tc.pz);

        std::cout << "Random mode enabled: px=" << tc.px
                  << ", pz=" << tc.pz
                  << ", seed=" << seed_to_use << "\n";
    } else {
        for (const auto &entry : tc.errors) {
            if (entry.row < 0 || entry.row >= n || entry.col < 0 || entry.col >= n) {
                std::cerr << "Invalid error coordinate in test '" << tc.name
                          << "': (" << entry.row << "," << entry.col << ") for N=" << n << "\n";
                return 1;
            }
            if (entry.value < 0 || entry.value > 3) {
                std::cerr << "Invalid error value in test '" << tc.name << "': " << entry.value << "\n";
                return 1;
            }
            error_matrix[entry.row][entry.col] = entry.value;
        }
    }

    std::cout << "Running test case #" << ACTIVE_TEST_CASE << ": " << tc.name
              << " (d=" << d << ", N=" << n << ")\n";

    generate_syndrome_matrices(error_matrix, syndrome_plaquette, syndrome_cross);
    print_matrix(error_matrix, "Error matrix before correction:");
    print_matrix(syndrome_plaquette, "Syndrome plaquette:");
    print_matrix(syndrome_cross, "Syndrome cross:");

    const std::string syndrome_path = "syndrome.json";
    const std::string correction_path = "correction.json";

    if (!write_syndromes_to_json(syndrome_path, syndrome_plaquette, syndrome_cross)) {
        std::cerr << "Failed to write syndrome JSON to " << syndrome_path << "\n";
        return 1;
    }

    if (run_python_decoder(syndrome_path, correction_path) != 0) {
        std::cerr << "Python decoder failed. Check pymatching installation.\n";
        return 1;
    }

    if (!read_correction_matrix_from_json(correction_path, correction_matrix, n)) {
        std::cerr << "Failed to read correction matrix from " << correction_path << "\n";
        return 1;
    }

    print_matrix(correction_matrix, "Correction matrix from pymatching:");
    apply_correction(error_matrix, correction_matrix);
    print_matrix(error_matrix, "Residual matrix after applying correction:");

    std::vector<std::vector<int>> residual_syndrome_plaquette(n_syndrome_rows, std::vector<int>(n, 0));
    std::vector<std::vector<int>> residual_syndrome_cross(n_syndrome_rows, std::vector<int>(n, 0));
    generate_syndrome_matrices(error_matrix, residual_syndrome_plaquette, residual_syndrome_cross);
    print_matrix(residual_syndrome_plaquette, "Residual syndrome plaquette:");
    print_matrix(residual_syndrome_cross, "Residual syndrome cross:");

    return 0;
}
