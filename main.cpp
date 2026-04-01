#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdio.h>
#include <string>
#include <tuple>
#include <vector>
#include <iostream>

#define D 4
#define N (D * 2)

int error_matrix[N][D];
int syndrome_plaquette[D][D];
int syndrome_cross[D][D];
int correction_matrix[N][D];
const double px = 0.1;
const double pz = 0.1;

int wrap_index(int index, int size) {
    int wrapped = index % size;
    if (wrapped < 0) {
        wrapped += size;
    }
    return wrapped;
}

void toggle_plaquette(int i, int j) {
    int row = wrap_index(i, D);
    int col = wrap_index(j, D);
    syndrome_plaquette[row][col] ^= 1;
}

void toggle_cross(int i, int j) {
    int row = wrap_index(i, D);
    int col = wrap_index(j, D);
    syndrome_cross[row][col] ^= 1;
}

void generate_syndrome_matrices() {
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            syndrome_plaquette[i][j] = 0;
            syndrome_cross[i][j] = 0;
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
            const bool even_row = (i % 2 == 0);

            if (even_row) {
                // X error on vertical qubit -> adjacent plaquettes.
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
                    toggle_plaquette(i / 2, j);
                    toggle_plaquette(i / 2, j + 1);
                }
                // Z error on vertical qubit -> adjacent stars/crosses.
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_cross(i / 2 - 1, j);
                    toggle_cross(i / 2, j);
                }
            } else {
                // X error on horizontal qubit -> adjacent plaquettes.
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
                    toggle_plaquette((i - 1) / 2, j);
                    toggle_plaquette((i + 1) / 2, j);
                }
                // Z error on horizontal qubit -> adjacent stars/crosses.
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_cross((i - 1) / 2, j - 1);
                    toggle_cross((i - 1) / 2, j);
                }
            }
        }
    }
}

void generate_syndrome_matrices_asymmetric() {
    for (int i = 0; i < D; i++) {
        for (int j = 0; j < D; j++) {
            syndrome_plaquette[i][j] = 0;
            syndrome_cross[i][j] = 0;
        }
    }
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
            const bool even_row = (i % 2 == 0);

            if (even_row) {
                // X error on vertical qubit -> adjacent plaquettes.
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
                    toggle_plaquette(i / 2, j);
                    toggle_plaquette(i / 2, j + 1);
                }
                // Z error on vertical qubit -> adjacent crosses.
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_cross(i / 2 - 1, j);
                    toggle_cross(i / 2, j);
                }
            } else {
                // X error on horizontal qubit -> adjacent crosses.
                if (error_matrix[i][j] == 1 || error_matrix[i][j] == 3) {
                    toggle_cross((i - 1) / 2, j - 1);
                    toggle_cross((i - 1) / 2, j);
                }
                // Z error on horizontal qubit -> adjacent plaquettes.
                if (error_matrix[i][j] == 2 || error_matrix[i][j] == 3) {
                    toggle_plaquette((i - 1) / 2, j);
                    toggle_plaquette((i + 1) / 2, j);
                }
            }
        }
    }
}

void generate_random_error_matrix(double px, double pz) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
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

template <size_t ROWS, size_t COLS>
void print_matrix(const int (&matrix)[ROWS][COLS], const char *title) {
    printf("%s\n", title);
    for (size_t i = 0; i < ROWS; i++) {
        for (size_t j = 0; j < COLS; j++) {
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
    values.reserve(N * D);

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

    if (values.size() != static_cast<size_t>(N * D)) {
        return false;
    }

    size_t idx = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
            correction_matrix[i][j] = values[idx++];
        }
    }
    return true;
}

void apply_correction() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
            error_matrix[i][j] ^= correction_matrix[i][j];
        }
    }
}

int run_python_decoder(const std::string &syndrome_path, const std::string &correction_path) {
    std::string script_path = "decode_with_pymatching.py";
    if (!std::filesystem::exists(script_path)) {
        script_path = "../decode_with_pymatching.py";
    }

    const std::string cmd =
        "MPLCONFIGDIR=/tmp python3 " + script_path + " " + syndrome_path + " " + correction_path;
    return std::system(cmd.c_str());
}


std::tuple<int, int> same_right(std::tuple<int, int> position){
    int i, j;
    std::tie(i, j) = position;
    return std::make_tuple(i, j+1);
}

std::tuple<int, int> same_left(std::tuple<int, int> position){
    int i, j;
    std::tie(i, j) = position;
    return std::make_tuple(i, j-1);
}

std::tuple<int, int> same_up(std::tuple<int, int> position){
    int i, j;
    std::tie(i, j) = position;
    return std::make_tuple(i-2, j);
}

std::tuple<int, int> same_down(std::tuple<int, int> position){
    int i, j;
    std::tie(i, j) = position;
    return std::make_tuple(i+2, j);
}

std::tuple<int, int> different_down_right(std::tuple<int, int> position){
    int i, j;
    std::tie(i, j) = position;
    if (i % 2 != 0){
        return std::make_tuple(i+1, j);
    }
    return std::make_tuple(i+1, j+1);
}

std::tuple<int, int> different_down_left(std::tuple<int, int> position){
    int i, j;
    std::tie(i, j) = position;
    if (i % 2 != 0){
        return std::make_tuple(i+1, j-1);
    }
    return std::make_tuple(i+1, j);
}

std::tuple<int, int> different_up_left(std::tuple<int, int> position){
    int i, j;
    std::tie(i, j) = position;
    if (i % 2 != 0){
        return std::make_tuple(i-1, j-1);
    }
    return std::make_tuple(i-1, j);
}

std::tuple<int, int> different_up_right(std::tuple<int, int> position){
    int i, j;
    std::tie(i, j) = position;
    if (i % 2 != 0){
        return std::make_tuple(i-1, j);
    }
    return std::make_tuple(i-1, j+1);
}


void insert_in_path(std::tuple<int, int> position, std::vector<std::tuple<int, int>> &path){
    //insert position in path, update error_matrix and syndrome matrices accordingly
    path.push_back(position);
}

void remove_errors(std::vector<std::tuple<int, int>> &path){
    for (int i = 0; i < path.size(); i++){
        //remove error in position path[i], update error_matrix and syndrome matrices accordingly
        const auto &[row, col] = path[i];
        error_matrix[row][col] = 0;

    }
}


int check_visited(std::tuple<int, int> position, std::vector<std::tuple<int, int>> &path){
    for (int i = path.size() - 1; i >= 0; i--){
        if (std::get<0>(path[i]) == std::get<0>(position) && std::get<1>(path[i]) == std::get<1>(position)){
            return i;
        }
    }
    return -1;
}



std::tuple<int, int> choose_next(std::tuple<int, int> position, int error_code, std::string &momentum, std::vector<std::tuple<int, int>> &path){
    //choose next position in path, according to the error code and the current position
    //return (-1, -1) if cannot find next position
    const auto matches_error = [error_code](const std::tuple<int, int> &candidate) {
        const auto &[row, col] = candidate;
        return error_matrix[row][col] == error_code;
    };

    if (momentum == "down") {
        std::tuple<int, int> guess;


        //1- down left
        guess = different_down_left(position);
        std::tuple<int, int> guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                momentum = "down";
                std::cout << "down left" << std::endl;
                return guess;
            }
        }

        //2. down_right
        guess = different_down_right(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                momentum = "down";
                std::cout << "down right" << std::endl;
                return guess;
            }
        }

        //3. down
        guess = same_down(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                momentum = "down";
                std::cout << "down" << std::endl;
                return guess;
            }
        }

        //4. left
        guess = same_left(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "left" << std::endl;
                return guess;
            }
        }

        //5. right
        guess = same_right(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "right" << std::endl;
                return guess;
            }
        }

        //6. up_left
        guess = different_up_left(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "up left" << std::endl;
                momentum = "up";
                return guess;
            }
        }

        //7. up_right
        guess = different_up_right(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                    
                std::cout << "up right" << std::endl;
                momentum = "up";
                return guess;
            }
        }

        //8. up
        guess = same_up(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "up" << std::endl;
                momentum = "up";
                return guess;
            }
        }


        return std::make_tuple(-1, -1);

        
    } else if (momentum == "up") {

        //1. down_right
        std::tuple<int, int> guess = different_down_right(position);
        std::tuple<int, int> guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "down right" << std::endl;
                momentum = "down";
                return guess;
            }
        }

        //2. down_left
        guess = different_down_left(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "down left" << std::endl;
                momentum = "down";
                return guess;
            }
        }

        //3. down
        guess = same_down(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );

        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "down" << std::endl;
                momentum = "down";
                return guess;
            }

        }

        //4. right
        guess = same_right(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "right" << std::endl;
                return guess;
            }
        }

        //5. left
        guess = same_left(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){

                std::cout << "left" << std::endl;
                return guess;
            }
        }

        //6 up_right
        guess = different_up_right(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "up right" << std::endl;
                momentum = "up";
                return guess;
            }
        }

        //7 up_left
        guess = different_up_left(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)){
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "up left" << std::endl;
                momentum = "up";
                return guess;
            }

        }

        //8. up
        guess = same_up(position);
        guess_wrapped = std::make_tuple(wrap_index(std::get<0>(guess), N), wrap_index(std::get<1>(guess), D) );
        if (matches_error(guess_wrapped)) {
            if (check_visited(guess_wrapped, path) == -1){
                std::cout << "up" << std::endl;
                momentum = "up";
                return guess;
            }

        }

        return std::make_tuple(-1, -1);

    } else throw std::runtime_error("Momentum not valid");

}

std::tuple<int, int, int> find_first_error(int error_code){
    bool first_found = false;
    int first_i = -1, first_j = -1;
    int number_of_errors = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++){
            if (error_matrix[i][j] == error_code) {
                if (!first_found) {                
                    first_i = i;
                    first_j = j;
                    first_found = true;
                }
                number_of_errors += 1;
            }
        }
    }
    return std::make_tuple(first_i, first_j, number_of_errors);
}


int check_trivial_loop(const int error_code){
    if (error_code <= 0 || error_code > 3) throw std::runtime_error("Error_code not valid");
    std::vector<std::tuple<int, int>> path;
    std::tuple<int, int> next_wrapped;
    std::tuple<int, int> first_error;
    int number_of_errors = 0;
    std::string momentum = "down"; 
    bool first_found = false;

    std::tie(std::get<0>(first_error), std::get<1>(first_error), number_of_errors) = find_first_error(error_code);
    next_wrapped = first_error;


    std::cout << "Number of errors: " << number_of_errors << std::endl;
    std::cout << "Starting position: (" << std::get<0>(next_wrapped) << ", " << std::get<1>(next_wrapped) << ")" << std::endl;

    bool found_loop = true;
    int iter = 0;
    int pacman_effect_count_vertical = 0;
    int pacman_effect_count_horizontal = 0;
    int non_trivial_vertical_loop_count = 0;
    int non_trivial_horizontal_loop_count = 0;

    while (number_of_errors > 0 && iter < 100) {
        iter ++;
        if (found_loop) {
            std::tie(std::get<0>(first_error), std::get<1>(first_error), number_of_errors) = find_first_error(error_code);

            found_loop = false;
        }



        std::tuple<int, int> next = choose_next(next_wrapped, error_code, momentum, path);
        //cannot procede, error in code, 
        next_wrapped = std::make_tuple(wrap_index(std::get<0>(next), N), wrap_index(std::get<1>(next), D) );

        if (std::get<0>(next_wrapped) == -1 || std::get<1>(next_wrapped) == -1) throw std::runtime_error("Error in code, cannot find next error in path, returned (-1,-1)");

        std::cout << "Chosen next: (" << std::get<0>(next_wrapped) << ", " << std::get<1>(next_wrapped) << ")" << std::endl;
        insert_in_path(next_wrapped, path);
        if (std::get<0>(next) != std::get<0>(next_wrapped)){
            pacman_effect_count_vertical ++;
        }

        if (std::get<1>(next) != std::get<1>(next_wrapped)){
            pacman_effect_count_horizontal ++;
        }
        if (std::get<0>(next_wrapped) == std::get<0>(first_error) && std::get<1>(next_wrapped) == std::get<1>(first_error)){
            //non trivial loop
            if (pacman_effect_count_vertical % 2 == 1) {
                non_trivial_vertical_loop_count ++;
            } 
            if (pacman_effect_count_horizontal % 2 == 1) {
                non_trivial_horizontal_loop_count ++;
            }

            remove_errors(path);
            number_of_errors -= path.size();
            path.clear();
            found_loop = true;
        }

    }

    if (non_trivial_vertical_loop_count % 2 == 1 && non_trivial_horizontal_loop_count % 2 == 1) {
        std::cout << "Non trivial vertical and horizontal loop found, horizontal count: " \
         << non_trivial_horizontal_loop_count << std::endl << "vertical count: " << non_trivial_vertical_loop_count << std::endl;
        return -1;
    }

    if (non_trivial_vertical_loop_count % 2 == 1) {
        std::cout << "Non trivial vertical loop found, count: " << non_trivial_vertical_loop_count << std::endl;
        return -1;
    }

    if (non_trivial_horizontal_loop_count % 2 == 1) {
        std::cout << "Non trivial horizontal loop found, count: " << non_trivial_horizontal_loop_count << std::endl;
        return -1;
    }

    return 0;
}

 
int main() {
    const std::string syndrome_path = "syndrome.json";
    const std::string correction_path = "correction.json";
    const bool use_asymmetric_weighted_mode = true; // false => syndrome normale + decoder normale
    const double horizontal_weight = 2.0, vertical_weight = 1.0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
            error_matrix[i][j] = 0;
            correction_matrix[i][j] = 0;
        }
    }

    // Esempio deterministico; puoi sostituire con generate_random_error_matrix(px, pz)
    error_matrix[0][0] = 1;
    error_matrix[2][0] = 1;
    error_matrix[3][1] = 1;
    error_matrix[4][1] = 1;

    error_matrix[0][2] = 1;
    error_matrix[2][2] = 1;
    error_matrix[3][3] = 1;
    error_matrix[4][3] = 1;



    // generate_random_error_matrix(px, pz);

    if (use_asymmetric_weighted_mode) {
        generate_syndrome_matrices_asymmetric();
        print_matrix(error_matrix, "Error matrix before correction asym:");
        print_matrix(syndrome_plaquette, "Syndrome plaquette asym:");
        print_matrix(syndrome_cross, "Syndrome cross asym:");
    } else {
        generate_syndrome_matrices();
        print_matrix(error_matrix, "Error matrix before correction sym:");
        print_matrix(syndrome_plaquette, "Syndrome plaquette sym:");
        print_matrix(syndrome_cross, "Syndrome cross sym:");
    }

    if (!write_syndromes_to_json(syndrome_path)) {
        fprintf(stderr, "Failed to write syndrome JSON to %s\n", syndrome_path.c_str());
        return 1;
    }

    int decoder_status = 0;
    if (use_asymmetric_weighted_mode) {
        std::string weighted_script_path = "decode_with_pymatching_weighted.py";
        if (!std::filesystem::exists(weighted_script_path)) {
            weighted_script_path = "../decode_with_pymatching_weighted.py";
        }
        const std::string weighted_cmd =
            "MPLCONFIGDIR=/tmp python3 " + weighted_script_path + " " + syndrome_path + " " + correction_path +
            " " + std::to_string(horizontal_weight) + " " + std::to_string(vertical_weight);
        decoder_status = std::system(weighted_cmd.c_str());
    } else {
        decoder_status = run_python_decoder(syndrome_path, correction_path);
    }

    if (decoder_status != 0) {
        fprintf(stderr, "Python decoder failed. Check pymatching installation.\n");
        return 1;
    }

    if (!read_correction_matrix_from_json(correction_path)) {
        fprintf(stderr, "Failed to read correction matrix from %s\n", correction_path.c_str());
        return 1;
    }

    print_matrix(correction_matrix, "Correction matrix from pymatching:");
    apply_correction();

    print_matrix(error_matrix, "Residual matrix after applying correction:");

    if (use_asymmetric_weighted_mode) {
        generate_syndrome_matrices_asymmetric();
    } else {
        generate_syndrome_matrices();
    }
    
    check_trivial_loop(1);
    check_trivial_loop(2);

    print_matrix(error_matrix, "Residual matrix after removing trivial loops:");

    generate_syndrome_matrices();
    print_matrix(syndrome_plaquette, "Residual syndrome plaquette:");
    print_matrix(syndrome_cross, "Residual syndrome cross:");

    return 0;
}
