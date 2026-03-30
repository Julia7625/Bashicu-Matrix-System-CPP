#include <iostream>
#include <string>
#include <cstdint>
#include <vector>
#include <initializer_list>
#include <algorithm>

typedef std::vector<std::vector<uint64_t>> mat_type;

template <typename T>
void concatenate(std::vector<T>& a, const std::vector<T>& b) {
    // in-place
    a.insert(a.end(), b.begin(), b.end());
}

void ascend(mat_type& matrix, const std::vector<uint64_t>& delta, const mat_type& ascension_matrix) {
    // adds to delta to each column of bad part where ascension_matrix is 1
    for (uint64_t col = 0; col < matrix.size(); col++) {
        for (uint64_t row = 0; row < matrix[col].size(); row++) {
            matrix[col][row] += delta[row] * ascension_matrix[col][row];
        }
    }
}

struct BMS {
    // Supports Bashicu Matrix System 2.3 (4) and 3.3
    mat_type matrix;
    uint64_t rows;
    uint64_t columns;

    BMS(std::initializer_list<std::initializer_list<uint64_t>> m) {
        matrix.reserve(columns = m.size());
        bool check = true;
        for (const auto& column : m) {
            if (check) {
                rows = column.size();
                check = false;
            }
            matrix.emplace_back(column);
        }
    }

    BMS(mat_type m) : matrix(m) {
        rows = m[0].size();
        columns = m.size();
    }

    std::string to_string() const noexcept {
        std::string result;
        result.reserve(columns*(rows*2+1));
        for (const std::vector<uint64_t>& column : matrix) {
            result += "(";
            for (uint64_t i = 0; i < rows - 1; i++) {
                result += std::to_string(column[i]) + ",";
            }
            result += std::to_string(column[rows - 1]) + ")";
        }
        return result;
    }

    std::pair<uint64_t, bool> find_parent(uint64_t column, uint64_t row) const {
        // returns the parent column and true if there is a parent, 0 and false otherwise
        uint64_t parent = column;
        while (true) {
            if (row == 0) {
                // row = 0: search every column
                if (parent == 0) {
                    return {0, false};  // parent not found
                }
                parent--;
            } else {
                // row > 0: only search on the ancestors of the upper row
                bool has_parent;
                std::tie(parent, has_parent) = find_parent(parent, row - 1);
                if (!has_parent) return {0, false};  // parent not found
            }

            if (matrix[parent][row] < matrix[column][row]) {
                return {parent, true};
            }
        }
    }

    void expand(uint64_t n, int version = 23) {
        // Expands with BM2.3 = BM4 if version = 23 otherwise uses BM3.3
        uint64_t first_zero = 0;
        for (const uint64_t& num : matrix.back()) {
            if (num == 0) break;
            first_zero++;
        }

        if (first_zero == 0) {
            // the last column is (0,0,...,0,0) so we just remove it
            return;
        }
        // first_zero is the index of the first zero in the last column of the matrix

        uint64_t bad_root = find_parent(columns - 1, first_zero - 1).first;

        std::vector<uint64_t> last_column = matrix.back();
        matrix.pop_back();
        columns--;

        mat_type bad_part = std::vector(matrix.begin()+bad_root, matrix.end());
        std::vector<uint64_t> delta(rows, 0);
        // delta = last column - bad root with all rows at or below the last nonzero row of the last column set to zero
        for (uint64_t i = 0; i < first_zero - 1; i++) {
            delta[i] = last_column[i] - matrix[bad_root][i];
        }
        mat_type ascension_matrix;
        ascension_matrix = calculate_ascension_matrix(bad_root, columns, rows, first_zero, version);
        for (uint64_t i = 1; i <= n; i++) {
            ascend(bad_part, delta, ascension_matrix);
            concatenate(matrix, bad_part);
        }
        columns += (columns - bad_root)*n;
    }

    mat_type calculate_ascension_matrix(uint64_t bad_root, uint64_t columns, uint64_t rows, uint64_t first_zero, int version=23) const {
        uint64_t bad_length = columns - bad_root;
        if (first_zero <= 2 || bad_length == 1) {
            // shortcut: ascension matrix doesn't matter until 3 rows and 2 columns
            return mat_type(bad_length, std::vector<uint64_t>(rows, 1));
        }
        mat_type ascension_matrix{bad_length, std::vector<uint64_t>(rows, 0)};
        ascension_matrix[0][0] = 1;  // edge case the loops won't encounter

        bool descendant_of_bad_root = true;
        uint64_t bad_root_bottom = matrix[bad_root][rows - 1];
        for (uint64_t col = bad_root+1; col < columns; col++) {
            // Search each column check if rows have parent in bad part
            if (descendant_of_bad_root && matrix[col][rows - 1] <= bad_root_bottom) {
                descendant_of_bad_root = false;
            }

            ascension_matrix[col - bad_root][0] = 1;  // top row cannot have 0 at ascension matrix
            
            for (uint64_t row = 1; row < first_zero - 1; row++) {
                ascension_matrix[0][row] = 1;  // ascension matrix is 1 at the bad root
                if (descendant_of_bad_root) {
                    ascension_matrix[col - bad_root][row] = 1;  // descendants of the bad root by the lowest nonzero row always ascend
                    continue;
                }

                bool has_parent;
                uint64_t parent;
                std::tie(parent, has_parent) = find_parent(col, row);
                if (!has_parent) continue;

                bool condition = version == 23 ? parent >= bad_root : parent > bad_root;
                if (condition && ascension_matrix[parent - bad_root][row] == 1) {
                    ascension_matrix[col - bad_root][row] = 1;
                }
            }

        }
        return ascension_matrix;
    }
};

uint64_t bms_value(BMS m, uint64_t n, int version=23) {
    while (!m.matrix.empty()) {
        n *= n;
        m.expand(n, version);
    }
    return n;
}
