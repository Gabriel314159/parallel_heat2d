#pragma once
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace heat2d {
class Grid {
public:
    int nx, ny, pitch;
    Grid(int columns, int rows): nx(columns), ny(rows), pitch(columns+2),
        data_(static_cast<std::size_t>(columns+2)*static_cast<std::size_t>(rows+2), 0.0) {
        if (nx<1 || ny<1) throw std::invalid_argument("Empty local grid.");
    }
    double* row(int i) { return data_.data()+static_cast<std::size_t>(i)*pitch; }
    const double* row(int i) const { return data_.data()+static_cast<std::size_t>(i)*pitch; }
    double& operator()(int i, int j) { return row(i)[j]; }
    double operator()(int i, int j) const { return row(i)[j]; }
    void swap_data(Grid& other) { data_.swap(other.data_); }
private:
    std::vector<double> data_;
};
}
