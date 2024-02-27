#include "include\NSMatrixAssembly.hpp"
#include <iostream>

int main() {
    // Define a sparse matrix of type double
    SpD mat(5, 5); // 5x5 sparse matrix

    // Fill the matrix with some values
    mat.insert(0, 1) = 10.0;
    mat.insert(1, 2) = 20.0;
    mat.insert(2, 3) = 30.0;
    mat.insert(3, 4) = 40.0;

    // Compress the matrix to remove any gaps in the data structure
    mat.makeCompressed();

    // Iterate over the outer dimension (columns for a column-major matrix)
    for (int k = 0; k < mat.outerSize(); ++k) {
        // Iterate over the non-zero elements of the column/row
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
            // it.row() and it.col() provide the row and column indices respectively
            // For a column-major matrix, it.row() is the inner index, and k is the outer index
            // it.value() provides the value of the non-zero element
            std::cout << "Element at (" << it.row() << "," << it.col() << ") = " << it.value() << std::endl;
        }
    }

    return 0;
}
