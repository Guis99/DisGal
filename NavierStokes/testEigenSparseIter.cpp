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

    int nNodes = 3;
    DvD ones = DvD::Ones(nNodes);
    DvD twos = 2*DvD::Ones(nNodes);
    DvD out(3*nNodes+1);
    DvD allZeros = DvD::Zero(nNodes+1);

    out << ones, twos, allZeros;

    std::cout<<out<<std::endl;

    SpD mat1(nNodes, nNodes);
    SpD mat2(nNodes, nNodes);

    std::vector<Eigen::Triplet<double>> tripletList1;
    std::vector<Eigen::Triplet<double>> tripletList2;

    int num=1;
    for (int j=0; j<nNodes; j++) {
        for (int i=0; i<nNodes; i++) {
            tripletList1.emplace_back(i,j,num);
            tripletList2.emplace_back(i,j,2*num);
            num++;
        }
    }

    mat1.setFromTriplets(tripletList1.begin(), tripletList1.end());
    mat2.setFromTriplets(tripletList2.begin(), tripletList2.end());

    SpD mat1T = (SpD)(mat1.transpose());
    SpD mat2T = (SpD)(mat2.transpose());


    std::cout<<"Testing matrix compose"<<std::endl;
    std::cout<<mat1<<"-----------------\n"<<mat2<<std::endl;

    SpD bigMat1(3*nNodes, 3*nNodes);
    SpD bigMat2(3*nNodes, 3*nNodes);

    std::vector<Eigen::Triplet<double>> btripletList1;
    std::vector<Eigen::Triplet<double>> btripletList2;

    int offset2 = 2*nNodes;
    int offset1 = nNodes;

    for (int k = 0; k < mat1.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat1, k); it; ++it) {
            btripletList1.emplace_back(it.row(), it.col()+offset2, it.value());
            btripletList1.emplace_back(it.col()+offset2, it.row(), it.value());
        }
    }

    for (int k = 0; k < mat2.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat2, k); it; ++it) {
            btripletList1.emplace_back(it.row()+offset1, it.col()+offset2, it.value());
            btripletList1.emplace_back(it.col()+offset2, it.row()+offset1, it.value());
        }
    }

    // transpose with 4 separate mats

    for (int k = 0; k < mat1.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat1, k); it; ++it) {
            btripletList2.emplace_back(it.row(), it.col()+offset2, it.value());
        }
    }

    for (int k = 0; k < mat2.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat2, k); it; ++it) {
            btripletList2.emplace_back(it.row()+offset1, it.col()+offset2, it.value());
        }
    }

    for (int k = 0; k < mat1T.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat1T, k); it; ++it) {
            btripletList2.emplace_back(it.row()+offset2, it.col(), it.value());
        }
    }

    for (int k = 0; k < mat2T.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat2T, k); it; ++it) {
            btripletList2.emplace_back(it.row()+offset2, it.col()+offset1, it.value());
        }
    }

    bigMat1.setFromTriplets(btripletList1.begin(), btripletList1.end());
    bigMat2.setFromTriplets(btripletList2.begin(), btripletList2.end());


    std::cout<<"Testing matrix compose"<<std::endl;
    std::cout<<bigMat1<<"-----------------\n"<<bigMat2<<"-----------------\n"<<bigMat1-bigMat2<<std::endl;

    

    return 0;
}
