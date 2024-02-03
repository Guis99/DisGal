#include "include\MatrixAssembly.hpp"

#include <iostream>
#include <fstream>
#include <vector>

using namespace QTM;

int main(int argc, char* argv[]) {
    int nx;
    int ny;
    int deg;
    double Lx;
    double Ly;

    std::cout<<argc<<std::endl;

    for (int i=0; i<argc; i++) {
        std::cout<<argv[i]<<std::endl;
    }


    // nx = std::stoi(argv[2]); ny = std::stoi(argv[3]); deg = std::stoi(argv[1]);
    // Lx = std::stod(argv[4]); Ly = std::stod(argv[5]);

    // std::string source = argv[6];
    // std::vector<std::string> bcs;

    // for (int i=7; i<11; i++) {
    //     bcs.push_back(argv[i]);
    // }

    // int widthX = nx*(deg+1);
    // int widthY = ny*(deg+1);

    // QuadTreeMesh mesh(deg, nx, ny, Lx, Ly);   

    // double c = 1;
    // double k = 1;
    // double f = 0;

    // DD x = PoissonSolve(mesh, c, k, source, bcs);

    // std::vector<double> xGrid;
    // std::vector<double> yGrid;
    // xGrid.reserve(widthX);
    // yGrid.reserve(widthY);

    // Eigen::Map<DD> xOffsets(xGrid.data(), widthX, 1);
    // Eigen::Map<DD> yOffsets(yGrid.data(), widthY, 1);

    // std::ofstream fileOut("output.txt");

    // if (fileOut.is_open()) {
    //     fileOut << x;
    // }

}