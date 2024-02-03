#include "include\MatrixAssembly.hpp"

#include <iostream>
#include <fstream>
#include <vector>

using namespace QTM;

int main(int argc, char* argv[]) {
    int nx;
    int ny;
    int deg;

    std::cout<<argc<<std::endl;

    for (int i=0; i<argc; i++) {
        std::cout<<argv[i]<<std::endl;
    }


    nx = std::stoi(argv[1]); ny = std::stoi(argv[2]); deg = std::stoi(argv[3]);

    // std::cin>>nxElem;
    // std::cin>>nyElem;
    // std::cin>>xdeg;
    // std::cin>>ydeg;
    // std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    int widthX = nx*(deg+1);
    int widthY = ny*(deg+1);

    QuadTreeMesh mesh(deg, nx, ny, Lx, Ly);   

    double c = 1;
    double k = 1;
    double f = 0;

    DD x = PoissonSolve(mesh, c, k, f);

    std::vector<double> xGrid;
    std::vector<double> yGrid;
    xGrid.reserve(widthX);
    yGrid.reserve(widthY);

    for (int i=0; i<widthX; i++) {
        auto xcoords = mesh.posOfNodes(std::vector<int>{i});
        xGrid.push_back(xcoords[0][0]);
    } 

    for (int i=0; i<widthY; i++) {
        int yIdx = i*widthX;
        auto ycoords = mesh.posOfNodes(std::vector<int>{yIdx});
        yGrid.push_back(ycoords[0][1]);
    } 

    Eigen::Map<DD> xOffsets(xGrid.data(), widthX, 1);
    Eigen::Map<DD> yOffsets(yGrid.data(), widthY, 1);

    std::ofstream fileOut("output.txt");

    if (fileOut.is_open()) {
        fileOut << x;
    }

}