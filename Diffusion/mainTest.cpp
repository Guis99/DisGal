#include "include/MatrixAssembly.hpp"

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

    std::cout<<argc<<"\n-------------\n";
    
    nx = std::stoi(argv[2]); ny = std::stoi(argv[3]); deg = std::stoi(argv[1]);
    Lx = std::stod(argv[4]); Ly = std::stod(argv[5]);

    std::string source = argv[6];
    std::vector<std::string> bcs;

    std::cout<<"Order (polynomial degree): "<<deg<<"\nnx: "<<nx<<"\nny: "<<ny<<"\nLx: "<<Lx<<"\nLx: "<<Ly<<std::endl;
    std::cout<<"Source term: "<<source<<std::endl;

    std::cout<<"Boundary conditions:\n";
    for (int i=7; i<11; i++) {
        std::cout<<argv[i]<<std::endl;
        bcs.push_back(argv[i]);
    }

    double penalty = std::stod(argv[11]); 

    QuadTreeMesh mesh(deg, nx, ny, Lx, Ly);   

    // mesh.Refine(mesh.leaves);
    // mesh.Refine(mesh.leaves);

    double c = 1;
    double k = 1;

    DD z = PMA::PoissonSolve(mesh, c, k, source, bcs, penalty);

    std::vector<std::array<double,2>> allNodePos = mesh.nodePositions;

    std::ofstream fileOut("output.txt");

    std::cout<<"Writing results to file"<<std::endl;
    if (fileOut.is_open()) {
        for (size_t i = 0; i < allNodePos.size(); ++i) {
            // Extract x and y from the coordinates vector
            double x = allNodePos[i][0];
            double y = allNodePos[i][1];

            // Write x, y, and z to the file separated by commas
            fileOut << x << "," << y << "," << z(i) << std::endl;
        }

        // Close the file
        fileOut.close();
    }

    std::ofstream outFile("quadtree.json");
    mesh.exportToJson(outFile);
}