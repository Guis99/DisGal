#include "include/MatrixAssembly.hpp"

#include <iostream>
#include <fstream>
#include <vector>

using namespace QTM;

int main(int argc, char* argv[]) {

    uint64_t numThreads = std::stoi(argv[argc-1]);
   
    int nx;
    int ny;
    int deg;
    Real_b Lx;
    Real_b Ly;

    // static_assert(sizeof(Real_b) == 4);

    std::cout<<"\n-----------------------------------\n";
    std::cout<<"| Initializing dG Poisson problem |";
    std::cout<<"\n-----------------------------------\n";

    nx = std::stoi(argv[2]); ny = std::stoi(argv[3]); deg = std::stoi(argv[1]);
    Lx = std::stod(argv[4]); Ly = std::stod(argv[5]);

    Real_b penalty = std::stod(argv[6]); 

    std::string source = argv[7];

    std::cout<<"Order (polynomial degree): "<<deg<<"\nnx: "<<nx<<"\nny: "<<ny<<"\nLx: "<<Lx<<"\nLx: "<<Ly<<std::endl;
    std::cout<<"Source term: "<<source<<std::endl;

    int numBoundaries = std::stoi(argv[8]);
    std::vector<std::string> bcs;

    std::vector<bool> essentialBC;
    std::vector<bool> naturalBC;

    int currIdx = 9;
    int numEss = 0;
    int numNat = 0;

    for (int i=0; i<numBoundaries; i++) {
        if (argv[currIdx+i][0] == '1') {
            numEss++;
            essentialBC.push_back(1);
        } else {
            essentialBC.push_back(0);
        }
    }
    currIdx += numBoundaries;

    for (int i=0; i<numBoundaries; i++) {
        if (argv[currIdx+i][0] == '1') {
            numNat++;
            naturalBC.push_back(1);
        } else {
            naturalBC.push_back(0);
        }
    }
    currIdx += numBoundaries;

    std::vector<std::string> dbcs;
    std::vector<std::string> nbcs;

    std::cout<<"Boundary conditions for Dirichlet boundary condition:\n";
    for (int i=0; i<numEss; i++) {
        std::cout<<argv[i+currIdx]<<std::endl;
        dbcs.push_back(argv[i+currIdx]);
    }
    currIdx += numEss; 
    std::cout<<"-------------------------\n";

    std::cout<<"Boundary conditions for Neumann boundary condition:\n";
    for (int i=0; i<numNat; i++) {
        std::cout<<argv[i+currIdx]<<std::endl;
        nbcs.push_back(argv[i+currIdx]);
    }
    currIdx += numNat; 
    std::cout<<"-------------------------\n";

    QuadTreeMesh mesh(deg, nx, ny, Lx, Ly);         

    Real_b c = 1;
    Real_b k = 1;

    DvD z = PMA::dgPoissonSolve(mesh, k, source, essentialBC, naturalBC, dbcs, nbcs, penalty, 0, numThreads);

    std::vector<std::array<Real_b,2>> allNodePos = mesh.nodePositions;

    std::ofstream fileOut("outputdG.txt");

    std::cout<<"Writing results to file"<<std::endl;
    if (fileOut.is_open()) {
        for (size_t i = 0; i < allNodePos.size(); ++i) {
            // Extract x and y from the coordinates vector
            Real_b x = allNodePos[i][0];
            Real_b y = allNodePos[i][1];

            // Write x, y, and z to the file separated by commas
            fileOut << x << "," << y << "," << z(i) << std::endl;
        }

        // Close the file
        fileOut.close();
    }

    std::ofstream outFile("quadtreedG.json");
    mesh.exportToJson(outFile);
}
