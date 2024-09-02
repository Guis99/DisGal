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

    std::cout<<argc<<"\n-----------------------------------\n";
    std::cout<<"| Initializing dG Poisson problem |";
    std::cout<<argc<<"\n-----------------------------------\n";

    nx = std::stoi(argv[2]); ny = std::stoi(argv[3]); deg = std::stoi(argv[1]);
    Lx = std::stod(argv[4]); Ly = std::stod(argv[5]);

    double penalty = std::stod(argv[6]); 

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

    double c = 1;
    double k = 1;

    DD z = PMA::dgPoissonSolve(mesh, k, source, essentialBC, naturalBC, dbcs, nbcs, penalty, 0);

    std::vector<std::array<double,2>> allNodePos = mesh.nodePositions;

    std::ofstream fileOut("outputdG.txt");

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

    std::ofstream outFile("quadtreedG.json");
    mesh.exportToJson(outFile);
}