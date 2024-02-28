#include "include\NSMatrixAssembly.hpp"

#include <iostream>
#include <fstream>
#include <vector>

using namespace QTM;

void exportToJson(QuadTreeMesh mesh, std::ostream& out) {
    std::vector<std::shared_ptr<Cell>> allLeaves;
    for (auto cell : mesh.topCells) {
        auto leaves = cell->traverse();
        allLeaves.insert(allLeaves.end(), leaves.begin(), leaves.end());
    }
    out << "{";
    out << "\"children\": [\n";

    out << "{";
    out << "\"x\": " << allLeaves[0]->center[0] << ", ";
    out << "\"y\": " << allLeaves[0]->center[1] << ", ";
    out << "\"width\": " << 2*allLeaves[0]->width << ", ";
    out << "\"level\": " << allLeaves[0]->level << ", ";
    out << "\"CID\": " << allLeaves[0]->CID << "}\n ";

    for (int i=1; i<allLeaves.size(); i++) {
        auto leaf = allLeaves[i];
        out <<",{";
        out << "\"x\": " << leaf->center[0] << ", ";
        out << "\"y\": " << leaf->center[1] << ", ";
        out << "\"width\": " << 2*leaf->width << ", ";
        out << "\"level\": " << leaf->level << ", ";
        out << "\"CID\": " << leaf->CID << "}\n ";
    }
    out << "]}";
}

int main(int argc, char* argv[]) {
    std::cout<<"-------------------------------------\n";
    std::cout<<"| Initializing Stokes problem solver |"<<std::endl;
    std::cout<<"-------------------------------------\n";
    int nx;
    int ny;
    int deg;
    double Lx;
    double Ly;

    nx = std::stoi(argv[2]); ny = std::stoi(argv[3]); deg = std::stoi(argv[1]);
    Lx = std::stod(argv[4]); Ly = std::stod(argv[5]); 
    double penalty = std::stod(argv[6]); 

    std::string source;
    std::vector<std::string> bcs;

    std::cout<<"Order (polynomial degree): "<<deg<<"\nnx: "<<nx<<", ny: "<<ny<<"\nLx: "<<Lx<<", Lx: "<<Ly<<std::endl;
    std::cout<<"-------------------------\n";
    std::cout<<"Source terms: "<<std::endl;
    std::vector<std::string> sources;
    std::vector<std::string> dirs = {"x: ","y: "};
    for (int i=7; i<9; i++) {
        source = argv[i];
        std::cout<<dirs[i-7]<<source<<std::endl;
        sources.push_back(source);
    }
    std::cout<<"-------------------------\n";

    std::vector<bool> velEss;
    std::vector<bool> pressEss;
    std::vector<bool> natBC;

    int numBoundaries = std::stoi(argv[9]);
    int currIdx = 10;
    int numVel = 0;
    int numPress = 0;
    int numNat = 0;

    for (int i=0; i<numBoundaries; i++) {
        if (argv[currIdx+i][0] == '1') {
            numVel++;
            velEss.push_back(1);
        } else {
            velEss.push_back(0);
        }
    }
    currIdx += numBoundaries;

    for (int i=0; i<numBoundaries; i++) {
        if (argv[currIdx+i][0] == '1') {
            numPress++;
            pressEss.push_back(1);
        } else {
            pressEss.push_back(0);
        }
    }
    currIdx += numBoundaries;

    for (int i=0; i<numBoundaries; i++) {
        if (argv[currIdx+i][0] == '1') {
            numNat++;
            natBC.push_back(1);
        } else {
            natBC.push_back(0);
        }
    }
    currIdx += numBoundaries;

    std::vector<std::string> Ubcs;
    std::vector<std::string> Vbcs;
    std::vector<std::string> Pbcs;
    std::vector<std::string> Nbcs;

    std::cout<<"Boundary conditions for x-velocity:\n";
    for (int i=0; i<numVel; i++) {
        std::cout<<argv[i+currIdx]<<std::endl;
        Ubcs.push_back(argv[i+currIdx]);
    }
    currIdx += numVel; 
    std::cout<<"-------------------------\n";

    std::cout<<"Boundary conditions for y-velocity:\n";
    for (int i=0; i<numVel; i++) {
        std::cout<<argv[i+currIdx]<<std::endl;
        Vbcs.push_back(argv[i+currIdx]);
    }
    currIdx += numVel; 
    std::cout<<"-------------------------\n";

    std::cout<<"Boundary conditions for pressure:\n";
    for (int i=0; i<numPress; i++) {
        std::cout<<argv[i+currIdx]<<std::endl;
        Pbcs.push_back(argv[i+currIdx]);
    }
    currIdx += numPress; 
    std::cout<<"-------------------------\n";

    std::cout<<"Neumann boundary conditions:\n";
    for (int i=0; i<numNat; i++) {
        std::cout<<argv[i+currIdx]<<std::endl;
        Nbcs.push_back(argv[i+currIdx]);
    }
    currIdx += numNat; 
    std::cout<<"-------------------------\n";
    std::cout<<"SIPG penalty param: "<<penalty<<std::endl;
    std::cout<<"-------------------------\n";

    QuadTreeMesh mesh(deg, nx, ny, Lx, Ly);   
    std::vector<std::vector<int>> bcNodes = mesh.boundaryNodes;

    double rho = 1;
    double mu = 1;

    DvD z = IncompressibleStokesSolve(mesh, rho, mu, penalty, 
                                    sources, Ubcs, Vbcs, Pbcs, Nbcs, 
                                    bcNodes, velEss, pressEss, natBC);

    std::vector<std::array<double,2>> allNodePos = mesh.nodePositions;

    std::ofstream fileOut("output.txt");

    std::cout<<"Writing results to file"<<std::endl;
    int nNodes = mesh.nNodes();
    int offset1 = 2*nNodes;
    if (fileOut.is_open()) {
        for (size_t i = 0; i < allNodePos.size(); ++i) {
            // Extract x and y from the coordinates vector
            double x = allNodePos[i][0];
            double y = allNodePos[i][1];

            // Write x, y, and z to the file separated by commas
            fileOut << x << "," << y << "," << z(i) << "," << z(i+nNodes) << "," << z(i+offset1)<< std::endl;
        }

        // Close the file
        fileOut.close();
    }

    std::ofstream outFile("quadtree.json");
    exportToJson(mesh, outFile);
}