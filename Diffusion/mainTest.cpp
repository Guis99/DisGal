#include "include\MatrixAssembly.hpp"

#include <iostream>
#include <fstream>
#include <vector>

using namespace QTM;

void exportToJson(std::shared_ptr<Cell> root, std::ostream& out, int indent);
void exportToJson(std::shared_ptr<Cell> root, std::ostream& out, int indent) {
    std::string tabs;
    for (int i=0; i<indent; i++) {
        tabs += "\t";
    }
    out << tabs << "{";
    out << "\"x\": " << root->center[0] << ",\n ";
    out << tabs << "\"y\": " << root->center[1] << ",\n ";
    out << tabs << "\"width\": " << 2*root->width << ",\n ";
    out << tabs << "\"level\": " << root->level << ",\n ";
    out << tabs << "\"isLeaf\": " << (root->isLeaf() ? "true" : "false") << ",\n ";
    out << tabs << "\"CID\": " << root->CID << ",\n ";
    out << tabs << "\"children\": [\n";
    for (int i = 0; i < 4; ++i) {
        if (root->children[i] != nullptr) {
            if (i > 0) out << ",\n ";
            exportToJson(root->children[i], out, indent+1);
        }
    }
    out << tabs << "]}";
}

int main(int argc, char* argv[]) {
    int nx;
    int ny;
    int deg;
    double Lx;
    double Ly;

    std::cout<<argc<<"\n-------------\n";

    // for (int i=0; i<argc; i++) {
    //     std::cout<<argv[i]<<std::endl;
    // }


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

    int penalty = std::stoi(argv[11]); 

    QuadTreeMesh mesh(deg, nx, ny, Lx, Ly);   

    double c = 1;
    double k = 1;

    DD z = PoissonSolve(mesh, c, k, source, bcs, penalty);

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
    outFile << "{";
    outFile << "\"x\": " << 0 << ",\n ";
    outFile << "\"y\": " << 0 << ",\n ";
    outFile << "\"width\": " << 0 << ",\n ";
    outFile << "\"level\": " << 0 << ",\n ";
    outFile << "\"isLeaf\": " << "false" << ",\n ";
    outFile << "\"CID\": " << 0 << ",\n ";
    outFile << "\"children\": [\n";

    std::cout<<"Writing mesh data to file\n";
    for (auto cell : mesh.topCells) {
        exportToJson(cell, outFile, 1);
        if (cell != mesh.topCells.back()) {
            outFile << ", ";
        }
    }
    outFile << "]}";
}