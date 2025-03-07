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

    double resNorm = PMA::dgComputeResidual(mesh, k, source, essentialBC, naturalBC, dbcs, nbcs, penalty, 0);

    std::vector<std::array<double,2>> allNodePos = mesh.nodePositions;

    std::ofstream fileOut("outputConv.txt");

    std::cout<<"Writing results to file"<<std::endl;
    if (fileOut.is_open()) {
        fileOut << resNorm << std::endl;

        // Close the file
        fileOut.close();
    }
}

double PMA::dgComputeResidual(QTM::QuadTreeMesh& inputMesh,
                double k,
                std::string source,
                std::vector<bool> isDirichletBC,
                std::vector<bool> isNeumannBC,
                std::vector<std::string> dbcs,
                std::vector<std::string> nbcs,
                double penaltyParam,
                double dirichletPenalty) {
    
    auto boundaryNodes = inputMesh.boundaryNodes;
    std::vector<int> freeNodes = inputMesh.freeNodes;
    int nNodes = inputMesh.nNodes();

    std::vector<int> allBoundaryNodes;
    for (auto bNodes : boundaryNodes) {
        allBoundaryNodes.insert(allBoundaryNodes.end(), bNodes.begin(), bNodes.end());
    }

    PMA::quadUtils package = PMA::GenerateAssemblyPackage(inputMesh);

    std::cout<<"Assembling stiffness matrix"<<std::endl;
    SpD KMatrix = PMA::StiffnessMatrix(inputMesh, k);
    std::cout<<"Assembling penalty matrix"<<std::endl;
    SpD PMatrix = PMA::PenaltyMatrix(inputMesh, k, penaltyParam, package);
    std::cout<<"Assembling flux matrix"<<std::endl;
    SpD SMatrix = PMA::FluxMatrix(inputMesh, k, package);
    SpD SMatrixT = (SpD)(SMatrix.transpose());
    std::cout<<"Assembling flux matrix"<<std::endl;
    SpD Miscellaneous = PMA::BoundaryMatrix(inputMesh, k, isDirichletBC, dbcs, penaltyParam, package);
    std::cout<<"Assembling overall system matrix"<<std::endl;
    SpD StiffnessMatrix = KMatrix + PMatrix - SMatrix - SMatrixT + Miscellaneous;

    // std::cout<<"All:/n"<<StiffnessMatrix<<std::endl;
    // std::cout<<"Boundary:/n"<<Miscellaneous<<std::endl;

    std::cout<<"Assembling RHS vector"<<std::endl;
    DvD FMatrix = PMA::AssembleFVec(inputMesh, 1.0, source);

    std::cout<<"Assembling dirichlet boundary condition vector"<<std::endl;
    DvD dirichletBound = PMA::IntegrateDirichlet(inputMesh, isDirichletBC, dbcs, penaltyParam, package);

    std::cout<<"Assembling neumann boundary condition vector"<<std::endl;
    DvD neumannBound = PMA::IntegrateNeumann(inputMesh, isNeumannBC, nbcs, package);

    DvD b = FMatrix + dirichletBound + neumannBound;

    std::cout<<"Solving system with "<<inputMesh.nNodes()<<" nodes"<<std::endl;
    DvD x = PMA::ComputeSolutionStationaryLinearNoElim(StiffnessMatrix, b);
    std::cout<<"System solved!"<<std::endl;

    DvD resVec = StiffnessMatrix * x - b;
    double resNorm = resVec.norm(); 

    return resNorm;
}