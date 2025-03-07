#include "common.hpp"

// std::vector<DvD> generalExplicitRK(SpD& A,
//                                     SpD& columnSpace, SpD& nullSpace,
//                                     DvD& boundaryVals,
//                                     DvD& initialCondition,
//                                     double timeStep,
//                                     int numTimeSteps,
//                                     int s, DD a,
//                                     std::vector<double> b, std::vector<double> c) {
//     std::vector<DvD> out;
//     std::vector<DvD> kOps; kOps.reserve(s);
//     SpD fPropOp;
//     kOps[0] = A;
//     for (int i=1; i<s; i++) {
//         for (int j=0; j<i-1; j++) {

//         }

//     }

//     DvD prevState = initialCondition;

//     DvD x;
//     for (int i=1; i<numTimeSteps; i++) {
//         x = fPropOp * x;
//         prevState = x;
//         out.push_back(columnSpace * x);
//     }
// }

// std::vector<DvD> generalImplicitRK(SpD& A,
//                                     SpD& columnSpace, SpD& nullSpace,
//                                     DvD& boundaryValues,
//                                     DvD& initialCondition,
//                                     double timeStep,
//                                     int numTimeSteps,
//                                     int s, DD a, std::vector<double> b, std::vector<double> c) {
//     std::vector<DvD> out;
//     std::vector<DvD> k; k.reserve(s);
//     SpD fPropOp;
//     k[0] = A * initialCondition;
//     DvD prevState = initialCondition;

//     DvD x;
//     for (int i=1; i<numTimeSteps; i++) {
//         x = fPropOp * x;
//         prevState = x;
//         out.push_back(columnSpace * x);
//     }
// }

int main(int argc, char* argv[]) {
    int deg;
    int numElem;
    double width;

    std::string IC;

    double timeLength;
    int numTimeSteps;

    std::cout<<"-----------------------------------------\n";
    std::cout<<"| Initializing 1D linear advection solver |"<<std::endl;
    std::cout<<"-----------------------------------------\n";

    deg = std::stoi(argv[1]); numElem = std::stoi(argv[2]); width = std::stod(argv[3]);
    IC = argv[4];
    timeLength = std::stod(argv[5]); numTimeSteps = std::stoi(argv[6]);
    double timeStep = timeLength / numTimeSteps; 
    
    int integratorType = std::stoi(argv[7]);
    auto baseline = argv[8];
    double cutoff = std::stod(argv[9]);

    std::vector<std::string> integrators = {"Forward Euler",
                                            "Crank-Nicholson",
                                            "RK4",
                                            "GL1",
                                            "GL2"};

    std::cout<<"ts: "<<timeStep<<", tl: "<<timeLength<<", nts: "<<numTimeSteps<<std::endl;

    std::cout<<"Order (polynomial degree): "<<deg<<"\nnx: "<<numElem<<std::endl;
    std::cout<<"-------------------------"<<std::endl;
    std::cout<<"Initial condition: "<<IC<<std::endl;

    std::cout<<"Simulating for duration "<<timeLength<<" over "<<numTimeSteps<<" timesteps using "<<integrators[integratorType]<<" integration"<<std::endl;
    std::cout<<"-------------------------"<<std::endl;
    
    double spacing = width / numElem;

    std::cout<<std::endl<<"Generating mesh..."<<std::endl;
    std::cout<<"-------------------------"<<std::endl;
    BasicMesh1D mesh = CreateMesh(deg, numElem, spacing);
    auto nodes = mesh.nodes;
    int nNodes = nodes.size();

    DEBUG_PRINT("ALL NODE POS");
    for (auto node : nodes) {
        DEBUG_PRINT(node.position);
    }

    DvD initialCondition = getICVec(mesh, nNodes, baseline, cutoff, IC);
    DEBUG_PRINT("ic size: ",initialCondition.rows()," x ", initialCondition.cols());

    std::vector<int> boundaryNodes = {};
    std::vector<int> freeNodes; freeNodes.reserve(nNodes);

    for (int i=0; i<nNodes; i++) {
        freeNodes.push_back(i);
    }

    int numDofs = numElem * (deg + 1);
    std::cout<<"Generating system matrices..."<<std::endl;
    SpD MassMatrix(numDofs, numDofs); SpD StiffnessMatrix(numDofs, numDofs);
    AssembleSystemMatrices(mesh, MassMatrix, StiffnessMatrix);

    SpD nullSpace(nNodes, boundaryNodes.size());
    SpD columnSpace(nNodes, freeNodes.size());
    std::cout<<"Assembling null-orth matrix..."<<std::endl;
    GetExtensionMatrices(mesh, boundaryNodes, freeNodes, nullSpace, columnSpace);

    std::vector<double> oneVal = {0};
    Eigen::Map<DvD> boundaryValsMap(oneVal.data(), 1,  1);
    DvD boundaryVals = (DvD) boundaryValsMap;

    std::cout<<"Computing solution..."<<std::endl;
    std::vector<DvD> solns = ComputeTransientSolution(StiffnessMatrix, MassMatrix, 
                                                        columnSpace, 
                                                        nullSpace, initialCondition, 
                                                        timeStep, numTimeSteps,
                                                        integratorType);

    std::vector<double> xPositions;
    xPositions.reserve(nNodes);
    for (int i=0; i<nNodes; i++) {
        xPositions.push_back(nodes[i].position);
    }

    Eigen::Map<DD> xOffsets(xPositions.data(), nNodes, 1);
    std::ofstream fileX("xt.txt");
    std::ofstream fileZ("zt.txt");

    if (fileZ.is_open())
    {
        for (auto x : solns) {
            fileZ << x;
            fileZ << "\n";
        }
    }
    if (fileX.is_open())
    {
        fileX << xOffsets;
    }

    return 0;
}