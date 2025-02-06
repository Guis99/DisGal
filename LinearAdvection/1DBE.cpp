#include "common.hpp"

std::vector<DvD> TimeStep::solver_RK4_NonLinear_Burger(SpD &A, 
                            SpD &columnSpace, SpD &nullSpace,
                            DvD &initialCondition,
                            double timeStep,
                            int numTimeSteps) {
    std::vector<DvD> out; out.reserve(numTimeSteps);
    out.push_back(initialCondition);
    DvD prevState = initialCondition;
    DvD x;
    DvD k1; DvD k2; DvD k3; DvD k4;

    // auto systemDeriv = [](SpD& A, DvD&& state) {{
        // return state;}
        // return (DvD)(A * state.asDiagonal() * state);
    // };

    
    for (int i=1; i<numTimeSteps; i++) {
        k1 = A * prevState.asDiagonal() * prevState;
        k2 = A * (prevState + .5 * timeStep * k1).asDiagonal() * (prevState + .5 * timeStep * k1);
        k3 = A * (prevState + .5 * timeStep * k2).asDiagonal() * (prevState + .5 * timeStep * k2);
        k4 = A * (prevState + timeStep * k3).asDiagonal() * (prevState + timeStep * k3);

        // k1 = systemDeriv(A, prevState);
        // k2 = systemDeriv(A, prevState + .5 * timeStep * k1);
        // k3 = systemDeriv(A, prevState + .5 * timeStep * k2);
        // k4 = systemDeriv(A, prevState + timeStep * k3);
        // k2 = A * (prevState + .5 * timeStep * k1);
        // k3 = A * (prevState + .5 * timeStep * k2);
        // k4 = A * (prevState + timeStep * k3);

        x = prevState + (timeStep / 6) * (k1 + 2*k2 + 2*k3 + k4);
        DvD diff = x - prevState;
        // std::cout<<"diff at iter "<<i<<": "<< diff.norm()<<std::endl;
        prevState = x;
        out.push_back(x);
    }

    return out;
}

std::vector<DvD> ComputeTransientSolutionBurger(SpD &K11, 
                                SpD &M11, SpD SIP11,
                                SpD &columnSpace, 
                                SpD &nullSpace, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps,
                                int integrator) {
    // Eliminate boundary rows and columns
    // SpD K11 = columnSpace.transpose() * StiffnessMatrix * columnSpace;
    // SpD M11 = columnSpace.transpose() * MassMatrix * columnSpace;
    // SpD SIP11 = columnSpace.transpose() * SIPMatrix * columnSpace;

    int system_size = K11.rows();
    SpD dummyId(system_size, system_size); 
    std::vector<Eigen::Triplet<double>> tripletListID;

    tripletListID.reserve(system_size);

    for (int i=0; i<system_size; i++) {
        tripletListID.emplace_back(i, i, 1.0);
    }

    dummyId.setFromTriplets(tripletListID.begin(), tripletListID.end());

    // Invert diagonal mass matrix
    Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolverMM;    
    LuSolverMM.analyzePattern(M11);
    LuSolverMM.factorize(M11);

    SpD M11_Inv = LuSolverMM.solve(dummyId);

    // std::cout<<"inverse"<<std::endl<<combinedMatsInv<<std::endl;

    std::vector<DvD> out;
    out.push_back(initialCondition);
    DvD prevState = initialCondition;
    DvD x;
    SpD K1 = .5 * M11_Inv * (K11 - M11);

    // std::cout<<"K1: "<<K1<<std::endl;
    // std::cout<<"M11_Inv: "<<M11_Inv<<std::endl;
    // std::cout<<"K11: "<<K11<<std::endl;
    // std::cout<<"M11: "<<M11<<std::endl;

    switch (integrator) {
        case 0: // Forward Euler
            return out;
            break;
        case 1: // Crank-Nicholson
            return out;
            break;
        case 2: // RK4
            return TimeStep::solver_RK4_NonLinear_Burger(K1, columnSpace, nullSpace, initialCondition, timeStep, numTimeSteps);
            break;
        case 3: // GL1
            return out;
            break;
        case 4: // GL2
            return TimeStep::solver_GL2(K1, columnSpace, nullSpace, initialCondition, timeStep, numTimeSteps);
            break;
    }
    return out;
}

int main(int argc, char* argv[]) {
    int deg;
    int numElem;
    double width;

    std::string IC;

    double timeLength;
    int numTimeSteps;

    std::cout<<"-----------------------------------------\n";
    std::cout<<"| Initializing 1D Burger's Equation solver |"<<std::endl;
    std::cout<<"-----------------------------------------\n";

    deg = std::stoi(argv[1]); numElem = std::stoi(argv[2]); width = std::stod(argv[3]);
    IC = argv[4];
    timeLength = std::stod(argv[5]); numTimeSteps = std::stoi(argv[6]);
    double timeStep = timeLength / numTimeSteps; 
    
    int integratorType = std::stoi(argv[7]);
    std::vector<std::string> integrators = {"Forward Euler",
                                            "Crank-Nicholson",
                                            "RK4",
                                            "GL1",
                                            "GL2"};

    int systemForm = std::stoi(argv[8]); // 0 = fully non-linear state vector, 1 = pseudo-linear matrix form

    char* baseline = argv[9];
    double cutoff = std::stod(argv[10]);

    std::cout<<"Order (polynomial degree): "<<deg<<"\nnx: "<<numElem<<std::endl;
    std::cout<<"-------------------------"<<std::endl;
    std::cout<<"Initial condition: "<<IC<<std::endl;

    std::cout<<"Simulating for duration "<<timeLength<<" over "<<numTimeSteps<<" timesteps using "<<integrators[integratorType]<<" integration"<<std::endl;
    std::cout<<"-------------------------"<<std::endl;
    
    double spacing = width / numElem;

    std::cout<<"Generating mesh..."<<std::endl;
    BasicMesh1D mesh = CreateMesh(deg, numElem, spacing);
    auto nodes = mesh.nodes;
    int nNodes = nodes.size();
    
    DvD initialCondition = getICVec(mesh, nNodes, baseline, cutoff, IC);
    DEBUG_PRINT(initialCondition);

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
    std::vector<DvD> solns = ComputeTransientSolutionBurger(StiffnessMatrix, MassMatrix, 
                                                        MassMatrix,
                                                        columnSpace, 
                                                        nullSpace, 
                                                        initialCondition, timeStep, numTimeSteps,
                                                        integratorType);

    
    std::cout<<"Success!"<<std::endl;
    std::vector<double> xPositions;
    xPositions.reserve(nNodes);
    for (int i=0; i<nNodes; i++) {
        xPositions.push_back(nodes[i].position);
    }

    Eigen::Map<DD> xOffsets(xPositions.data(), nNodes, 1);
    std::ofstream fileX("x_be.txt");
    std::ofstream fileZ("z_be.txt");

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