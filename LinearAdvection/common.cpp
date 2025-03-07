#include "common.hpp"

BasicMesh1D CreateMesh(int interpolantOrder, int numElems, 
                        double elemWidth) {
    BasicMesh1D mesh;

    // Populate basic mesh props
    mesh.interpolantOrder = interpolantOrder;
    mesh.numElems = numElems;
    mesh.elemWidth = elemWidth;

    // Generate nodes for basis func
    auto elemIntervals = Utils::genGaussPoints(interpolantOrder);
    mesh.spacing = elemIntervals;

    std::vector<Node> nodes;
    std::vector<Element> elements;

    // Populate node list
    Node node;
    int NID = 0;
    nodes.reserve(numElems * interpolantOrder + 1);
    for (int i=0; i<numElems; i++) { // iter over elements
        double offset = elemWidth * i;
        for (int j=0; j<interpolantOrder + 1; j++) { // iter from element left bound
            // Calculate pos and push
            node.ID = NID++;
            node.position = Utils::TransformPoint(elemIntervals[j], offset, offset+elemWidth);
            nodes.push_back(node);
        }
    }

    std::vector<int> nodeIds; std::vector<int> dofs;
    nodeIds.reserve(interpolantOrder + 1);
    dofs.reserve(interpolantOrder + 1);
    int dofId = 0;
    Element elem; int EID = 0; elements.reserve(numElems);
    for (int i=0; i<numElems; i++) {
        elem.ID = EID++;
        for (int j=0; j<interpolantOrder + 1; j++) {
            nodeIds.push_back(i*(interpolantOrder+1) + j);
            dofs.push_back(dofId++);
        }
        elem.dofs = dofs;
        elements.push_back(elem);
        dofs.clear();
    }

    mesh.nodes = nodes;
    mesh.elements = elements;

    return mesh;
}

std::vector<double> GetPosOfElemNOdes(BasicMesh1D& mesh, int elem) {
    std::vector<double> out;
    auto elemNodes = mesh.elements[elem].dofs;
    out.reserve(mesh.interpolantOrder + 1);
    for (int i=0; i<mesh.interpolantOrder + 1; i++) {
        out.push_back(mesh.nodes[elemNodes[i]].position);
    }
    return out;
}

DD GenerateQuadWeights(std::vector<double> &gpX, int numXNodes) {
    // Generate quadrature weight matrices
    // Get integrals of basis functions
    std::vector<double> integX = Utils::integrateLagrange(gpX);

    Eigen::Map<DvD> integXMat(integX.data(),numXNodes);

    DD weightMat = integXMat.asDiagonal();

    return weightMat;
}

void GetExtensionMatrices(BasicMesh1D &inputMesh,
                                        std::vector<int> &boundaryNodes, 
                                        std::vector<int> &freeNodes,
                                        SpD &nullSpace,
                                        SpD &columnSpace) {
    int nNodes = (inputMesh.nodes).size();

    std::sort(boundaryNodes.begin(), boundaryNodes.end());

    std::vector<Eigen::Triplet<double>> tripletListNS;
    std::vector<Eigen::Triplet<double>> tripletListCS;
    tripletListNS.reserve(boundaryNodes.size());
    tripletListCS.reserve(freeNodes.size());

    for (int i=0; i<boundaryNodes.size(); i++) {
        tripletListNS.emplace_back(boundaryNodes[i], i, 1.0);
    }

    for (int i=0; i<freeNodes.size(); i++) {
        tripletListCS.emplace_back(freeNodes[i], i, 1.0);
    }
    
    nullSpace.setFromTriplets(tripletListNS.begin(), tripletListNS.end());
    columnSpace.setFromTriplets(tripletListCS.begin(), tripletListCS.end());
}

DvD getICVec(const BasicMesh1D& mesh, int nNodes, char* baseline, double cutoff, std::string IC) {
    int stopIdx = 0;
    std::vector<double> x; x.reserve(nNodes);

    DEBUG_PRINT("X POS LOOP");
    while (mesh.nodes[stopIdx].position <= cutoff && stopIdx < nNodes) { // if cutoff is greater than last position we assume no discon
        // DEBUG_PRINT(mesh.nodes[stopIdx].position);
        // DEBUG_PRINT(stopIdx);
        x.push_back(mesh.nodes[stopIdx++].position); 
        // DEBUG_PRINT("-------");
    }

    DEBUG_PRINT("ALL x");
    for (auto vecElem : x) {
        DEBUG_PRINT(vecElem);
    }
    // x is NOT the initial values, but the positions of each node
    std::string ICMod = "(" + IC + ")+";
    ICMod += baseline;
    DEBUG_PRINT("ic string: ", ICMod);
    std::vector<double> ICVec = Utils::EvalSymbolicBC1D(x.data(), stopIdx, ICMod);

    for (auto vecElem : ICVec) {
        DEBUG_PRINT(vecElem);
    }

    ICVec.reserve(nNodes);

    double baselineD = std::stod(baseline);
    for (stopIdx; stopIdx<nNodes; stopIdx++) {
        ICVec.push_back(baselineD);
    }

    DEBUG_PRINT("extended");
    for (auto vecElem : ICVec) {
        DEBUG_PRINT(vecElem);
    }

    Eigen::Map<DvD> InitialCondition(ICVec.data(), nNodes, 1);
    DvD out = ((DvD)InitialCondition);

    return out;
}

std::vector<DvD> TimeStep::solver_RK4(SpD &A, 
                            SpD &columnSpace, SpD &nullSpace, 
                            DvD &initialCondition,
                            double timeStep,
                            int numTimeSteps) {
    std::vector<DvD> out; out.reserve(numTimeSteps);
    out.push_back(initialCondition);
    DvD prevState = initialCondition;
    DvD x;
    DvD k1; DvD k2; DvD k3; DvD k4;
    
    for (int i=1; i<numTimeSteps; i++) {
        k1 = A * prevState;
        k2 = A * (prevState + .5 * timeStep * k1);
        k3 = A * (prevState + .5 * timeStep * k2);
        k4 = A * (prevState + timeStep * k3);

        x = prevState + (timeStep / 6) * (k1 + 2*k2 + 2*k3 + k4);
        prevState = x;
        out.push_back(x);
    }

    return out;
}

std::vector<DvD> TimeStep::solver_GL1(SpD &A, 
                            SpD &columnSpace, SpD &nullSpace, 
                            DvD &initialCondition,
                            double timeStep,
                            int numTimeSteps) {
    std::vector<DvD> out; out.reserve(numTimeSteps);
    out.push_back(initialCondition);
    DvD prevState = initialCondition;
    DvD x;
    DvD k1; DvD k2; DvD k3; DvD k4;
    
    for (int i=1; i<numTimeSteps; i++) {
        k1 = A * prevState;
        k2 = A * (prevState + .5 * timeStep * k1);
        k3 = A * (prevState + .5 * timeStep * k2);
        k4 = A * (prevState + timeStep * k3);

        x = prevState + (timeStep / 6) * (k1 + 2*k2 + 2*k3 + k4);
        prevState = x;
        out.push_back(x);
    }
    return out;
}

std::vector<DvD> TimeStep::solver_GL2(SpD &A, 
                            SpD &columnSpace, SpD &nullSpace, 
                            DvD &initialCondition,
                            double timeStep,
                            int numTimeSteps) {
    std::vector<DvD> out; out.reserve(numTimeSteps);
    out.push_back(initialCondition);
    DEBUG_PRINT("ic gl2 begin size: ",initialCondition.rows()," x ", initialCondition.cols());
    DvD prevState = initialCondition;
    DvD x;

    // Generate identity matrix
    int n = A.rows();
    SpD I_n(n, n); 
    SpD I_2n(2*n, 2*n);
    std::vector<Eigen::Triplet<double>> tripletListID;
    std::vector<Eigen::Triplet<double>> tripletListIDTwo;
    std::cout<<"time1"<<std::endl;

    tripletListID.reserve(n);

    for (int i=0; i<n; i++) {
        tripletListID.emplace_back(i, i, 1.0);
    }
    I_n.setFromTriplets(tripletListID.begin(), tripletListID.end());

    for (int i=0; i<2*n; i++) {
        tripletListIDTwo.emplace_back(i, i, 1.0);
    }
    I_2n.setFromTriplets(tripletListIDTwo.begin(), tripletListIDTwo.end());

    // Constants
    double SQRT3 = std::sqrt(3);

    double a11 = .25; double a12 = .25 - SQRT3/6;
    double a21 = .25 + SQRT3/6; double a22 = .25;

    double c1 = .5 - SQRT3/6; double c2 = .5 + SQRT3/6;
    double b1 = .5; double b2 = .5;

    SpD blockRow(n, 2*n);
    SpD blockA(2*n, 2*n);
    SpD blockCol(2*n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (I_n.coeff(j, i) != 0) {
                blockRow.insert(j, i) = b1 * I_n.coeff(j, i);
            }
        }
    }

    // Top-right block: b2 * I_n
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (I_n.coeff(j, i) != 0) {
                blockRow.insert(j, n + i) = b2 * I_n.coeff(j, i);
            }
        }
    }

    // Initialize blockA (2n x 2n)
    blockA.resize(2 * n, 2 * n);
    blockA.setZero(); // Set all entries to zero

    // Populate blockA
    // Top-left block: I_n - timeStep * a11 * A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double value = I_n.coeff(j, i) - timeStep * a11 * A.coeff(j, i);
            if (value != 0) {
                blockA.insert(j, i) = value;
            }
        }
    }

    // Top-right block: -timeStep * a12 * A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double value = -timeStep * a12 * A.coeff(j, i);
            if (value != 0) {
                blockA.insert(j, n + i) = value;
            }
        }
    }

    // Bottom-left block: I_n - timeStep * a21 * A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double value = -timeStep * a21 * A.coeff(j, i);
            if (value != 0) {
                blockA.insert(n + j, i) = value;
            }
        }
    }

    // Bottom-right block: -timeStep * a22 * A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double value = I_n.coeff(j, i) - timeStep * a22 * A.coeff(j, i);
            if (value != 0) {
                blockA.insert(n + j, n + i) = value;
            }
        }
    }

    // Initialize blockCol (2n x n)
    blockCol.resize(2 * n, n);
    blockCol.setZero(); // Set all entries to zero

    // Populate blockCol
    // Top block: A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (A.coeff(j, i) != 0) {
                blockCol.insert(j, i) = A.coeff(j, i);
                blockCol.insert(n + j, i) = A.coeff(j, i);
            }
        }
    }

    // Finalize the matrices
    blockRow.makeCompressed();
    blockA.makeCompressed();
    blockCol.makeCompressed();
    
    Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolverMM;    
    LuSolverMM.analyzePattern(blockA);
    LuSolverMM.factorize(blockA);

    SpD blockA_Inv = LuSolverMM.solve(I_2n);

    SpD fPropOp = I_n + timeStep * blockRow * blockA_Inv * blockCol;

    DEBUG_PRINT("fprop size: ",fPropOp.rows()," x ", fPropOp.cols());
    DEBUG_PRINT("pstate size: ",prevState.rows()," x ", prevState.cols());
    
    for (int i=1; i<numTimeSteps; i++) {
        x = fPropOp * prevState;
        prevState = x;
        out.push_back(x);
    }
    return out;
}

std::vector<DvD> TimeStep::solver_GL2_nonlinear(SpD &A, 
                            SpD &columnSpace, SpD &nullSpace, 
                            DvD &initialCondition,
                            double timeStep,
                            int numTimeSteps) {
    std::vector<DvD> out; out.reserve(numTimeSteps);
    out.push_back(initialCondition);
    DvD prevState = initialCondition;
    DvD x;
    DvD k1; DvD k2; DvD k3; DvD k4;

    // Generate identity matrix
    int system_size = A.rows();
    SpD dummyId(system_size, system_size); 
    std::vector<Eigen::Triplet<double>> tripletListID;

    tripletListID.reserve(system_size);

    for (int i=0; i<system_size; i++) {
        tripletListID.emplace_back(i, i, 1.0);
    }
    dummyId.setFromTriplets(tripletListID.begin(), tripletListID.end());

    // Constants
    double SQRT3 = std::sqrt(3);

    double a11 = .25; double a12 = .25 - SQRT3/6;
    double a21 = .25 + SQRT3/6; double a22 = .25;

    double c1 = .5 - SQRT3/6; double c2 = .5 + SQRT3/6;
    double b1 = .5; double b2 = .5;
    
    for (int i=1; i<numTimeSteps; i++) {
        // rhs = A * prevState;
    

        x = prevState + (timeStep / 6) * (k1 + 2*k2 + 2*k3 + k4);
        prevState = x;
        out.push_back(x);
    }

    return out;
}

void AssembleSystemMatrices(BasicMesh1D& mesh, SpD &MassMatrix, SpD &StiffnessMatrix) {
    int xdeg = mesh.interpolantOrder;
    int numElemNodes = xdeg + 1;
    int nNodes = mesh.nodes.size(); int nElements = mesh.elements.size();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);

    // Generate derivative matrices
    std::vector<double> AInitializer;
    AInitializer.reserve(numElemNodes * numElemNodes);

    double *AxInitIdx = AInitializer.data(); 

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numElemNodes; k++) { 
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPointsX, gaussPointsX);
        std::copy(xPartials.begin(), xPartials.end(), AxInitIdx);
        AxInitIdx += numElemNodes;
    }
    Eigen::Map<DD> Ax(AInitializer.data(), numElemNodes, numElemNodes); 

    // Generate quadrature weight matrices
    DD weightMat = GenerateQuadWeights(gaussPointsX, numElemNodes);

    std::vector<Eigen::Triplet<double>> tripletListM;
    std::vector<Eigen::Triplet<double>> tripletListK;
    tripletListM.reserve(nElements * numElemNodes);
    tripletListK.reserve(nElements * numElemNodes * numElemNodes);

    DvD massVec = weightMat.diagonal();
    DvD localElemVecMass(numElemNodes);
    DD localElemMatK(numElemNodes, numElemNodes);

    // Add contributions from basis functions (stiffness matrix is block-diagonal)
    for (auto &elm : mesh.elements) {
        double Lx = mesh.elemWidth; // Jacobian factors
        // calculate local matrix
        localElemVecMass = massVec*Lx/2;
        localElemMatK = Ax.transpose() * weightMat;
        
        // Get nodes in element
        std::vector<int> dofsInElem = elm.dofs;
        // Generate i,j,v triplets
        for (int j=0; j<numElemNodes; j++) {
            tripletListM.emplace_back(dofsInElem[j], dofsInElem[j], localElemVecMass(j));
            for (int i=0; i<numElemNodes; i++) {
                tripletListK.emplace_back(dofsInElem[i],dofsInElem[j],localElemMatK(i,j));
            }
        }
    }

    // Add contributions from flux terms
    for (auto &elm : mesh.elements) {
        int rightVal = (elm.dofs).back(); 
        if (elm.ID == nElements - 1) {
            tripletListK.emplace_back(rightVal, rightVal, -1.0); // right boundary
            tripletListK.emplace_back(0, rightVal, 1.0); // periodic boundary
        } else {
            tripletListK.emplace_back(rightVal, rightVal, -1.0); // upwind flux, outward jump
            tripletListK.emplace_back(rightVal + 1, rightVal, 1.0); // upwind flux, inward jump
        }
    }

    // Declare and construct sparse matrix from triplets
    MassMatrix.setFromTriplets(tripletListM.begin(), tripletListM.end());
    StiffnessMatrix.setFromTriplets(tripletListK.begin(), tripletListK.end());
}

std::vector<DvD> ComputeTransientSolution(SpD &K11, 
                                SpD &M11, SpD &columnSpace, 
                                SpD &nullSpace,  
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps,
                                int integrator) {
    // Eliminate boundary rows and columns

    SpD combinedMats = M11 - timeStep * K11 / 2;
    int system_size = combinedMats.rows();
    SpD dummyId(system_size, system_size); 
    std::vector<Eigen::Triplet<double>> tripletListID;

    tripletListID.reserve(combinedMats.rows());

    for (int i=0; i<combinedMats.rows(); i++) {
        tripletListID.emplace_back(i, i, 1.0);
    }
    
    // Crank-Nicholson-specifc
    dummyId.setFromTriplets(tripletListID.begin(), tripletListID.end());
    Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolver;    
    LuSolver.analyzePattern(combinedMats);
    LuSolver.factorize(combinedMats);
    SpD combinedMatsInv = LuSolver.solve(dummyId);

    // Invert diagonal mass matrix
    Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolverMM;    
    LuSolverMM.analyzePattern(M11);
    LuSolverMM.factorize(M11);

    SpD A = LuSolverMM.solve(dummyId) * K11;

    // std::cout<<"A: "<<A<<std::endl;

    // std::cout<<"inverse"<<std::endl<<combinedMatsInv<<std::endl;

    std::vector<DvD> out;
    out.push_back(initialCondition);
    DvD prevState = initialCondition;
    DvD x;
    DEBUG_PRINT("ic trasol begin size: ",initialCondition.rows()," x ", initialCondition.cols());


    switch (integrator) {
        case 0: // Forward Euler
            break;
        case 1: // Crank-Nicholson
            break;
        case 2: // RK4
            return TimeStep::solver_RK4(A, columnSpace, nullSpace, initialCondition, timeStep, numTimeSteps);
            break;
        case 3: // GL1
            return TimeStep::solver_GL1(A, columnSpace, nullSpace, initialCondition, timeStep, numTimeSteps);
            break;
        case 4: // GL2
            return TimeStep::solver_GL2(A, columnSpace, nullSpace, initialCondition, timeStep, numTimeSteps);
            break;
    }

    // Time-stepping
    for (int i=1; i<numTimeSteps; i++) {
        x = combinedMatsInv * (M11 + timeStep * K11 / 2) * prevState;
        prevState = x;
        out.push_back(x);
    }
    return out;
}