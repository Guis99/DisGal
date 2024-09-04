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

std::vector<DvD> TimeStep::solver_RK4(SpD &A, 
                            SpD &columnSpace, SpD &nullSpace, 
                            DvD &boundaryVals, 
                            DvD &initialCondition,
                            double timeStep,
                            int numTimeSteps) {
    std::vector<DvD> out; out.reserve(numTimeSteps);
    out.push_back(columnSpace * initialCondition);
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
        out.push_back(columnSpace * x);
    }

    return out;
}

std::vector<DvD> TimeStep::solver_GL1(SpD &A, 
                            SpD &columnSpace, SpD &nullSpace, 
                            DvD &boundaryVals, 
                            DvD &initialCondition,
                            double timeStep,
                            int numTimeSteps) {
    std::vector<DvD> out; out.reserve(numTimeSteps);
    out.push_back(columnSpace * initialCondition);
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
        out.push_back(columnSpace * x);
    }

    return out;
}

std::vector<DvD> TimeStep::solver_GL2(SpD &A, 
                            SpD &columnSpace, SpD &nullSpace, 
                            DvD &boundaryVals, 
                            DvD &initialCondition,
                            double timeStep,
                            int numTimeSteps) {
    std::vector<DvD> out; out.reserve(numTimeSteps);
    out.push_back(columnSpace * initialCondition);
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
        out.push_back(columnSpace * x);
    }

    return out;
}

std::vector<DvD> TimeStep::solver_GL2_nonlinear(SpD &A, 
                            SpD &columnSpace, SpD &nullSpace, 
                            DvD &boundaryVals, 
                            DvD &initialCondition,
                            double timeStep,
                            int numTimeSteps) {
    std::vector<DvD> out; out.reserve(numTimeSteps);
    out.push_back(columnSpace * initialCondition);
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
        out.push_back(columnSpace * x);
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
        } else {
            tripletListK.emplace_back(rightVal, rightVal, -1.0); // upwind flux, outward jump
            tripletListK.emplace_back(rightVal + 1, rightVal, 1.0); // upwind flux, inward jump
        }
    }

    // Declare and construct sparse matrix from triplets
    MassMatrix.setFromTriplets(tripletListM.begin(), tripletListM.end());
    StiffnessMatrix.setFromTriplets(tripletListK.begin(), tripletListK.end());
}

std::vector<DvD> ComputeTransientSolution(SpD &StiffnessMatrix, 
                                SpD &MassMatrix, SpD &columnSpace, 
                                SpD &nullSpace, DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps,
                                int integrator) {
    // Eliminate boundary rows and columns
    SpD K11 = columnSpace.transpose() * StiffnessMatrix * columnSpace;
    SpD M11 = columnSpace.transpose() * MassMatrix * columnSpace;

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


    // std::cout<<"inverse"<<std::endl<<combinedMatsInv<<std::endl;

    std::vector<DvD> out;
    out.push_back(columnSpace * initialCondition);
    DvD prevState = initialCondition;
    DvD x;

    switch (integrator) {
        case 0: // Forward Euler
            break;
        case 1: // Crank-Nicholson
            break;
        case 2: // RK4
            return TimeStep::solver_RK4(A, columnSpace, nullSpace, boundaryVals, initialCondition, timeStep, numTimeSteps);
            break;
        case 3: // GL1
            return TimeStep::solver_GL1(A, columnSpace, nullSpace, boundaryVals, initialCondition, timeStep, numTimeSteps);
            break;
        case 4: // GL2
            return TimeStep::solver_GL2(A, columnSpace, nullSpace, boundaryVals, initialCondition, timeStep, numTimeSteps);
            break;
    }

    // Time-stepping
    for (int i=1; i<numTimeSteps; i++) {
        x = combinedMatsInv * (M11 + timeStep * K11 / 2) * prevState;
        prevState = x;
        out.push_back(columnSpace * x);
    }
    return out;
}