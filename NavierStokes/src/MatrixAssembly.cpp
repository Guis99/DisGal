#include <iostream>
#include <algorithm>
#include "..\include\MatrixAssembly.hpp"

using namespace Solvers;

// typedef Eigen::SparseMatrix<double> SpD;
// typedef Eigen::MatrixXd DD;
// typedef Eigen::VectorXd DvD;

DD MatrixAssembly::GenerateQuadWeights(std::vector<double> &gpX, std::vector<double> &gpY, int numXNodes, int numYNodes, int numElemNodes) {
    // Generate quadrature weight matrices
    // Get integrals of basis functions
    std::vector<double> integX = Utils::integrateLagrange(gpX);
    std::vector<double> integY = Utils::integrateLagrange(gpY);

    Eigen::Map<DvD> integXMat(integX.data(),numXNodes);
    Eigen::Map<DvD> integYMat(integY.data(),numYNodes);

    DD weightMat(numElemNodes, numElemNodes);
    weightMat << Eigen::kroneckerProduct((DD)integXMat.asDiagonal(),(DD)integYMat.asDiagonal());

    return weightMat;
}

SpD MatrixAssembly::MassMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double c) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXNodes = xdeg + 1; int numYNodes = ydeg + 1;
    int numElemNodes = numXNodes * numYNodes;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);

    // get base elem matrix
    DD weightMat = MatrixAssembly::GenerateQuadWeights(gaussPointsX, gaussPointsY, numXNodes, numYNodes, numElemNodes);
    weightMat *= c;

    // Turn to vector because matrix is diagonal
    DvD weightVector = weightMat.diagonal();

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    // Integrate over all elements
    DvD localElemMat(numElemNodes);
    for (auto &elm : inputMesh.Elements) {
        double Lx = elm.getWidth(); double Ly = elm.getHeight(); // Jacobian factors
        // calculate local matrix
        localElemMat = weightVector*Lx*Ly/4;
        // Get nodes in element
        std::vector<int> nodesInElm = elm.Nodes;
        // Generate i,j,v triplets
        for (int i=0; i<numElemNodes; i++) {
            tripletList.emplace_back(nodesInElm[i], nodesInElm[i], localElemMat(i));
        }     
    }

    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

SpD MatrixAssembly::StiffnessMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double k) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXNodes = xdeg + 1; int numYNodes = ydeg + 1;
    int numElemNodes = numXNodes * numYNodes;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);

    // Generate derivative matrices
    std::vector<double> AxInitializer; std::vector<double> AyInitializer;
    AxInitializer.reserve(numXNodes*numXNodes); AyInitializer.reserve(numYNodes*numYNodes);

    double *AxInitIdx = AxInitializer.data(); 
    double *AyInitIdx = AyInitializer.data(); 

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numXNodes; k++) { 
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPointsX, gaussPointsX);
        std::copy(xPartials.begin(), xPartials.end(), AxInitIdx);
        AxInitIdx += numXNodes;
    }
    
    for (int k=0; k<numYNodes; k++) {
        std::vector<double> yPartials = Utils::numDeriv(.00001, k, gaussPointsY, gaussPointsY);
        std::copy(yPartials.begin(), yPartials.end(), AyInitIdx);
        AyInitIdx += numYNodes;
    }

    // map derivative values to matrix
    Eigen::Map<DD> Ax(AxInitializer.data(), numXNodes, numXNodes); 
    Eigen::Map<DD> Ay(AyInitializer.data(), numYNodes, numYNodes);

    // Generate quadrature weight matrices
    DD weightMat = MatrixAssembly::GenerateQuadWeights(gaussPointsX, gaussPointsY, numXNodes, numYNodes, numElemNodes);

    // Generate mass matrices
    DD Bx; Bx.setIdentity(numXNodes, numXNodes);
    DD By; By.setIdentity(numYNodes, numYNodes);

    DD coeffMat(numElemNodes, numElemNodes);
    coeffMat.setIdentity(); coeffMat *= k;
    
    // Get element-wise matrix intermediates
    DD combinedX(numElemNodes, numElemNodes);
    combinedX << Eigen::kroneckerProduct(By, Ax);
    DD combinedY(numElemNodes, numElemNodes);
    combinedY << Eigen::kroneckerProduct(Ay, Bx);

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    auto combineXT = (DD)combinedX.transpose();
    auto combineYT = (DD)combinedY.transpose();

    // Integrate over all elements
    DD localElemMat(numElemNodes, numElemNodes);
    for (auto &elm : inputMesh.Elements) {
        double Lx = elm.getWidth(); double Ly = elm.getHeight(); // Jacobian factors
        // calculate local matrix
        // auto int1 = (DD)(combineXT*coeffMat);
        // auto int2 = (DD)(combineYT*coeffMat);

        // auto int5 = (DD)(int1*weightMat);
        // auto int6 = (DD)(int2*weightMat);
        // std::cout<<"here9"<<std::endl;
        // auto int3 = (DD)(int5*combinedX*Ly/Lx);
        // auto int4 = (DD)(int6*combinedY*Lx/Ly);
        localElemMat = combinedX.transpose()*coeffMat*weightMat*combinedX*Ly/Lx +
                        combinedY.transpose()*coeffMat*weightMat*combinedY*Lx/Ly;
        // std::cout<<"here10"<<std::endl;
        // localElemMat = int3+int4;
        // std::cout<<"here11"<<std::endl;
        
        // Get nodes in element
        std::vector<int> nodesInElm = elm.Nodes;
        // Generate i,j,v triplets
        for (int j=0; j<numElemNodes; j++) {
            for (int i=0; i<numElemNodes; i++) {
                tripletList.emplace_back(nodesInElm[i],nodesInElm[j],localElemMat(i,j));
            }
        }
    }
    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

SpD MatrixAssembly::AssembleFVec(Meshing::BasicMesh::BasicMesh2D &inputMesh, double f) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXNodes = xdeg + 1; int numYNodes = ydeg + 1;
    int numElemNodes = numXNodes * numYNodes;
    int nNodes = inputMesh.nNodes(); int nElements = inputMesh.nElements();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);
    std::vector<double> gaussPointsY = Utils::genGaussPoints(ydeg);

    // Generate quadrature weight matrices
    DD weightMat = MatrixAssembly::GenerateQuadWeights(gaussPointsX, gaussPointsY, numXNodes, numYNodes, numElemNodes);

    // Turn weight mat int vector and mult. by source since diagonal
    DvD sourceVec = f * weightMat.diagonal();

    auto allNodesPos = inputMesh.allNodePos();
    auto startpoint = allNodesPos.data(); auto allocSize = allNodesPos.size();

    auto fEval = Utils::EvalSymbolicBC(startpoint, allocSize, "Specify source term");

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes);

    // Integrate over all elements
    DvD localElemMat(numElemNodes);
    for (auto &elm : inputMesh.Elements) {
        std::vector<int> nodesInElm = elm.Nodes;
        double Lx = elm.getWidth(); double Ly = elm.getHeight(); // Jacobian factors
        // calculate local matrix
        std::vector<double> collectSourceVals; collectSourceVals.reserve(numElemNodes);

        for (int i : nodesInElm) {
            collectSourceVals.push_back(fEval[i]);
        }
        Eigen::Map<DvD> sourceVector(collectSourceVals.data(), numElemNodes, 1);

        localElemMat = weightMat*sourceVector*Lx*Ly/4;
        // Get nodes in element
        
        // Generate i,j,v triplets
        for (int i=0; i<numElemNodes; i++) {
            tripletList.emplace_back(nodesInElm[i], 0, localElemMat(i));
        }  
    }

    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,1);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

void MatrixAssembly::GetExtensionMatrices(Meshing::BasicMesh::BasicMesh2D &inputMesh,
                                        std::vector<int> &boundaryNodes, 
                                        std::vector<int> &freeNodes,
                                        SpD &nullSpace,
                                        SpD &columnSpace) {
    int nNodes = inputMesh.nNodes();

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

DvD MatrixAssembly::EvalDirichletBoundaryCond(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &boundaryNodes) {
    int xdeg = inputMesh.xdeg; int ydeg = inputMesh.ydeg;
    int numXElems = inputMesh.xOffsets.size()-1; int numYElems = inputMesh.yOffsets.size()-1;
    int xWidth = xdeg*numXElems; int yWidth = ydeg*numYElems;
    
    int numBoundaryNodes = boundaryNodes.size();

    std::vector<std::array<double,2>> boundaryNodePos = inputMesh.posOfNodes(boundaryNodes);
    // Boundary nodes are given in clockwise order, not column-major order
    std::vector<double> boundaryNodeValues; 
    boundaryNodeValues.resize(numBoundaryNodes);

    std::vector<double> boundaryCalc;
    std::array<double,2> *currPointer = boundaryNodePos.data();
    auto *currNodeValPointer = boundaryNodeValues.data();
    int ptrIncr;
    std::string prompt = "Specify boundary condition";

    for (int i=0; i<4; i++) {
        // Alternate between x and y-iteration
        ptrIncr = i % 2 == 0 ? xWidth : yWidth;
        // Take bcFunc and evaluate it
        boundaryCalc = Utils::EvalSymbolicBC(currPointer, ptrIncr, prompt);
        std::copy(boundaryCalc.begin(), boundaryCalc.end(), currNodeValPointer); 
        currPointer += ptrIncr; currNodeValPointer += ptrIncr;
    }

    auto RmOrder = boundaryNodes;
    std::sort(RmOrder.begin(), RmOrder.end());

    auto BcValsSorted = Utils::ReshuffleNodeVals(RmOrder, boundaryNodes, boundaryNodeValues);

    Eigen::Map<DvD> boundaryNodeValuesVec(BcValsSorted.data(), numBoundaryNodes, 1);
    return (DvD)boundaryNodeValuesVec;
}

DvD MatrixAssembly::ComputeSolutionStationary(SpD &StiffnessMatrix, SpD &fVec, SpD &columnSpace, SpD &nullSpace, DvD &boundaryVals) {
    // Eliminate rows and columns corr. to boundary nodes
    SpD A11 = columnSpace.transpose() * StiffnessMatrix * columnSpace;
    // Eliminate boundary rows and free columns
    SpD A12 = columnSpace.transpose() * StiffnessMatrix * nullSpace;
    // Eliminate boundary rows
    SpD F11 = columnSpace.transpose() * fVec;

    Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int> > LuSolver;    
    LuSolver.analyzePattern(A11);
    LuSolver.factorize(A11);

    DvD x = LuSolver.solve(F11 - A12 * boundaryVals);
    x = columnSpace * x + nullSpace * boundaryVals;
    return x;
}

std::vector<DvD> MatrixAssembly::ComputeSolutionTimeDependent1stOrder(SpD &StiffnessMatrix, 
                                SpD &MassMatrix, 
                                SpD &fVec, SpD &columnSpace, 
                                SpD &nullSpace, DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps) {
    std::cout<<"here1"<<std::endl;
    SpD K11 = columnSpace.transpose() * StiffnessMatrix * columnSpace;
    std::cout<<"here2"<<std::endl;
    SpD M11 = columnSpace.transpose() * MassMatrix * columnSpace;
    std::cout<<"here3"<<std::endl;
    // Eliminate boundary rows and free columns
    SpD K12 = columnSpace.transpose() * StiffnessMatrix * nullSpace;
    // Eliminate boundary rows
    SpD F11 = columnSpace.transpose() * fVec;
    std::cout<<"here4"<<std::endl;

    SpD combinedMats = timeStep * K11 + M11;
    std::cout<<"here5"<<std::endl;
    int system_size = combinedMats.rows();
    SpD dummyId(system_size, system_size); 
    std::vector<Eigen::Triplet<double>> tripletListID;
    std::cout<<"here6"<<std::endl;

    tripletListID.reserve(combinedMats.rows());

    for (int i=0; i<combinedMats.rows(); i++) {
        std::cout<<"iter "<<i<<std::endl;
        tripletListID.emplace_back(i, i, 1.0);
    }
    std::cout<<"here100"<<std::endl;
    
    dummyId.setFromTriplets(tripletListID.begin(), tripletListID.end());
    std::cout<<"here7"<<std::endl;
    Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int> > LuSolver;    
    LuSolver.analyzePattern(combinedMats);
    LuSolver.factorize(combinedMats);
    std::cout<<"here8"<<std::endl;
    
    SpD combinedMatsInv = LuSolver.solve(dummyId);

    std::vector<DvD> out;
    out.reserve(numTimeSteps);
    std::cout<<"here9"<<std::endl;
    out.push_back(columnSpace * initialCondition + nullSpace * boundaryVals);
    std::cout<<"here10"<<std::endl;
    DvD prevState = initialCondition;
    DvD x;

    // Time-stepping
    for (int i=1; i<numTimeSteps; i++) {
        std::cout<<"iter "<<i<<std::endl;
        std::cout<<"here11"<<std::endl;
        x = combinedMatsInv * (M11 * prevState + timeStep * F11 - K12 * boundaryVals);
        std::cout<<"here12"<<std::endl;
        prevState = x;
        std::cout<<"here13"<<std::endl;
        out.push_back(columnSpace * x + nullSpace * boundaryVals);
    }
    return out;
}