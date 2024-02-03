#include "..\include\MatrixAssembly.hpp"

DD GenerateQuadWeights(std::vector<double> &gpX, std::vector<double> &gpY, int numXNodes, int numYNodes, int numElemNodes) {
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


SpD StiffnessMatrix(QTM::QuadTreeMesh mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);

    // Generate derivative matrix
    std::vector<double> AInitializer; 
    AInitializer.reserve(numElemNodes * numElemNodes);

    double *AInitIdx = AInitializer.data(); 

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPoints, gaussPoints);
        std::copy(xPartials.begin(), xPartials.end(), AInitIdx);
        AInitIdx += numNodes;
    }
    
    // map derivative values to matrix
    Eigen::Map<DD> A(AInitializer.data(), numNodes, numNodes); 

    // Generate quadrature weight matrices
    DD weightMat = GenerateQuadWeights(gaussPoints, gaussPoints, numNodes, numNodes, numElemNodes);

    // Generate mass matrices
    DD B; B.setIdentity(numNodes, numNodes);

    DD coeffMat(numElemNodes, numElemNodes);
    coeffMat.setIdentity(); coeffMat *= k;
    
    // Get element-wise matrix intermediates
    DD combinedX(numElemNodes, numElemNodes);
    combinedX << Eigen::kroneckerProduct(B, A);
    DD combinedY(numElemNodes, numElemNodes);
    combinedY << Eigen::kroneckerProduct(A, B);

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    auto combineXT = (DD)combinedX.transpose();
    auto combineYT = (DD)combinedY.transpose();

    // Integrate over all elements
    DD localElemMat(numElemNodes, numElemNodes);
    auto leaves = mesh.GetAllCells();
    for (auto &elm : leaves) {
        // calculate local matrix
        localElemMat = combineXT*coeffMat*weightMat*combinedX +
                        combineYT*coeffMat*weightMat*combinedY;
 
        // Get nodes in element
        auto nodeBound = elm->nodes[0];
        // std::vector<int> nodesInElm = elm->Nodes;
        // Generate i,j,v triplets
        for (int j=0; j<numElemNodes; j++) {
            for (int i=0; i<numElemNodes; i++) {
                tripletList.emplace_back(nodeBound+i,nodeBound+j,localElemMat(i,j));
            }
        }
    }
    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

SpD PenaltyMatrix(QTM::QuadTreeMesh mesh, double k, double alpha) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;

    // basis func to node mapping
    DD B; B.setIdentity(numElemNodes, numElemNodes);
    DD Bs; Bs.setIdentity(2*numNodes-1, 2*numNodes-1);

    // get basis func vals for split cell quad
    std::vector<double> BhInitializer; 
    BhInitializer.reserve(numNodes * numNodes);

    // get unit vecs
    DvD topVec(numNodes); topVec(0) = 1;
    DvD bottomVec(numNodes); bottomVec(deg) = 1;

    for (int i=0; i<numNodes-1; i++) {
        topVec(i+1) = 0;
        bottomVec(i) = 0;
    }
    
    // Generate values for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<double> xVals = Utils::evalLagrangeInterp(k, mesh.halfGaussPoints, mesh.gaussPoints);
        BhInitializer.insert(BhInitializer.end(), xVals.begin(), xVals.end());
    }
    
    // map basis values to matrix
    Eigen::Map<DD> Bh(BhInitializer.data(), 2*numNodes-1, numNodes); // eval points are traversed first
    Bh.transposeInPlace();

    std::array<DD,4> splitCellVals;

    DD splitCellPlaceholder(numElemNodes, mesh.halfGaussPoints.size());
    splitCellPlaceholder << Eigen::kroneckerProduct(bottomVec, Bh); splitCellVals[0] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(Bh, bottomVec); splitCellVals[1] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(topVec, Bh); splitCellVals[2] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(Bh, topVec); splitCellVals[3] = splitCellPlaceholder;
    
    // get quad weights in 1D
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
    std::vector<double> integX = Utils::integrateLagrange(gaussPoints);

    // Convert 1D quad weights to diag matrix
    Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    DD quadWeights1D = (DD)integXMat.asDiagonal();

    auto leaves = mesh.GetAllCells();
    std::vector<QTM::Direction> directions = {QTM::Direction::N, QTM::Direction::E, 
                                                QTM::Direction::S, QTM::Direction::W};

    std::vector<std::vector<int>> localNodes = {mesh.GetLocalBoundaryNodes(QTM::Direction::N),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::E),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::S),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::W)};
    // Signs for element edges
    int edgeSign[4] = {1,1,-1-1};
    int cellSign; int neighborSign;
    double a; // penalty parameters

    std::vector<int> boundaryNodes;
    std::vector<int> neighborNodes;
    std::vector<int> elemNodes;
    std::vector<std::shared_ptr<QTM::Cell>> neighbors;
    std::vector<Eigen::Triplet<double>> tripletList; tripletList.reserve(nNodes);
    for (auto &elm : leaves) {
        // get neighbors
        for (auto dir : directions) {
            cellSign = edgeSign[dir]; neighborSign = -cellSign;
            neighbors = mesh.GetCellNeighbors(dir, elm->CID);
            elemNodes = mesh.GetGlobalBoundaryNodes(dir, elm->CID);
            for (auto neighbor : neighbors) { 
                if (neighbor == nullptr) { // case curr boundary is exterior
                    
                } 
                else if (elm->level < neighbor->level) {

                } else if (neighbor->CID > elm->CID || elm->level > neighbor->level) { // wrong side of boundary
                    continue;
                } else { // neighbor exists 
                    neighborNodes = mesh.GetGlobalBoundaryNodes((QTM::Direction)((dir+2)%4), neighbor->CID);
                    // calculate penalty param
                    a = alpha/std::min(elm->width, neighbor->width)/2;
                    // calculate jump matrix
                    DD topJump = cellSign*B(Eigen::all, localNodes[dir]);
                    DD bottomJump = neighborSign*B(Eigen::all, localNodes[(QTM::Direction)((dir+2)%4)]);
                    DD jumpMatrix(2*numElemNodes, numNodes);
                    jumpMatrix << topJump, bottomJump;

                    boundaryNodes.reserve(2*numNodes);
                    boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                    boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                    auto jumpMatrixT = (DD)jumpMatrix.transpose();
                    DD localElemMat = (DD)(jumpMatrix * elm->width*quadWeights1D * jumpMatrixT);

                    for (int j=0; j<boundaryNodes.size(); j++) {
                        for (int i=0; i<boundaryNodes.size(); i++) {
                            if (localElemMat(i,j) != 0) {
                                tripletList.emplace_back(boundaryNodes[i],boundaryNodes[j],localElemMat(i,j));
                            }
                        }
                    }    
                }
                boundaryNodes.clear();
            }
        }
    }
    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

SpD FluxMatrix(QTM::QuadTreeMesh mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);

    // get quad weights in 1D
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
    std::vector<double> integX = Utils::integrateLagrange(gaussPoints);

    // Convert 1D quad weights to diag matrix
    Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    DD quadWeights1D = (DD)integXMat.asDiagonal();

    // basis func to node mapping
    DD B; B.setIdentity(numElemNodes, numElemNodes);
    DD Bs; Bs.setIdentity(numNodes, numNodes);

    // Generate derivative matrix
    std::vector<double> AInitializer; 
    AInitializer.reserve(numElemNodes * numElemNodes);

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPoints, gaussPoints);
        AInitializer.insert(AInitializer.end(), xPartials.begin(), xPartials.end());
    }
    
    // map derivative values to matrix
    Eigen::Map<DD> A(AInitializer.data(), numNodes, numNodes); 

    // Get element-wise matrix intermediates
    DD combinedX(numElemNodes, numElemNodes);
    combinedX << Eigen::kroneckerProduct(Bs, A);
    DD combinedY(numElemNodes, numElemNodes);
    combinedY << Eigen::kroneckerProduct(A, Bs);

    // Signs for element edges
    int edgeSign[4] = {1,1,-1-1};
    int cellSign; int neighborSign;

    std::vector<int> boundaryNodes;
    std::vector<int> neighborNodes;
    std::vector<int> elemNodes;
    std::vector<std::shared_ptr<QTM::Cell>> neighbors;
    std::vector<Eigen::Triplet<double>> tripletList; tripletList.reserve(nNodes);
    auto leaves = mesh.GetAllCells();
    std::vector<QTM::Direction> directions = {QTM::Direction::N, QTM::Direction::E, 
                                                QTM::Direction::S, QTM::Direction::W};

    std::vector<std::vector<int>> localNodes = {mesh.GetLocalBoundaryNodes(QTM::Direction::N),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::E),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::S),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::W)};
    for (auto &elm : leaves) {
        // get neighbors
        for (auto dir : directions) {
            cellSign = edgeSign[dir]; neighborSign = -cellSign;
            neighbors = mesh.GetCellNeighbors(dir, elm->CID);
            elemNodes = mesh.GetGlobalBoundaryNodes(dir, elm->CID);
            for (auto neighbor : neighbors) { 
                if (neighbor == nullptr) { // case curr boundary is exterior
                    continue;
                } else if (neighbor->CID > elm->CID) { // wrong side of boundary
                    
                } else { // neighbor exists 
                    neighborNodes = mesh.GetGlobalBoundaryNodes((QTM::Direction)((dir+2)%4), neighbor->CID);
                    // calculate jump matrix
                    DD topJump = cellSign*B(Eigen::all, localNodes[dir]);
                    DD bottomJump = neighborSign*B(Eigen::all, localNodes[(QTM::Direction)((dir+2)%4)]);
                    DD jumpMatrix(2*numElemNodes, elemNodes.size());
                    jumpMatrix << topJump, bottomJump;

                    // calculate flux matrix
                    DD topGradX = cellSign*combinedX(localNodes[dir], Eigen::all);
                    DD bottomGradX = neighborSign*combinedX(localNodes[(QTM::Direction)((dir+2)%4)], Eigen::all );

                    DD topGradY = cellSign*combinedY(localNodes[dir], Eigen::all);
                    DD bottomGradY = neighborSign*combinedY(localNodes[(QTM::Direction)((dir+2)%4)], Eigen::all);

                    DD fluxMatrixX(elemNodes.size(), 2*numElemNodes);
                    DD fluxMatrixY(elemNodes.size(), 2*numElemNodes);

                    fluxMatrixX << topGradX, bottomGradX;
                    fluxMatrixY << topGradY, bottomGradY;

                    boundaryNodes.reserve(elemNodes.size() + neighborNodes.size());
                    boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                    boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                    // assemble local matrix
                    DD localElemMat = (DD)(jumpMatrix * elm->width*quadWeights1D * fluxMatrixX + 
                                        jumpMatrix * elm->width*quadWeights1D * fluxMatrixY);

                    for (int j=0; j<boundaryNodes.size(); j++) {
                        for (int i=0; i<boundaryNodes.size(); i++) {
                            if (localElemMat(i,j) != 0) {
                            tripletList.emplace_back(boundaryNodes[i],boundaryNodes[j],localElemMat(i,j));
                            }
                        }
                    }
                }
                boundaryNodes.clear();
            }
        }
    }
    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

SpD AssembleFVec(QTM::QuadTreeMesh mesh, double f, std::string evalStr) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);

    // Generate quadrature weight matrices
    DD weightMat = GenerateQuadWeights(gaussPoints, gaussPoints, numNodes, numNodes, numElemNodes);

    // Turn weight mat int vector and mult. by source since diagonal
    // DvD sourceVec = f * weightMat.diagonal();

    auto allNodesPos = mesh.AllNodePos();
    auto startpoint = allNodesPos.data(); auto allocSize = allNodesPos.size();
    auto fEval = Utils::EvalSymbolicBC(startpoint, allocSize, evalStr);

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes);

    // Integrate over all elements
    DvD localElemMat(numElemNodes);
    auto leaves = mesh.GetAllCells();
    for (auto &elm : leaves) {
        auto nodesInElm = elm->nodes;
        double Lx = elm->width; // Jacobian factors
        // calculate local matrix
        std::vector<double> collectSourceVals; collectSourceVals.reserve(numElemNodes);

        for (int i=nodesInElm[0]; i<=nodesInElm[1]; i++) {
            collectSourceVals.push_back(fEval[i]);
        }
        Eigen::Map<DvD> sourceVector(collectSourceVals.data(), numElemNodes, 1);

        localElemMat = weightMat*sourceVector*Lx*Lx/4;
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

double ComputeResidual() {

}

std::vector<std::shared_ptr<QTM::Cell>> TestResiduals(DvD &solution, QTM::QuadTreeMesh& mesh, double residualLimit) {
    std::vector<std::shared_ptr<QTM::Cell>> out; out.reserve(mesh.leaves.size());
    for (auto leaf : mesh.leaves) {
        double res = ComputeResidual();
        if (leaf->level < 5 && ComputeResidual() > residualLimit) {
             out.push_back(leaf);
        }
    }

    return out;
}

DvD EvalDirichletBoundaryCond(QTM::QuadTreeMesh &inputMesh, std::vector<int> &boundaryNodes, std::vector<std::string> &strs) {
    int deg = inputMesh.deg;
    int numXElems = inputMesh.xOffsets.size()-1; int numYElems = inputMesh.yOffsets.size()-1;
    int xWidth = deg*numXElems;
    
    int numBoundaryNodes = boundaryNodes.size();

    std::vector<std::array<double,2>> boundaryNodePos = inputMesh.posOfNodes(boundaryNodes);
    // Boundary nodes are given in clockwise order, not column-major order
    std::vector<double> boundaryNodeValues; 
    boundaryNodeValues.resize(numBoundaryNodes);

    std::vector<double> boundaryCalc;
    std::array<double,2> *currPointer = boundaryNodePos.data();
    auto *currNodeValPointer = boundaryNodeValues.data();
    int ptrIncr;
    std::string prompt;

    for (int i=0; i<4; i++) {
        // Alternate between x and y-iteration
        ptrIncr = i % 2 == 0 ? xWidth : yWidth;
        prompt = strs[i];
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

void GetExtensionMatrices(QTM::QuadTreeMesh &inputMesh,
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

DvD ComputeSolutionStationary(SpD &StiffnessMatrix, SpD &PenaltyMatrix, SpD &FluxMatrix, SpD &fVec, SpD &columnSpace, SpD &nullSpace, DvD &boundaryVals) {
    // Eliminate rows and columns corr. to boundary nodes
    SpD CombinedMatrix = StiffnessMatrix - FluxMatrix - FluxMatrix.transpose() + PenaltyMatrix;
    SpD columnSpaceT = static_cast<SpD>(columnSpace.transpose());
    SpD A11 = columnSpaceT * CombinedMatrix * columnSpace;
    // Eliminate boundary rows and free columns
    SpD A12 = columnSpaceT * CombinedMatrix * nullSpace;
    // Eliminate boundary rows
    SpD F11 = columnSpaceT * fVec;

    Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int> > LuSolver;    
    LuSolver.analyzePattern(A11);
    LuSolver.factorize(A11);

    DvD x = LuSolver.solve(F11 - A12 * boundaryVals);
    x = columnSpace * x + nullSpace * boundaryVals;
    return x;
}

DD PoissonSolve(QTM::QuadTreeMesh &inputMesh,
                double c,
                double k,
                std::string source,
                std::vector<std::string> bcs) {
    
    std::vector<int> boundaryNodes = inputMesh.getBoundaryNodes();
    std::vector<int> freeNodes = inputMesh.getFreeNodes();
    int nNodes = inputMesh.nNodes();

    std::cout<<"Assembling stiffness matrix"<<std::endl;
    SpD KMatrix = StiffnessMatrix(inputMesh, k);
    std::cout<<"Assembling penalty matrix"<<std::endl;
    SpD PMatrix = PenaltyMatrix(inputMesh, k, 1);
    std::cout<<"Assembling flux matrix"<<std::endl;
    SpD SMatrix = FluxMatrix(inputMesh, k);

    std::cout<<"Assembling RHS vector"<<std::endl;
    SpD FMatrix = AssembleFVec(inputMesh, 1.0, source);
    std::cout<<"Assembling boundary condition vector"<<std::endl;
    DvD boundaryVals = EvalDirichletBoundaryCond(inputMesh, boundaryNodes, bcs);

    SpD nullSpace(nNodes, boundaryNodes.size());
    SpD columnSpace(nNodes, freeNodes.size());
    std::cout<<"Assembling null-orth matrix"<<std::endl;
    GetExtensionMatrices(inputMesh, boundaryNodes, freeNodes, nullSpace, columnSpace);
    std::cout<<"Solving system with "<<freeNodes.size()<<" nodes"<<std::endl;
    DD x = ComputeSolutionStationary(KMatrix, PMatrix, SMatrix, FMatrix, columnSpace, nullSpace, boundaryVals);
    std::cout<<"System solved!"<<std::endl;
    return x;
}

