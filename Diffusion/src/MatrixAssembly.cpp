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


SpD StiffnessMatrix(QTM::QuadTreeMesh& mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();
    
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
        // std::cout<<elm->CID<<"\n----------------\n";
        // std::cout<<combinedX.rows()<<", "<<combinedX.cols()<<std::endl;
        // std::cout<<combinedY.rows()<<", "<<combinedY.cols()<<std::endl;
        // std::cout<<combineXT.rows()<<", "<<combineXT.cols()<<std::endl;
        // std::cout<<combineYT.rows()<<", "<<combineYT.cols()<<std::endl;
        // std::cout<<coeffMat.rows()<<", "<<coeffMat.cols()<<std::endl;
        // std::cout<<weightMat.rows()<<", "<<weightMat.cols()<<std::endl;
        
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

SpD PenaltyMatrix(QTM::QuadTreeMesh& mesh, double k, double alpha) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = mesh.numElemNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();

    // basis func to node mapping
    DD B; B.setIdentity(numElemNodes, numElemNodes);

    // get basis func vals for split cell quad
    std::vector<double> BhInitializer; 
    BhInitializer.reserve((mesh.halfGaussPoints.size()) * numNodes);

    // get unit vecs
    DvD topVec = DvD::Zero(numNodes); topVec(0) = 1;
    DvD bottomVec = DvD::Zero(numNodes); bottomVec(deg) = 1;

    // Generate values for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<double> xVals = Utils::evalLagrangeInterp(k, mesh.halfGaussPoints, mesh.gaussPoints);
        BhInitializer.insert(BhInitializer.end(), xVals.begin(), xVals.end());
    }

    // map basis values to matrix
    Eigen::Map<DD> BhT(BhInitializer.data(), mesh.halfGaussPoints.size(), numNodes); // eval points are traversed first
    DD Bh = (DD)(BhT.transpose());

    std::array<DD,4> splitCellVals;

    DD splitCellPlaceholder(numElemNodes, mesh.halfGaussPoints.size());
    splitCellPlaceholder << Eigen::kroneckerProduct(bottomVec, Bh); splitCellVals[0] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(Bh, bottomVec); splitCellVals[1] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(topVec, Bh); splitCellVals[2] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(Bh, topVec); splitCellVals[3] = splitCellPlaceholder;

    // index vectors for split cell gauss points
    std::vector<int> frontIdxs(numNodes, 0);  
    std::vector<int> backIdxs(numNodes, 0);

    for (int i=0; i<numNodes; i++) {
        frontIdxs[i] = i;
        backIdxs[i] = i+deg;
    }

    std::vector<int> splitIdx[2] = { frontIdxs, backIdxs };
    
    // get quad weights in 1D
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
    std::vector<double> integX = Utils::integrateLagrange(gaussPoints);

    // Convert 1D quad weights to diag matrix
    Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    DD quadWeights1D = (DD)integXMat.asDiagonal();

    auto leaves = mesh.leaves;
    std::vector<QTM::Direction> directions = {QTM::Direction::N, QTM::Direction::E, 
                                                QTM::Direction::S, QTM::Direction::W};

                                                
    std::vector<QTM::Direction> oppdirs = {QTM::Direction::S, QTM::Direction::W, 
                                            QTM::Direction::N, QTM::Direction::E};

    std::vector<std::vector<int>> localNodes = {mesh.GetLocalBoundaryNodes(QTM::Direction::N),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::E),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::S),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::W)};

    double a; // penalty parameters

    std::vector<int> boundaryNodes;
    std::vector<int> neighborNodes;
    std::vector<int> neighborLocals;
    std::vector<int> elemNodes;
    std::vector<int> elemLocals;
    std::vector<std::shared_ptr<QTM::Cell>> neighbors;
    std::vector<Eigen::Triplet<double>> tripletList; tripletList.reserve(nNodes);
    for (auto &elm : leaves) {
        elemNodes = mesh.GetGlobalElemNodes(elm->CID);
        elemLocals = mesh.GetTrimmedLocalNodes(elm->CID, elemNodes);
        // std::cout<<"-----------"<<std::endl;
        // std::cout<<"elem: "<<elm->CID<<std::endl;
        // for (auto i : elemNodes) {
        //     std::cout<<i<<std::endl;
        // }
        // std::cout<<"neighbors:\n";
        // get neighbors
        for (auto dir : directions) {
            // std::cout<<"neighbor in direction "<<dir<<std::endl;
            QTM::Direction oppdir = oppdirs[dir];
            neighbors = mesh.GetCellNeighbors(dir, elm->CID);
            if (neighbors[0] == nullptr) {
                continue;
            }
            for (int NI = 0; NI < neighbors.size(); NI++) { 
                auto neighbor = neighbors[NI];
                auto jac = std::min(elm->width, neighbor->width);
                // calculate penalty param
                a = alpha/jac/2;
                if (neighbor->CID < elm->CID || elm->level < neighbor->level) { // case appropriate neighbor exists
                    neighborNodes = mesh.GetGlobalElemNodes(neighbor->CID);
                    neighborLocals = mesh.GetTrimmedLocalNodes(neighbor->CID, neighborNodes);
                } else { 
                    continue;
                }
                // calculate jump matrix
                DD topJump;
                DD bottomJump = -B(neighborLocals, localNodes[oppdir]);
                if (elm->level == neighbor->level) {
                    topJump = B(elemLocals, localNodes[dir]);
                } else {
                    topJump = splitCellVals[dir](elemLocals, splitIdx[NI]);
                }
            
                DD jumpMatrix(elemNodes.size() + neighborNodes.size(), numNodes);
                // std::cout<<topJump.rows()<<", "<<topJump.cols()<<std::endl;
                // std::cout<<bottomJump.rows()<<", "<<bottomJump.cols()<<std::endl;
                // std::cout<<elemNodes.size()<<", "<<neighborNodes.size()<<", "<<neighborLocals.size()<<std::endl;
                // std::cout<<"-----------"<<std::endl;
                jumpMatrix << topJump, bottomJump;

                DD jumpMatrixT = (DD)(jumpMatrix.transpose());
                DD localElemMat = (DD)(a * jumpMatrix * jac*quadWeights1D * jumpMatrixT);

                boundaryNodes.reserve(elemNodes.size() + neighborNodes.size());
                boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                for (int j=0; j<boundaryNodes.size(); j++) {
                    for (int i=0; i<boundaryNodes.size(); i++) {
                        if (localElemMat(i,j) != 0) {
                            tripletList.emplace_back(boundaryNodes[i],boundaryNodes[j],localElemMat(i,j));
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

SpD FluxMatrix(QTM::QuadTreeMesh& mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = mesh.numElemNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);

    // get quad weights in 1D
    std::vector<double> integX = Utils::integrateLagrange(gaussPoints);

    // Convert 1D quad weights to diag matrix
    Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    DD quadWeights1D = (DD)integXMat.asDiagonal();

    // basis func to node mapping
    DD B; B.setIdentity(numElemNodes, numElemNodes);
    DD Bs; Bs.setIdentity(numNodes, numNodes);

    // Generate derivative matrix
    std::vector<double> AInitializer; 
    AInitializer.reserve(numNodes * numNodes);

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

    // get basis func vals for split cell quad
    std::vector<double> BhInitializer; 
    BhInitializer.reserve((mesh.halfGaussPoints.size()) * numNodes);

    // get unit vecs
    DvD topVec = DvD::Zero(numNodes); topVec(0) = 1;
    DvD bottomVec = DvD::Zero(numNodes); bottomVec(deg) = 1;

    DD topRowVec = DD::Zero(1,numNodes); topRowVec(0,0) = 1;
    DD bottomRowVec = DD::Zero(1,numNodes); bottomRowVec(0,deg) = 1;

    DD fullRowVec = DD::Ones(1,numNodes);

    // Generate values for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<double> xVals = Utils::evalLagrangeInterp(k, mesh.halfGaussPoints, mesh.gaussPoints);
        BhInitializer.insert(BhInitializer.end(), xVals.begin(), xVals.end());
    }

    // map basis values to matrix
    Eigen::Map<DD> BhT(BhInitializer.data(), mesh.halfGaussPoints.size(), numNodes); // eval points are traversed first
    DD Bh = (DD)(BhT.transpose());

    std::array<DD,4> splitCellVals;

    DD splitCellPlaceholder(numElemNodes, mesh.halfGaussPoints.size());
    splitCellPlaceholder << Eigen::kroneckerProduct(bottomVec, Bh); splitCellVals[0] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(Bh, bottomVec); splitCellVals[1] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(topVec, Bh); splitCellVals[2] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(Bh, topVec); splitCellVals[3] = splitCellPlaceholder;

    // get basis func gradients for split cell quad
    std::vector<double> AhInitializer; 
    AhInitializer.reserve((mesh.halfGaussPoints.size()) * numNodes);

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<mesh.halfGaussPoints.size(); k++) { 
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, mesh.halfGaussPoints, gaussPoints);
        AhInitializer.insert(AhInitializer.end(), xPartials.begin(), xPartials.end());
    }
    
    // map derivative values to matrix
    Eigen::Map<DD> Ah(AhInitializer.data(), mesh.halfGaussPoints.size(), numNodes); 

    DD topGrads = A(0, Eigen::all);
    DD bottomGrads = A(numNodes-1, Eigen::all);

    std::array<DD,4> splitCellGradsX;
    
    DD splitCellGradXPlaceholder(mesh.halfGaussPoints.size(), numElemNodes);
    splitCellGradXPlaceholder << Eigen::kroneckerProduct(bottomRowVec, Ah); splitCellGradsX[0] = splitCellGradXPlaceholder;
    splitCellGradXPlaceholder << Eigen::kroneckerProduct(BhT, bottomGrads); splitCellGradsX[1] = splitCellGradXPlaceholder;
    splitCellGradXPlaceholder << Eigen::kroneckerProduct(topRowVec, Ah); splitCellGradsX[2] = splitCellGradXPlaceholder;
    splitCellGradXPlaceholder << Eigen::kroneckerProduct(BhT, topGrads); splitCellGradsX[3] = splitCellGradXPlaceholder;

    std::array<DD,4> splitCellGradsY;
    DD splitCellGradYPlaceholder(mesh.halfGaussPoints.size(), numElemNodes); 
    splitCellGradYPlaceholder << Eigen::kroneckerProduct(bottomGrads, BhT); splitCellGradsY[0] = splitCellGradYPlaceholder;
    splitCellGradYPlaceholder << Eigen::kroneckerProduct(Ah, bottomRowVec); splitCellGradsY[1] = splitCellGradYPlaceholder;
    splitCellGradYPlaceholder << Eigen::kroneckerProduct(topGrads, BhT); splitCellGradsY[2] = splitCellGradYPlaceholder;
    splitCellGradYPlaceholder << Eigen::kroneckerProduct(Ah, topRowVec); splitCellGradsY[3] = splitCellGradYPlaceholder;

    // index vectors for split cell gauss points
    std::vector<int> frontIdxs(numNodes, 0);  
    std::vector<int> backIdxs(numNodes, 0);

    for (int i=0; i<numNodes; i++) {
        frontIdxs[i] = i;
        backIdxs[i] = i+deg;
    }

    std::vector<int> splitIdx[2] = { frontIdxs, backIdxs };

    double normalX[4] = {0,1,0,-1};
    double normalY[4] = {1,0,-1,0};

    std::vector<int> boundaryNodes;
    std::vector<int> neighborNodes;
    std::vector<int> neighborLocals;
    std::vector<int> elemNodes;
    std::vector<int> elemLocals;
    std::vector<std::shared_ptr<QTM::Cell>> neighbors;
    std::vector<Eigen::Triplet<double>> tripletList; tripletList.reserve(nNodes);
    auto leaves = mesh.GetAllCells();
    std::vector<QTM::Direction> directions = {QTM::Direction::N, QTM::Direction::E, 
                                                QTM::Direction::S, QTM::Direction::W};

    std::vector<QTM::Direction> oppdirs = {QTM::Direction::S, QTM::Direction::W, 
                                            QTM::Direction::N, QTM::Direction::E};

    std::vector<std::vector<int>> localNodes = {mesh.GetLocalBoundaryNodes(QTM::Direction::N),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::E),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::S),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::W)};
    for (auto &elm : leaves) {
        elemNodes = mesh.GetGlobalElemNodes(elm->CID);
        elemLocals = mesh.GetTrimmedLocalNodes(elm->CID, elemNodes);
        // get neighbors
        for (auto dir : directions) {
            neighbors = mesh.GetCellNeighbors(dir, elm->CID);
            if (neighbors[0] == nullptr) {
                continue;
            }
            QTM::Direction oppdir = oppdirs[dir];
            for (int NI = 0; NI < neighbors.size(); NI++) { 
                auto neighbor = neighbors[NI];
                auto jac = std::min(elm->width, neighbor->width);
                if (neighbor->CID < elm->CID || elm->level < neighbor->level) { // case appropriate neighbor exists
                    neighborNodes = mesh.GetGlobalElemNodes(neighbor->CID);
                    neighborLocals = mesh.GetTrimmedLocalNodes(neighbor->CID, neighborNodes);
                } else { 
                    continue;
                }
                // jump matrix setup
                DD topJump;
                DD bottomJump = -B(neighborLocals, localNodes[oppdir]);

                // flux matrix setup
                DD topGradX;
                DD topGradY;
                DD bottomGradX = normalX[oppdir] * combinedX(localNodes[oppdir], neighborLocals)/jac;
                DD bottomGradY = normalY[oppdir] * combinedY(localNodes[oppdir], neighborLocals)/jac;

                if (elm->level == neighbor->level) {
                    topJump = B(elemLocals, localNodes[dir]);
                    topGradX = normalX[dir] * combinedX(localNodes[dir], elemLocals)/jac;
                    topGradY = normalY[dir] * combinedY(localNodes[dir], elemLocals)/jac;
                } else {
                    topJump = splitCellVals[dir](elemLocals, splitIdx[NI]);
                    topGradX = normalX[dir] * splitCellGradsX[dir](splitIdx[NI], elemLocals)/jac/2;
                    topGradY = normalY[dir] * splitCellGradsY[dir](splitIdx[NI], elemLocals)/jac/2;
                }

                // calculate jump matrix
                DD jumpMatrix(elemNodes.size() + neighborNodes.size(), numNodes);
                jumpMatrix << topJump, bottomJump;

                // calculate flux matrix
                DD fluxMatrixX(numNodes, elemNodes.size() + neighborNodes.size());
                DD fluxMatrixY(numNodes, elemNodes.size() + neighborNodes.size());

                // place partial derivatives in combined mat
                fluxMatrixX << topGradX, bottomGradX;
                fluxMatrixY << topGradY, bottomGradY;

                boundaryNodes.reserve(elemNodes.size() + neighborNodes.size()); // aggregated nodes in current and neighbor cell
                boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                // assemble local matrix 
                DD localElemMat = (DD)(.5 * jac * jumpMatrix * quadWeights1D * fluxMatrixX + 
                                    .5 * jac * jumpMatrix * quadWeights1D * fluxMatrixY);

                for (int j=0; j<boundaryNodes.size(); j++) {
                    for (int i=0; i<boundaryNodes.size(); i++) {
                        if (localElemMat(i,j) != 0) {
                            tripletList.emplace_back(boundaryNodes[i],boundaryNodes[j],localElemMat(i,j));
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

SpD AssembleFVec(QTM::QuadTreeMesh& mesh, double f, std::string evalStr) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();
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
        double jac = elm->width; // Jacobian factors
        // calculate local matrix
        std::vector<double> collectSourceVals; collectSourceVals.reserve(numElemNodes);
        auto nodes = elm->nodes;
        for (int i=nodes[0]; i<=nodes[1]; i++) {
            collectSourceVals.push_back(fEval[i]);
        }
        Eigen::Map<DvD> sourceVector(collectSourceVals.data(), numElemNodes, 1);

        localElemMat = weightMat*sourceVector*jac*jac;
        // Get nodes in element
        
        // Generate i,j,v triplets
        for (int i=0; i<numElemNodes; i++) {
            tripletList.emplace_back(nodes[0]+i, 0, localElemMat(i));
        }  
    }

    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,1);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

std::vector<double> ComputeResiduals(QTM::QuadTreeMesh& mesh, DvD& solution, SpD& source) {
    std::vector<double> out;

    // get quadrature weights

    // second derivatives


    for (auto leaf : mesh.leaves) {

    }
    return out;
}

// std::vector<std::shared_ptr<QTM::Cell>> TestResiduals(DvD& solution, QTM::QuadTreeMesh& mesh, double residualLimit) {
//     std::vector<std::shared_ptr<QTM::Cell>> out; out.reserve(mesh.leaves.size());
//     for (auto leaf : mesh.leaves) {
//         double res = ComputeResidual();
//         if (leaf->level < 5 && ComputeResidual() > residualLimit) {
//              out.push_back(leaf);
//         }
//     }

//     return out;
// }

DvD EvalDirichletBoundaryCond(QTM::QuadTreeMesh& inputMesh, std::vector<std::vector<int>>& boundaryNodes, std::vector<int>& allBoundaryNodes, std::vector<std::string>& strs) {
    std::vector<std::array<double,2>> boundaryNodePos;
    for (auto bNodes : boundaryNodes) {
        auto nodePos = inputMesh.GetNodePos(bNodes);
        boundaryNodePos.insert(boundaryNodePos.end(), nodePos.begin(), nodePos.end());
    }
    int numBoundaryNodes = allBoundaryNodes.size();

    // Boundary nodes are given in clockwise order, not column-major order
    std::vector<double> boundaryNodeValues; 
    boundaryNodeValues.reserve(numBoundaryNodes);

    std::vector<double> boundaryCalc;
    std::array<double,2> *currPointer = boundaryNodePos.data();
    int ptrIncr;
    std::string prompt;

    for (int i=0; i<4; i++) {
        ptrIncr = boundaryNodes[i].size();
        // Take bcFunc and evaluate it
        boundaryCalc = Utils::EvalSymbolicBC(currPointer, ptrIncr, strs[i]);
        boundaryNodeValues.insert(boundaryNodeValues.end(), boundaryCalc.begin(), boundaryCalc.end());
        currPointer += ptrIncr;
    }

    auto RmOrder = allBoundaryNodes;
    std::sort(RmOrder.begin(), RmOrder.end());

    auto BcValsSorted = Utils::ReshuffleNodeVals(RmOrder, allBoundaryNodes, boundaryNodeValues);

    Eigen::Map<DvD> boundaryNodeValuesVec(BcValsSorted.data(), numBoundaryNodes, 1);
    return (DvD)boundaryNodeValuesVec;
}

void GetExtensionMatrices(QTM::QuadTreeMesh& inputMesh,
                                        std::vector<int>& boundaryNodes, 
                                        std::vector<int>& freeNodes,
                                        SpD& nullSpace,
                                        SpD& columnSpace) {
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

DvD ComputeSolutionStationary(SpD& StiffnessMatrix, SpD& PenaltyMatrix, SpD& FluxMatrix, SpD& fVec, SpD& columnSpace, SpD& nullSpace, DvD& boundaryVals) {
    SpD FluxMatrixT = (SpD)(FluxMatrix.transpose());
    SpD CombinedMatrix = StiffnessMatrix - FluxMatrix - FluxMatrixT + PenaltyMatrix;
    // Eliminate rows and columns corr. to boundary nodes
    SpD columnSpaceT = (SpD)(columnSpace.transpose());
    SpD nullSpaceT = (SpD)(nullSpace.transpose());

    // Eliminate boundary rows and boundary columns
    SpD A11 = columnSpaceT * CombinedMatrix * columnSpace;
    // Eliminate boundary rows and free columns
    SpD A12 = columnSpaceT * CombinedMatrix * nullSpace;
    // Eliminate boundary rows
    SpD F11 = columnSpaceT * fVec;

    // // get Kbar
    // Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolver1;    
    // LuSolver1.analyzePattern(A22);
    // LuSolver1.factorize(A22);
    // SpD KBarInterm = LuSolver1.solve(A21);
    // SpD Kbar = A11 - A12*KBarInterm;

    // // get Fbar
    // Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolver2;    
    // LuSolver2.analyzePattern(A22);
    // LuSolver2.factorize(A22);
    // SpD FBarInterm = LuSolver2.solve(F22);
    // SpD Fbar = F11 - A12 * FBarInterm;

    // Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolver;    
    // LuSolver.analyzePattern(Kbar);
    // LuSolver.factorize(Kbar);
    // DvD x = LuSolver.solve(Fbar);

    using namespace Eigen;

    Eigen::ConjugateGradient<SpD, Lower|Upper> cg;
    cg.compute(A11);
    DvD x = cg.solve(F11 - A12 * boundaryVals);


    // Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolver;    
    // LuSolver.analyzePattern(A11);
    // LuSolver.factorize(A11);
    // DvD x = LuSolver.solve(F11 - A12 * boundaryVals);

    x = columnSpace * x + nullSpace * boundaryVals;
    return x;
}

DD PoissonSolve(QTM::QuadTreeMesh& inputMesh,
                double c,
                double k,
                std::string source,
                std::vector<std::string> bcs,
                double penaltyParam) {
    
    auto boundaryNodes = inputMesh.boundaryNodes;
    std::vector<int> freeNodes = inputMesh.freeNodes;
    int nNodes = inputMesh.nNodes();

    std::vector<int> allBoundaryNodes;
    for (auto bNodes : boundaryNodes) {
        allBoundaryNodes.insert(allBoundaryNodes.end(), bNodes.begin(), bNodes.end());
    }

    std::cout<<"Assembling stiffness matrix"<<std::endl;
    SpD KMatrix = StiffnessMatrix(inputMesh, k);
    std::cout<<"Assembling penalty matrix"<<std::endl;
    SpD PMatrix = PenaltyMatrix(inputMesh, k, penaltyParam);
    std::cout<<"Assembling flux matrix"<<std::endl;
    SpD SMatrix = FluxMatrix(inputMesh, k);

    std::cout<<"Assembling RHS vector"<<std::endl;
    SpD FMatrix = AssembleFVec(inputMesh, 1.0, source);
    std::cout<<"Assembling boundary condition vector"<<std::endl;
    DvD boundaryVals = EvalDirichletBoundaryCond(inputMesh, boundaryNodes, allBoundaryNodes, bcs);

    SpD nullSpace(nNodes, allBoundaryNodes.size());
    SpD columnSpace(nNodes, freeNodes.size());
    std::cout<<"Assembling null-orth matrix"<<std::endl;
    GetExtensionMatrices(inputMesh, allBoundaryNodes, freeNodes, nullSpace, columnSpace);
    std::cout<<"Solving system with "<<freeNodes.size()<<" nodes"<<std::endl;
    DD x = ComputeSolutionStationary(KMatrix, PMatrix, SMatrix, FMatrix, columnSpace, nullSpace, boundaryVals);
    std::cout<<"System solved!"<<std::endl;
    return x;
}