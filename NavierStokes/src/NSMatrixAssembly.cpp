#include "..\include\NSMatrixAssembly.hpp"

SpD PVMatrix(QTM::QuadTreeMesh& mesh, int diffDir) {
    // Oriented as V-P.
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
    Eigen::Map<DD> AT(AInitializer.data(), numNodes, numNodes); 
    DD A = (DD)(AT.transpose());
    // Generate quadrature weight matrices
    DD weightMat = GenerateQuadWeights(gaussPoints, gaussPoints, numNodes, numNodes, numElemNodes);
    // Generate mass matrices
    DD B; B.setIdentity(numNodes, numNodes);
    
    // Get element-wise matrix intermediates
    DD combined(numElemNodes, numElemNodes);
    switch (diffDir) {
        case 1: {
            combined << Eigen::kroneckerProduct(B, A);
            break;
        }
        case 2: {
            combined << Eigen::kroneckerProduct(A, B); 
            break; 
        }
    }

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    auto combineT = (DD)combined.transpose();

    // Integrate over all elements
    DD localElemMat(numElemNodes, numElemNodes);
    auto leaves = mesh.GetAllCells();
    for (auto &elm : leaves) {
        auto jac = elm->width;
        // calculate local matrix
        localElemMat = jac*combined*weightMat;
 
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
    mat.makeCompressed();
    return mat;
}

SpD PressureJumpMatrix(QTM::QuadTreeMesh& mesh, int diffDir) {
    // pressure jump, velocity flux average, U-Q
    // <u dot n> [q] 
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

    std::array<DD,4> splitCellValsT;

    DD splitCellPlaceholderT(mesh.halfGaussPoints.size(), numElemNodes);
    splitCellPlaceholderT << Eigen::kroneckerProduct(bottomRowVec, BhT); splitCellValsT[0] = splitCellPlaceholderT;
    splitCellPlaceholderT << Eigen::kroneckerProduct(BhT, bottomRowVec); splitCellValsT[1] = splitCellPlaceholderT;
    splitCellPlaceholderT << Eigen::kroneckerProduct(topRowVec, BhT); splitCellValsT[2] = splitCellPlaceholderT;
    splitCellPlaceholderT << Eigen::kroneckerProduct(BhT, topRowVec); splitCellValsT[3] = splitCellPlaceholderT;

    // index vectors for split cell gauss points
    std::vector<int> frontIdxs(numNodes, 0);  
    std::vector<int> backIdxs(numNodes, 0);

    for (int i=0; i<numNodes; i++) {
        frontIdxs[i] = i;
        backIdxs[i] = i+deg;
    }

    std::vector<int> splitIdx[2] = { frontIdxs, backIdxs };

    double normals[4] = {0,0,0,0};
    double normalsX[4] = {0,1,0,-1};
    double normalsY[4] = {1,0,-1,0};
    switch (diffDir) {
        case 1: {
            std::copy(normalsX, normalsX+4, normals);
            break;
        }
        case 2: {
            std::copy(normalsY, normalsY+4, normals);
            break;
        }
    }

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
    double fac;
    for (auto &elm : leaves) {
        elemNodes = mesh.GetGlobalElemNodes(elm->CID);
        elemLocals = mesh.GetTrimmedLocalNodes(elm->CID, elemNodes);
        // get neighbors
        for (auto dir : directions) {
            neighbors = mesh.GetCellNeighbors(dir, elm->CID);
            QTM::Direction oppdir = oppdirs[dir];
            for (int NI = 0; NI < neighbors.size(); NI++) { 
                auto neighbor = neighbors[NI];
                if (!neighbor) {
                    // continue;
                    neighborNodes = {};
                    neighborLocals = {};
                    neighbor = elm;
                    fac = 1;
                }
                else if (neighbor->CID < elm->CID || elm->level < neighbor->level) { // case appropriate neighbor exists
                    neighborNodes = mesh.GetGlobalElemNodes(neighbor->CID);
                    neighborLocals = mesh.GetTrimmedLocalNodes(neighbor->CID, neighborNodes);
                    fac = .5;
                } else { 
                    continue;
                }
                auto jac = std::min(elm->width, neighbor->width);
                // jump matrix setup
                DD topJump;
                DD bottomJump = -B(neighborLocals, localNodes[oppdir]);

                DD topFlux;
                DD bottomFlux = normals[oppdir] * B(localNodes[oppdir], neighborLocals);

                if (elm->level == neighbor->level) {
                    topJump = B(elemLocals, localNodes[dir]);
                    topFlux = normals[dir] * B(localNodes[dir], elemLocals);
                } else {
                    topJump = splitCellVals[dir](elemLocals, splitIdx[NI]);
                    topFlux = normals[dir] * splitCellValsT[dir](splitIdx[NI], elemLocals);
                }

                // calculate jump matrix
                DD jumpMatrix(elemNodes.size() + neighborNodes.size(), numNodes);
                jumpMatrix << topJump, bottomJump;

                // calculate flux matrix
                DD fluxMatrix(numNodes, elemNodes.size() + neighborNodes.size());
                fluxMatrix << topFlux, bottomFlux;

                boundaryNodes.reserve(elemNodes.size() + neighborNodes.size()); // aggregated nodes in current and neighbor cell
                boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                // assemble local matrix 
                DD localElemMat = (DD)(fac * jac * jumpMatrix *  quadWeights1D * fluxMatrix);

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
    mat.makeCompressed();
    return mat;
}

SpD PressureAvgMatrix(QTM::QuadTreeMesh& mesh, int diffDir) {
    // pressure average, velocity flux jump, V-P
    // [v dot n] <p>
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

    std::array<DD,4> splitCellValsT;

    DD splitCellPlaceholderT(mesh.halfGaussPoints.size(), numElemNodes);
    splitCellPlaceholderT << Eigen::kroneckerProduct(bottomRowVec, BhT); splitCellValsT[0] = splitCellPlaceholderT;
    splitCellPlaceholderT << Eigen::kroneckerProduct(BhT, bottomRowVec); splitCellValsT[1] = splitCellPlaceholderT;
    splitCellPlaceholderT << Eigen::kroneckerProduct(topRowVec, BhT); splitCellValsT[2] = splitCellPlaceholderT;
    splitCellPlaceholderT << Eigen::kroneckerProduct(BhT, topRowVec); splitCellValsT[3] = splitCellPlaceholderT;

    // index vectors for split cell gauss points
    std::vector<int> frontIdxs(numNodes, 0);  
    std::vector<int> backIdxs(numNodes, 0);

    for (int i=0; i<numNodes; i++) {
        frontIdxs[i] = i;
        backIdxs[i] = i+deg;
    }

    std::vector<int> splitIdx[2] = { frontIdxs, backIdxs };

    double normals[4] = {0,0,0,0};
    double normalsX[4] = {0,1,0,-1};
    double normalsY[4] = {1,0,-1,0};
    switch (diffDir) {
        case 1: {
            std::copy(normalsX, normalsX+4, normals);
            break;
        }
        case 2: {
            std::copy(normalsY, normalsY+4, normals);
            break;
        }
    }

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
    double fac;
    for (auto &elm : leaves) {
        elemNodes = mesh.GetGlobalElemNodes(elm->CID);
        elemLocals = mesh.GetTrimmedLocalNodes(elm->CID, elemNodes);
        // get neighbors
        for (auto dir : directions) {
            neighbors = mesh.GetCellNeighbors(dir, elm->CID);
            QTM::Direction oppdir = oppdirs[dir];
            for (int NI = 0; NI < neighbors.size(); NI++) { 
                auto neighbor = neighbors[NI];
                if (!neighbor) {
                    neighbor = elm;
                    neighborNodes = {};
                    neighborLocals = {};
                    fac = 1;
                    // continue;
                }
                else if (neighbor->CID < elm->CID || elm->level < neighbor->level) { // case appropriate neighbor exists
                    neighborNodes = mesh.GetGlobalElemNodes(neighbor->CID);
                    neighborLocals = mesh.GetTrimmedLocalNodes(neighbor->CID, neighborNodes);
                    fac = .5;
                } else { 
                    continue;
                }
                auto jac = std::min(elm->width, neighbor->width);
                // jump matrix setup
                DD topJump;
                DD bottomJump = normals[oppdir] * B(neighborLocals, localNodes[oppdir]);

                DD topAvg;
                DD bottomAvg =  B(localNodes[oppdir], neighborLocals);

                if (elm->level == neighbor->level) {
                    topJump = normals[dir] * B(elemLocals, localNodes[dir]);
                    topAvg = B(localNodes[dir], elemLocals);
                } else {
                    topJump = normals[dir] * splitCellVals[dir](elemLocals, splitIdx[NI]);
                    topAvg = splitCellValsT[dir](splitIdx[NI], elemLocals);
                }

                // calculate jump matrix
                DD jumpMatrix(elemNodes.size() + neighborNodes.size(), numNodes);
                jumpMatrix << topJump, bottomJump;

                // calculate flux matrix
                DD fluxMatrix(numNodes, elemNodes.size() + neighborNodes.size());
                fluxMatrix << topAvg, bottomAvg;

                boundaryNodes.reserve(elemNodes.size() + neighborNodes.size()); // aggregated nodes in current and neighbor cell
                boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                // assemble local matrix 
                DD localElemMat = (DD)(fac * jac * jumpMatrix * quadWeights1D * fluxMatrix);

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
    mat.makeCompressed();
    return mat;
}

SpD PressurePenaltyMatrix(QTM::QuadTreeMesh& mesh, double mu) {
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
            if (!neighbors[0]) {
                continue;
            }
            for (int NI = 0; NI < neighbors.size(); NI++) { 
                auto neighbor = neighbors[NI];
                auto jac = std::min(elm->width, neighbor->width);
                // calculate penalty param
                a = 2*jac;
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
                        tripletList.emplace_back(boundaryNodes[i],boundaryNodes[j],localElemMat(i,j));
                    }
                }    
                boundaryNodes.clear();
            }
        }
    }
    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    mat.makeCompressed();
    return mat;
}

SpD PressureFluxMatrix(QTM::QuadTreeMesh& mesh, int diffDir) {
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

    std::array<DD,4> splitCellValsT;

    DD splitCellPlaceholderT(mesh.halfGaussPoints.size(), numElemNodes);
    splitCellPlaceholderT << Eigen::kroneckerProduct(bottomRowVec, BhT); splitCellValsT[0] = splitCellPlaceholderT;
    splitCellPlaceholderT << Eigen::kroneckerProduct(BhT, bottomRowVec); splitCellValsT[1] = splitCellPlaceholderT;
    splitCellPlaceholderT << Eigen::kroneckerProduct(topRowVec, BhT); splitCellValsT[2] = splitCellPlaceholderT;
    splitCellPlaceholderT << Eigen::kroneckerProduct(BhT, topRowVec); splitCellValsT[3] = splitCellPlaceholderT;

    // index vectors for split cell gauss points
    std::vector<int> frontIdxs(numNodes, 0);  
    std::vector<int> backIdxs(numNodes, 0);

    for (int i=0; i<numNodes; i++) {
        frontIdxs[i] = i;
        backIdxs[i] = i+deg;
    }

    std::vector<int> splitIdx[2] = { frontIdxs, backIdxs };

    double normals[4] = {0,0,0,0};
    double normalsX[4] = {0,1,0,-1};
    double normalsY[4] = {1,0,-1,0};
    switch (diffDir) {
        case 1: {
            std::copy(normalsX, normalsX+4, normals);
            break;
        }
        case 2: {
            std::copy(normalsY, normalsY+4, normals);
            break;
        }
    }

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
            if (!neighbors[0]) {
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
                DD bottomJump = -B(localNodes[oppdir], neighborLocals);

                DD topFlux;
                DD bottomFlux = normals[oppdir] * B(neighborLocals, localNodes[oppdir]);

                if (elm->level == neighbor->level) {
                    topJump = B(localNodes[dir], elemLocals);
                    topFlux = normals[dir] * B(elemLocals, localNodes[dir]);
                } else {
                    topJump = splitCellValsT[dir](splitIdx[NI], elemLocals);
                    topFlux = normals[dir] * splitCellVals[dir](elemLocals, splitIdx[NI]);
                }

                // calculate jump matrix
                DD jumpMatrix(numNodes, elemNodes.size() + neighborNodes.size());
                jumpMatrix << topJump, bottomJump;

                // calculate flux matrix
                DD fluxMatrix(elemNodes.size() + neighborNodes.size(), numNodes);
                fluxMatrix << topFlux, bottomFlux;

                boundaryNodes.reserve(elemNodes.size() + neighborNodes.size()); // aggregated nodes in current and neighbor cell
                boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                // assemble local matrix 
                DD localElemMat = (DD)(.5 * jac * fluxMatrix *  quadWeights1D * jumpMatrix);

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
    mat.makeCompressed();
    return mat;
}

// SpD ConvectionMatrix(QTM::QuadTreeMesh& mesh, double rho, DvD&& U, DvD&& V, DvD&& P) {

// }

// SpD ConvectionMatrix(QTM::QuadTreeMesh& mesh, double rho, DvD& state) {

// }

// SpD GeneralizedNeumannBC(QTM::QuadTreeMesh& mesh, double mu) {
//     // Implements boundary condition of form mu * n dot grad(u) - p*n = 0
// }

SpD AssembleStokesSystem(QTM::QuadTreeMesh& mesh, 
                            SpD& dgPoissonMat,
                            SpD& mat1, 
                            SpD& mat2,
                            SpD& mat1T,
                            SpD& mat2T) {
    int nNodes = mesh.nNodes();
    int offset1 = nNodes; 
    int offset2 = 2*nNodes;
    int offset3 = 3*nNodes;
    int numElemNodes = mesh.numElemNodes;
    std::vector<double> gaussPoints = mesh.gaussPoints;
    SpD PressurePenalty = PressurePenaltyMatrix(mesh, 1);
    int totalNNZ = 2*dgPoissonMat.nonZeros() + PressurePenalty.nonZeros() +
                    mat1T.nonZeros() + mat2T.nonZeros() + 
                    mat1.nonZeros() + mat2.nonZeros() + 2*nNodes;
    std::vector<Eigen::Triplet<double>> tripletList; tripletList.reserve(totalNNZ);

    // insert diffusive terms
    for (int k = 0; k < dgPoissonMat.outerSize(); ++k) {
        for (SpD::InnerIterator it(dgPoissonMat, k); it; ++it) {
            tripletList.emplace_back(it.row(), it.col(), it.value());
            tripletList.emplace_back(it.row()+offset1, it.col()+offset1, it.value());
        }
    }

    for (int k = 0; k < mat1.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat1, k); it; ++it) {
            tripletList.emplace_back(it.row(), it.col()+offset2, it.value());
        }
    }

    for (int k = 0; k < mat2.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat2, k); it; ++it) {
            tripletList.emplace_back(it.row()+offset1, it.col()+offset2, it.value());
        }
    }

    for (int k = 0; k < mat1T.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat1T, k); it; ++it) {
            tripletList.emplace_back(it.row()+offset2, it.col(), it.value());
        }
    }

    for (int k = 0; k < mat2T.outerSize(); ++k) {
        for (SpD::InnerIterator it(mat2T, k); it; ++it) {
            tripletList.emplace_back(it.row()+offset2, it.col()+offset1, it.value());
        }
    }

    // for (int k = 0; k < PressurePenalty.outerSize(); ++k) {
    //     for (SpD::InnerIterator it(PressurePenalty, k); it; ++it) {
    //         tripletList.emplace_back(it.row()+offset2, it.col()+offset2, it.value());
    //     }
    // }

    // add zero mean pressure constraint 
    int currDofIndx = offset2;
    auto weightMat = GenerateQuadWeights(gaussPoints, gaussPoints, mesh.deg+1, mesh.deg+1, numElemNodes);
    std::vector<double> weights(weightMat.diagonal().begin(), weightMat.diagonal().end());
    for (auto elm : mesh.leaves) {
        double area = (elm->width * elm->width);
        for (int i=0; i<numElemNodes; i++) {
            tripletList.emplace_back(offset3, currDofIndx, area * weights[i]);
            tripletList.emplace_back(currDofIndx, offset3, area * weights[i]);
            currDofIndx++;
        }
    }

    SpD mat(3*nNodes+1,3*nNodes+1);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

DvD AssembleStokesSource(QTM::QuadTreeMesh& mesh, 
                            DvD& XSource,
                            DvD& YSource) {
    int nNodes = mesh.nNodes();
    int offset1 = nNodes;
    int totalNNZ = XSource.nonZeros() + YSource.nonZeros();
    // std::vector<Eigen::Triplet<double>> tripletList; tripletList.reserve(totalNNZ);

    std::vector<double> outInit;
    outInit.insert(outInit.begin(), XSource.begin(), XSource.end());
    outInit.insert(outInit.begin(), YSource.begin(), YSource.end());

    DvD out(3*nNodes+1);
    DvD allZeros = DvD::Zero(nNodes+1);

    out << XSource, YSource, allZeros;

    // Eigen::Map<DvD> outP(outInit.data(),3*nNodes+1);
    // DvD out = (DvD)outP;

    // // insert X source
    // for (int k = 0; k < XSource.outerSize(); ++k) {
    //     for (SpD::InnerIterator it(XSource, k); it; ++it) {
    //         tripletList.emplace_back(it.row(), it.col(), it.value());
    //     }
    // }

    // // insert Y source
    // for (int k = 0; k < YSource.outerSize(); ++k) {
    //     for (SpD::InnerIterator it(YSource, k); it; ++it) {
    //         tripletList.emplace_back(it.row()+offset1, it.col(), it.value());
    //     }
    // }

    // SpD mat(3*nNodes+1,1);
    // mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return out;
}

DvD EvalPartialSymbolicBoundaryCond(QTM::QuadTreeMesh& inputMesh, 
                                    std::vector<std::vector<int>>& boundaryNodes, 
                                    std::vector<std::string>& strs, 
                                    std::vector<int>& allBoundaryNodes, 
                                    int offset) {
    std::vector<std::array<double,2>> boundaryNodePos;
    int nNodes = inputMesh.nNodes();
    int totalOffset = offset*nNodes;
    std::vector<int> allDirBoundary;
    for (auto bNodes : boundaryNodes) {
        allDirBoundary.insert(allDirBoundary.end(), bNodes.begin(), bNodes.end());
    }
    for (int bnd : allDirBoundary) {
        allBoundaryNodes.push_back(totalOffset+bnd);
    }
    for (auto bNodes : boundaryNodes) {
        auto nodePos = inputMesh.GetNodePos(bNodes);
        boundaryNodePos.insert(boundaryNodePos.end(), nodePos.begin(), nodePos.end());
    }

    int numBoundaryNodes = allDirBoundary.size();

    // Boundary nodes are given in clockwise order, not column-major order
    std::vector<double> boundaryNodeValues; 
    boundaryNodeValues.reserve(numBoundaryNodes);

    std::vector<double> boundaryCalc;
    std::array<double,2> *currPointer = boundaryNodePos.data();
    int ptrIncr;

    for (int i=0; i<boundaryNodes.size(); i++) {
        ptrIncr = boundaryNodes[i].size();
        // Take bcFunc and evaluate it
        boundaryCalc = Utils::EvalSymbolicBC(currPointer, ptrIncr, strs[i]);
        boundaryNodeValues.insert(boundaryNodeValues.end(), boundaryCalc.begin(), boundaryCalc.end());
        currPointer += ptrIncr;
    }

    auto RmOrder = allDirBoundary;
    std::sort(RmOrder.begin(), RmOrder.end());

    auto BcValsSorted = Utils::ReshuffleNodeVals(RmOrder, allDirBoundary, boundaryNodeValues);

    Eigen::Map<DvD> boundaryNodeValuesVec(BcValsSorted.data(), numBoundaryNodes, 1);
    return (DvD)boundaryNodeValuesVec;
}

SpD AssembleStokesSystem(QTM::QuadTreeMesh& mesh, 
                            SpD& dgPoissonMat,
                            SpD& QVMatrixX, 
                            SpD& QVMatrixY) {
    int nNodes = mesh.nNodes();
    int offset1 = nNodes; 
    int offset2 = 2*nNodes;
    int offset3 = 3*nNodes;
    int numElemNodes = mesh.numElemNodes;
    std::vector<double> gaussPoints = mesh.gaussPoints;
    int totalNNZ = 2*dgPoissonMat.nonZeros() + 2*QVMatrixX.nonZeros() + 2*QVMatrixY.nonZeros() + 2*nNodes;
    std::vector<Eigen::Triplet<double>> tripletList; tripletList.reserve(totalNNZ);

    // insert diffusive terms
    for (int k = 0; k < dgPoissonMat.outerSize(); ++k) {
        for (SpD::InnerIterator it(dgPoissonMat, k); it; ++it) {
            tripletList.emplace_back(it.row(), it.col(), it.value());
            tripletList.emplace_back(it.row()+offset1, it.col()+offset1, it.value());
        }
    }

    // insert pressure/divergence terms
    for (int k = 0; k < QVMatrixX.outerSize(); ++k) {
        for (SpD::InnerIterator it(QVMatrixX, k); it; ++it) {
            tripletList.emplace_back(it.row(), it.col()+offset2, it.value());
            tripletList.emplace_back(it.col()+offset2, it.row(), it.value());
        }
    }

    for (int k = 0; k < QVMatrixY.outerSize(); ++k) {
        for (SpD::InnerIterator it(QVMatrixY, k); it; ++it) {
            tripletList.emplace_back(it.row()+offset1, it.col()+offset2, it.value());
            tripletList.emplace_back(it.col()+offset2, it.row()+offset1, it.value());
        }
    }

    // add zero mean pressure constraint 
    int currDofIndx = offset2;
    auto weightMat = GenerateQuadWeights(gaussPoints, gaussPoints, mesh.deg+1, mesh.deg+1, numElemNodes);
    std::vector<double> weights(weightMat.diagonal().begin(), weightMat.diagonal().end());
    for (auto elm : mesh.leaves) {
        double area = (elm->width * elm->width);
        for (int i=0; i<numElemNodes; i++) {
            tripletList.emplace_back(offset3, currDofIndx, area * weights[i]);
            tripletList.emplace_back(currDofIndx, offset3, area * weights[i]);
            currDofIndx++;
        }
    }

    SpD mat(3*nNodes+1,3*nNodes+1);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

DvD IncompressibleStokesSolve(QTM::QuadTreeMesh& inputMesh,
                double rho,
                double mu,
                double penaltyParam,
                std::vector<std::string> source,
                std::vector<std::string> Ubcs,
                std::vector<std::string> Vbcs,
                std::vector<std::string> Pbcs,
                std::vector<std::string> Nbcs,
                std::vector<std::vector<int>>& boundaryNodes,
                std::vector<bool> velEss,
                std::vector<bool> pressEss,
                std::vector<bool> natBC) {
    std::vector<int> freeNodes = inputMesh.freeNodes;
    int nNodes = inputMesh.nNodes();
    int offset1 = 2*nNodes;

    for (auto i : velEss) {
        std::cout<<i<<std::endl;
    }

    std::cout<<"-------------"<<std::endl;

    for (auto i : pressEss) {
        std::cout<<i<<std::endl;
    }

    std::cout<<"Assembling stiffness matrix...";
    SpD KMatrix = StiffnessMatrix(inputMesh, mu); 
    std::cout<<" finished!"<<std::endl;
    
    std::cout<<"Assembling divergence matrix...";
    SpD PVMatrixX = PVMatrix(inputMesh, 1);
    SpD PVMatrixY = PVMatrix(inputMesh, 2); 
    std::cout<<" finished!"<<std::endl;

    std::cout<<"Assembling penalty matrix...";
    SpD PMatrix = PenaltyMatrix(inputMesh, mu, penaltyParam); 
    std::cout<<" finished!"<<std::endl;

    std::cout<<"Assembling flux matrix...";
    SpD SMatrix = FluxMatrix(inputMesh, mu); 
    std::cout<<" finished!"<<std::endl;

    std::cout<<"Assembling pressure flux matrix...";
    SpD PAMatrixX = PressureAvgMatrix(inputMesh, 1); // pressure average, velocity jump
    SpD PAMatrixY = PressureAvgMatrix(inputMesh, 2);

    SpD PJMatrixX = PressureJumpMatrix(inputMesh, 1); // pressure jump, velocity average
    SpD PJMatrixY = PressureJumpMatrix(inputMesh, 2); 
    std::cout<<" finished!"<<std::endl;

    std::cout<<"Assembling combined LHS matrix...";
    SpD SMatrixT = (SpD)(SMatrix.transpose()); 
    SpD dgPoissonMat = KMatrix + PMatrix - SMatrix - SMatrixT;

    // SpD mat1 = -PVMatrixX + PAMatrixX; // pressure gradient
    // SpD mat2 = -PVMatrixY + PAMatrixY;

    // SpD mat1T = PVMatrixX - PJMatrixX; // continuity
    // SpD mat2T = PVMatrixY - PJMatrixY;

    // SpD StiffnessMatrix = AssembleStokesSystem(inputMesh, dgPoissonMat, 
    //                         combinedPVMatrixX, combinedPVMatrixY,
    //                         combinedQUMatrixX, combinedQUMatrixY); std::cout<<" finished!"<<std::endl;

    SpD PJMatrixXT = PressureFluxMatrix(inputMesh, 1);
    SpD PJMatrixYT = PressureFluxMatrix(inputMesh, 2);

    SpD mat1 = -PVMatrixX + PAMatrixX; // pressure gradient
    SpD mat2 = -PVMatrixY + PAMatrixY;

    // SpD mat1 = -PVMatrixX + PJMatrixXT; // pressure gradient
    // SpD mat2 = -PVMatrixY + PJMatrixYT;

    // SpD mat1T = PVMatrixX - PJMatrixX; // continuity
    // SpD mat2T = PVMatrixY - PJMatrixY;

    // SpD mat1T = (SpD)(PVMatrixX.transpose());
    // SpD mat2T = (SpD)(PVMatrixY.transpose());

    // SpD mat1Int = -PVMatrixX + PJMatrixXT; // pressure gradient
    // SpD mat2Int = -PVMatrixY + PJMatrixYT;

    SpD mat1T = (SpD)(mat1.transpose()); // continuity
    SpD mat2T = (SpD)(mat2.transpose());

    SpD StiffnessMatrix = AssembleStokesSystem(inputMesh, dgPoissonMat, 
                            mat1, mat2,
                            mat1T, mat2T); std::cout<<" finished!"<<std::endl;

    std::cout<<"Assembling RHS vector...";
    DvD FMatrixX = AssembleFVec(inputMesh, 1.0, source[0]);
    DvD FMatrixY = AssembleFVec(inputMesh, 1.0, source[1]);
    DvD FMatrix = AssembleStokesSource(inputMesh, FMatrixX, FMatrixY);
    std::cout<<" finished!"<<std::endl;

    std::cout<<"Assembling boundary condition vector...";
    std::vector<int> allBoundaryNodes; // boundary dofs, *NOT* mesh nodes, where dirichlet bcs are applied 
    std::vector<std::vector<int>> VDN; // mesh nodes where velocity is fixed
    std::vector<std::vector<int>> PDN; // mesh nodes where pressure is fixed
    std::vector<int> freeNodesAll; freeNodesAll.reserve(3*freeNodes.size()); // interior and neumann dofs, *NOT* mesh nodes
    for (int i=0; i<boundaryNodes.size(); i++) {
        if (velEss[i]) {
            VDN.push_back(boundaryNodes[i]);
            continue;
        } else {
            for (int vn : boundaryNodes[i]) {
                freeNodesAll.push_back(vn);
                freeNodesAll.push_back(vn + nNodes);
            }
        }

        if (pressEss[i]) {
            PDN.push_back(boundaryNodes[i]);
        } else {
            for (int vn : boundaryNodes[i]) {
                freeNodesAll.push_back(vn + offset1);
            }
        }
    }
    DvD UDirichlet = EvalPartialSymbolicBoundaryCond(inputMesh, VDN, Ubcs, allBoundaryNodes, 0);
    DvD VDirichlet = EvalPartialSymbolicBoundaryCond(inputMesh, VDN, Vbcs, allBoundaryNodes, 1);
    DvD PDirichlet = EvalPartialSymbolicBoundaryCond(inputMesh, PDN, Pbcs, allBoundaryNodes, 2);
    DvD dirichletBoundaryVals(2*UDirichlet.rows() + PDirichlet.rows());
    dirichletBoundaryVals << UDirichlet, 
                            VDirichlet, 
                            PDirichlet;

    std::cout<<" finished!"<<std::endl;

    std::cout<<"Assembling null-orth matrix...";
    for (int i=0; i<3; i++) {
        int offset = i*nNodes;
        for (int fn : freeNodes) {
            freeNodesAll.push_back(offset + fn);
        }
    }
    std::sort(freeNodesAll.begin(), freeNodesAll.end());

    SpD nullSpace(3*nNodes+1, allBoundaryNodes.size());
    SpD columnSpace(3*nNodes+1, freeNodesAll.size());
    GetExtensionMatrices(inputMesh, allBoundaryNodes, freeNodesAll, nullSpace, columnSpace);
    std::cout<<" finished!"<<std::endl;

    std::cout<<"Solving system with "<<freeNodesAll.size()<<" degrees of freedom... ";
    DvD x = ComputeSolutionStationaryLinear(StiffnessMatrix, FMatrix, columnSpace, nullSpace, dirichletBoundaryVals);
    std::cout<<"System solved!"<<std::endl;
    std::cout<<"multiplier: "<<x(3*nNodes)<<std::endl;
    return x;
}

// DD IncompressibleNavierStokesSolve(QTM::QuadTreeMesh& inputMesh,
//                 double rho,
//                 double mu,
//                 std::vector<std::string> source,
//                 std::vector<std::string> Ubcs,
//                 std::vector<std::string> Vbcs,
//                 std::vector<std::string> Pbcs,
//                 double penaltyParam) {
//     auto boundaryNodes = inputMesh.boundaryNodes;
//     std::vector<int> freeNodes = inputMesh.freeNodes;
//     int nNodes = inputMesh.nNodes();

//     std::vector<int> allBoundaryNodes;
//     for (auto bNodes : boundaryNodes) {
//         allBoundaryNodes.insert(allBoundaryNodes.end(), bNodes.begin(), bNodes.end());
//     }

//     std::cout<<"Assembling stiffness matrix..."<<std::endl;
//     SpD KMatrix = StiffnessMatrix(inputMesh, mu);
//     std::cout<<"Assembling divergence matrix..."<<std::endl;
//     SpD QUMatrixX = QUMatrix(inputMesh, 1);
//     SpD QUMatrixY = QUMatrix(inputMesh, 2);
//     std::cout<<"Assembling penalty matrix..."<<std::endl;
//     SpD PMatrix = PenaltyMatrix(inputMesh, mu, penaltyParam);
//     std::cout<<"Assembling flux matrix..."<<std::endl;
//     SpD SMatrix = FluxMatrix(inputMesh, mu);
//     std::cout<<"Assembling pressure flux matrix..."<<std::endl;
//     SpD PJMatrixX = PressureJumpMatrix(inputMesh, 1);
//     SpD PJMatrixY = PressureJumpMatrix(inputMesh, 2);


//     std::cout<<"Assembling combined LHS matrix..."<<std::endl;
//     SpD StiffnessMatrix(3*nNodes, 3*nNodes);

//     std::cout<<"Assembling RHS vector..."<<std::endl;
//     DvD FMatrixX = AssembleFVec(inputMesh, 1.0, source[0]);
//     DvD FMatrixY = AssembleFVec(inputMesh, 1.0, source[1]);
//     SpD FMatrix(3*nNodes, 1);
//     std::cout<<"Assembling boundary condition vector..."<<std::endl;
//     DvD UDirichlet = EvalDirichletBoundaryCond(inputMesh, boundaryNodes, allBoundaryNodes, Ubcs);
//     DvD VDirichlet = EvalDirichletBoundaryCond(inputMesh, boundaryNodes, allBoundaryNodes, Vbcs);
//     DvD PDirichlet = EvalDirichletBoundaryCond(inputMesh, boundaryNodes, allBoundaryNodes, Pbcs);
//     DvD dirichletBoundaryVals(3*allBoundaryNodes.size());

//     SpD nullSpace(nNodes, allBoundaryNodes.size());
//     SpD columnSpace(nNodes, freeNodes.size());
//     SpD nullSpaceAll(3*nNodes, 3*allBoundaryNodes.size());
//     SpD columnSpaceAll(3*nNodes, 3*freeNodes.size());
//     std::cout<<"Assembling null-orth matrix..."<<std::endl;
//     GetExtensionMatrices(inputMesh, allBoundaryNodes, freeNodes, nullSpace, columnSpace);
// }