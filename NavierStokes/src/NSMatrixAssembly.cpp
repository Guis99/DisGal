#include "..\include\NSMatrixAssembly.hpp"

SpD QVMatrix(QTM::QuadTreeMesh& mesh, int diffDir) {
    // Oriented as V-P. For Q-U, take transpose
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

    DD topRowVec = DD::Zero(1,numNodes); topRowVec(0,0) = 1;
    DD bottomRowVec = DD::Zero(1,numNodes); bottomRowVec(0,deg) = 1;
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
    return mat;
}