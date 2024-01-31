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

DD PenaltyMatrix(QTM::QuadTreeMesh mesh, double k, double alpha) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;

    // get quad weights in 1D
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
    std::vector<double> integX = Utils::integrateLagrange(gaussPoints);

    // Convert 1D quad weights to diag matrix
    Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    DD quadWeights1D = (DD)integXMat.asDiagonal();

    auto leaves = mesh.GetAllCells();
    std::vector<QTM::Direction> directions = {QTM::Direction::N, QTM::Direction::E, 
                                                QTM::Direction::S, QTM::Direction::W};
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
            elemNodes = mesh.GetBoundaryNodes(dir, elm->CID);
            for (auto neighbor : neighbors) { 
                if (neighbor == nullptr) { // case curr boundary is exterior
                    continue;
                } else if (neighbor->CID > elm->CID) { // wrong side of boundary
                    
                } else { // neighbor exists 
                    neighborNodes = mesh.GetBoundaryNodes((QTM::Direction)((dir+2)%4), neighbor->CID);
                    // calculate penalty param
                    a = alpha/std::min(elm->width, neighbor->width)/2;
                    // calculate jump matrix
                    DD jumpMatrix(elemNodes.size() + neighborNodes.size(), elemNodes.size());


                    // TODO

                    DD jumpMatrixT = (DD)jumpMatrix.transpose();


                    boundaryNodes.reserve(elemNodes.size() + neighborNodes.size());
                    boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                    boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                    DD localElemMat(elemNodes.size() + neighborNodes.size(), elemNodes.size() + neighborNodes.size());
                    localElemMat = (DD)(jumpMatrix * elm->width*quadWeights1D * jumpMatrixT);

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
}

DD FluxMatrix(QTM::QuadTreeMesh mesh, double k) {
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
    for (auto &elm : leaves) {
        // get neighbors
        for (auto dir : directions) {
            cellSign = edgeSign[dir]; neighborSign = -cellSign;
            neighbors = mesh.GetCellNeighbors(dir, elm->CID);
            elemNodes = mesh.GetBoundaryNodes(dir, elm->CID);
            for (auto neighbor : neighbors) { 
                if (neighbor == nullptr) { // case curr boundary is exterior
                    continue;
                } else if (neighbor->CID > elm->CID) { // wrong side of boundary
                    
                } else { // neighbor exists 
                    neighborNodes = mesh.GetBoundaryNodes((QTM::Direction)((dir+2)%4), neighbor->CID);
                    // calculate jump matrix
                    DD jumpMatrix(elemNodes.size() + neighborNodes.size(), elemNodes.size());

                    // TODO


                    // calculate flux matrix
                    DD fluxMatrixX(elemNodes.size(), elemNodes.size() + neighborNodes.size());
                    DD fluxMatrixY(elemNodes.size(), elemNodes.size() + neighborNodes.size());

                    // TODO

                    boundaryNodes.reserve(elemNodes.size() + neighborNodes.size());
                    boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                    boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                    // assemble local matrix
                    DD localElemMat(elemNodes.size() + neighborNodes.size(), elemNodes.size() + neighborNodes.size());
                    localElemMat = (DD)(jumpMatrix * elm->width*quadWeights1D * fluxMatrixX + 
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
    DvD sourceVec = f * weightMat.diagonal();

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

