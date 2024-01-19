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


DD StiffnessMatrix(QTM::QuadTreeMesh mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);

    // Generate derivative matrices
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
        localElemMat = combinedX.transpose()*coeffMat*weightMat*combinedX +
                        combinedY.transpose()*coeffMat*weightMat*combinedY;
 
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

DD StiffnessMatrix(QTM::QuadTreeMesh mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);

    // Generate derivative matrices
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
        localElemMat = combinedX.transpose()*coeffMat*weightMat*combinedX +
                        combinedY.transpose()*coeffMat*weightMat*combinedY;
 
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

DD StiffnessMatrix(QTM::QuadTreeMesh mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);

    // Generate derivative matrices
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
        localElemMat = combinedX.transpose()*coeffMat*weightMat*combinedX +
                        combinedY.transpose()*coeffMat*weightMat*combinedY;
 
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

DD PenaltyMatrix(QTM::QuadTreeMesh mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
}

DD FluxMatrix(QTM::QuadTreeMesh mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
}