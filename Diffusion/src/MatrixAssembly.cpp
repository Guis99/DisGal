#include "..\include\MatrixAssembly.hpp"

DD StiffnessMatrix(QTM::QuadTreeMesh mesh, double k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nNodes = mesh.numLeaves;
    std::vector<double> gaussPointsX = Utils::genGaussPoints(deg);

    // Generate derivative matrices
    std::vector<double> AInitializer; 
    AInitializer.reserve(numElemNodes * numElemNodes);

    double *AInitIdx = AInitializer.data(); 

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPointsX, gaussPointsX);
        std::copy(xPartials.begin(), xPartials.end(), AInitIdx);
        AInitIdx += numNodes;
    }
    

    // map derivative values to matrix
    Eigen::Map<DD> A(AInitializer.data(), numNodes, numNodes); 

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