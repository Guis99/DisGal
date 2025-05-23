#include "../include/MatrixAssembly.hpp"
#include "../../Dependencies/Eigen/SVD"

DD PMA::GenerateQuadWeights(std::vector<Real_b> &gpX, std::vector<Real_b> &gpY, int numXNodes, int numYNodes, int numElemNodes) {
    // Generate quadrature weight matrices
    // Get integrals of basis functions
    std::vector<Real_b> integX = Utils::integrateLagrange(gpX);
    std::vector<Real_b> integY = Utils::integrateLagrange(gpY);

    Eigen::Map<DvD> integXMat(integX.data(),numXNodes);
    Eigen::Map<DvD> integYMat(integY.data(),numYNodes);

    DD weightMat(numElemNodes, numElemNodes);
    weightMat << Eigen::kroneckerProduct((DD)integXMat.asDiagonal(),(DD)integYMat.asDiagonal());

    return weightMat;
}


SpD PMA::StiffnessMatrix(QTM::QuadTreeMesh& mesh, Real_b k) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();
    
    std::vector<Real_b> gaussPoints = Utils::genGaussPoints(deg);

    // Generate derivative matrix
    std::vector<Real_b> AInitializer; 
    AInitializer.reserve(numElemNodes * numElemNodes);

    Real_b *AInitIdx = AInitializer.data(); 

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<Real_b> xPartials = Utils::numDeriv(.00001, k, gaussPoints, gaussPoints);
        std::copy(xPartials.begin(), xPartials.end(), AInitIdx);
        AInitIdx += numNodes;
    }
    // map derivative values to matrix
    Eigen::Map<DD> A(AInitializer.data(), numNodes, numNodes); 
    // Generate quadrature weight matrices
    DD weightMat = PMA::GenerateQuadWeights(gaussPoints, gaussPoints, numNodes, numNodes, numElemNodes);
    // Generate mass matrices
    DD B; B.setIdentity(numNodes, numNodes);
    SpD coeffMat(numElemNodes, numElemNodes);
    coeffMat.setIdentity(); coeffMat *= k;
    
    // Get element-wise matrix intermediates
    DD combinedX(numElemNodes, numElemNodes);
    combinedX << Eigen::kroneckerProduct(B, A);
    DD combinedY(numElemNodes, numElemNodes);
    combinedY << Eigen::kroneckerProduct(A, B);
    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<Real_b>> tripletList;
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
    mat.makeCompressed();
    return mat;
}

SpD PMA::PenaltyMatrix(QTM::QuadTreeMesh& mesh, Real_b k, Real_b alpha, PMA::quadUtils& package) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = mesh.numElemNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();

    //  load package data
    DD weightMat = package.weightMat;
    DD quadWeights1D = package.quadWeights1D;

    DD nodalVals = package.nodalVals;

    DD combinedX = package.combinedX;
    DD combinedY = package.combinedY;

    std::array<DD,4> splitCellVals = package.splitCellVals;

    std::array<DD,4> splitCellGradsX = package.splitCellGradsX;
    std::array<DD,4> splitCellGradsY = package.splitCellGradsY;

    std::vector<QTM::Direction> directions = package.directions;        
    std::vector<QTM::Direction> oppdirs = package.oppdirs;
    std::vector<std::vector<int>> localNodes = package.localNodes;

    // basis func to node mapping
    DD B; B.setIdentity(numElemNodes, numElemNodes);

    // get basis func vals for split cell quad
    std::vector<Real_b> BhInitializer; 
    BhInitializer.reserve((mesh.halfGaussPoints.size()) * numNodes);

    // get unit vecs
    DvD topVec = DvD::Zero(numNodes); topVec(0) = 1;
    DvD bottomVec = DvD::Zero(numNodes); bottomVec(deg) = 1;

    // index vectors for split cell gauss points
    std::vector<int> frontIdxs(numNodes, 0);  
    std::vector<int> backIdxs(numNodes, 0);

    for (int i=0; i<numNodes; i++) {
        frontIdxs[i] = i;
        backIdxs[i] = i+deg;
    }

    std::vector<int> splitIdx[2] = { frontIdxs, backIdxs };

    auto leaves = mesh.leaves;
    Real_b a; // penalty parameters

    std::vector<int> boundaryNodes;
    std::vector<int> neighborNodes;
    std::vector<int> neighborLocals;
    std::vector<int> elemNodes;
    std::vector<int> elemLocals;
    std::vector<std::shared_ptr<QTM::Cell>> neighbors;
    std::vector<Eigen::Triplet<Real_b>> tripletList; tripletList.reserve(nNodes);

    for (auto &elm : leaves) {
        elemNodes = mesh.GetGlobalElemNodes(elm->CID);
        elemLocals = mesh.GetTrimmedLocalNodes(elm->CID, elemNodes);
        // std::cout<<"here5"<<std::endl;
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
            // std::cout<<"here6"<<std::endl;
            // if (neighbors.size() == 0) {
            //     // neighborNodes = {};
            //     // neighborLocals = {};
            //     continue;
            // }
            for (int NI = 0; NI < neighbors.size(); NI++) { 
                auto neighbor = neighbors[NI];
                Real_b jac;
                if (!neighbor) { // skip exterior edges
                    // std::cout<<"here7"<<std::endl;
                    continue;
                }
                
                // calculate penalty param
                else if (neighbor->CID < elm->CID || elm->level < neighbor->level) { // case appropriate neighbor exists
                // std::cout<<"here8"<<std::endl;
                    jac = std::min(elm->width, neighbor->width);
                    neighborNodes = mesh.GetGlobalElemNodes(neighbor->CID);
                    neighborLocals = mesh.GetTrimmedLocalNodes(neighbor->CID, neighborNodes);
                    // std::cout<<"here9"<<std::endl;
                } else { 
                    continue;
                }
                a = alpha/jac/2;
                // calculate jump matrix
                DD topJump;
                // std::cout<<nodalVals.cols()<<", "<<nodalVals.rows()<<std::endl;
                // Utils::printVec(neighborLocals);
                // Utils::printVec(localNodes[oppdir]);
                DD bottomJump = -nodalVals(neighborLocals, localNodes[oppdir]);
                // std::cout<<"here10"<<std::endl;
                if (elm->level == neighbor->level) {
                    topJump = nodalVals(elemLocals, localNodes[dir]);
                } else {
                    topJump = splitCellVals[dir](elemLocals, splitIdx[NI]);
                }
                // std::cout<<"here11"<<std::endl;
                DD jumpMatrix(elemNodes.size() + neighborNodes.size(), numNodes);
                // std::cout<<topJump.rows()<<", "<<topJump.cols()<<std::endl;
                // std::cout<<bottomJump.rows()<<", "<<bottomJump.cols()<<std::endl;
                // std::cout<<elemNodes.size()<<", "<<neighborNodes.size()<<", "<<neighborLocals.size()<<std::endl;
                // std::cout<<"-----------"<<std::endl;
                jumpMatrix << topJump, bottomJump;

                DD jumpMatrixT = (DD)(jumpMatrix.transpose());
                // std::cout<<"here12"<<std::endl;
                DD localElemMat = (DD)(a * jumpMatrix * jac*quadWeights1D * jumpMatrixT);

                boundaryNodes.reserve(elemNodes.size() + neighborNodes.size());
                boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                // std::cout<<"here13"<<std::endl;

                for (int j=0; j<boundaryNodes.size(); j++) {
                    for (int i=0; i<boundaryNodes.size(); i++) {
                        tripletList.emplace_back(boundaryNodes[i],boundaryNodes[j],localElemMat(i,j));
                    }
                }    
                // std::cout<<"here14"<<std::endl;
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

SpD PMA::FluxMatrix(QTM::QuadTreeMesh& mesh, Real_b k, PMA::quadUtils& package) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = mesh.numElemNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();
    std::vector<Real_b> gaussPoints = Utils::genGaussPoints(deg);

    // get quad weights in 1D
    std::vector<Real_b> integX = Utils::integrateLagrange(gaussPoints);

    // load package data
    DD quadWeights1D = package.quadWeights1D;
    DD nodalValues = package.nodalVals;

    DD combinedX = package.combinedX;
    DD combinedY = package.combinedY;
    
    std::array<DD, 4> splitCellVals = package.splitCellVals;
    std::array<DD, 4> splitCellGradsX = package.splitCellGradsX;
    std::array<DD, 4> splitCellGradsY = package.splitCellGradsY;

    std::vector<QTM::Direction> directions = package.directions;        
    std::vector<QTM::Direction> oppdirs = package.oppdirs;
    std::vector<std::vector<int>> localNodes = package.localNodes;

    // Convert 1D quad weights to diag matrix
    // Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    // DD quadWeights1D = (DD)integXMat.asDiagonal();

    // // basis func to node mapping
    // DD B; B.setIdentity(numElemNodes, numElemNodes);
    // DD Bs; Bs.setIdentity(numNodes, numNodes);

    // // Generate derivative matrix
    // std::vector<Real_b> AInitializer; 
    // AInitializer.reserve(numNodes * numNodes);

    // // Generate derivatives for each basis function, copy to full array
    // for (int k=0; k<numNodes; k++) { 
    //     std::vector<Real_b> xPartials = Utils::numDeriv(.00001, k, gaussPoints, gaussPoints);
    //     AInitializer.insert(AInitializer.end(), xPartials.begin(), xPartials.end());
    // }
    
    // // map derivative values to matrix
    // Eigen::Map<DD> A(AInitializer.data(), numNodes, numNodes); 
    // // Get element-wise matrix intermediates
    // DD combinedX(numElemNodes, numElemNodes);
    // combinedX << Eigen::kroneckerProduct(Bs, A);
    // DD combinedY(numElemNodes, numElemNodes);
    // combinedY << Eigen::kroneckerProduct(A, Bs);

    // // get basis func vals for split cell quad
    // std::vector<Real_b> BhInitializer; 
    // BhInitializer.reserve((mesh.halfGaussPoints.size()) * numNodes);

    // // get unit vecs
    // DvD topVec = DvD::Zero(numNodes); topVec(0) = 1;
    // DvD bottomVec = DvD::Zero(numNodes); bottomVec(deg) = 1;

    // DD topRowVec = DD::Zero(1,numNodes); topRowVec(0,0) = 1;
    // DD bottomRowVec = DD::Zero(1,numNodes); bottomRowVec(0,deg) = 1;

    // DD fullRowVec = DD::Ones(1,numNodes);

    // // Generate values for each basis function, copy to full array
    // for (int k=0; k<numNodes; k++) { 
    //     std::vector<Real_b> xVals = Utils::evalLagrangeInterp(k, mesh.halfGaussPoints, mesh.gaussPoints);
    //     BhInitializer.insert(BhInitializer.end(), xVals.begin(), xVals.end());
    // }

    // // map basis values to matrix
    // Eigen::Map<DD> BhT(BhInitializer.data(), mesh.halfGaussPoints.size(), numNodes); // eval points are traversed first
    // DD Bh = (DD)(BhT.transpose());

    // std::array<DD,4> splitCellVals;

    // DD splitCellPlaceholder(numElemNodes, mesh.halfGaussPoints.size());
    // splitCellPlaceholder << Eigen::kroneckerProduct(bottomVec, Bh); splitCellVals[0] = splitCellPlaceholder;
    // splitCellPlaceholder << Eigen::kroneckerProduct(Bh, bottomVec); splitCellVals[1] = splitCellPlaceholder;
    // splitCellPlaceholder << Eigen::kroneckerProduct(topVec, Bh); splitCellVals[2] = splitCellPlaceholder;
    // splitCellPlaceholder << Eigen::kroneckerProduct(Bh, topVec); splitCellVals[3] = splitCellPlaceholder;

    // // get basis func gradients for split cell quad
    // std::vector<Real_b> AhInitializer; 
    // AhInitializer.reserve((mesh.halfGaussPoints.size()) * numNodes);

    // // Generate derivatives for each basis function, copy to full array
    // for (int k=0; k<numNodes; k++) { 
    //     std::vector<Real_b> xPartials = Utils::numDeriv(.00001, k, mesh.halfGaussPoints, gaussPoints);
    //     AhInitializer.insert(AhInitializer.end(), xPartials.begin(), xPartials.end());
    // }
    
    // // map derivative values to matrix
    // Eigen::Map<DD> Ah(AhInitializer.data(), mesh.halfGaussPoints.size(), numNodes); 

    // DD topGrads = A(0, Eigen::all);
    // DD bottomGrads = A(numNodes-1, Eigen::all);

    // std::array<DD,4> splitCellGradsX;
    
    // DD splitCellGradXPlaceholder(mesh.halfGaussPoints.size(), numElemNodes);
    // splitCellGradXPlaceholder << Eigen::kroneckerProduct(bottomRowVec, Ah); splitCellGradsX[0] = splitCellGradXPlaceholder;
    // splitCellGradXPlaceholder << Eigen::kroneckerProduct(BhT, bottomGrads); splitCellGradsX[1] = splitCellGradXPlaceholder;
    // splitCellGradXPlaceholder << Eigen::kroneckerProduct(topRowVec, Ah); splitCellGradsX[2] = splitCellGradXPlaceholder;
    // splitCellGradXPlaceholder << Eigen::kroneckerProduct(BhT, topGrads); splitCellGradsX[3] = splitCellGradXPlaceholder;
    // std::array<DD,4> splitCellGradsY;
    // DD splitCellGradYPlaceholder(mesh.halfGaussPoints.size(), numElemNodes); 
    // splitCellGradYPlaceholder << Eigen::kroneckerProduct(bottomGrads, BhT); splitCellGradsY[0] = splitCellGradYPlaceholder;
    // splitCellGradYPlaceholder << Eigen::kroneckerProduct(Ah, bottomRowVec); splitCellGradsY[1] = splitCellGradYPlaceholder;
    // splitCellGradYPlaceholder << Eigen::kroneckerProduct(topGrads, BhT); splitCellGradsY[2] = splitCellGradYPlaceholder;
    // splitCellGradYPlaceholder << Eigen::kroneckerProduct(Ah, topRowVec); splitCellGradsY[3] = splitCellGradYPlaceholder;

    // index vectors for split cell gauss points
    std::vector<int> frontIdxs(numNodes, 0);  
    std::vector<int> backIdxs(numNodes, 0);

    for (int i=0; i<numNodes; i++) {
        frontIdxs[i] = i;
        backIdxs[i] = i+deg;
    }

    std::vector<int> splitIdx[2] = { frontIdxs, backIdxs };

    Real_b normalX[4] = {0,1,0,-1};
    Real_b normalY[4] = {1,0,-1,0};

    std::vector<int> boundaryNodes;
    std::vector<int> neighborNodes;
    std::vector<int> neighborLocals;
    std::vector<int> elemNodes;
    std::vector<int> elemLocals;
    std::vector<std::shared_ptr<QTM::Cell>> neighbors;
    std::vector<Eigen::Triplet<Real_b>> tripletList; tripletList.reserve(nNodes);
    auto leaves = mesh.GetAllCells();
                                                
    for (auto &elm : leaves) {
        elemNodes = mesh.GetGlobalElemNodes(elm->CID);
        elemLocals = mesh.GetTrimmedLocalNodes(elm->CID, elemNodes);
        // get neighbors
        for (auto dir : directions) {
            neighbors = mesh.GetCellNeighbors(dir, elm->CID);
            QTM::Direction oppdir = oppdirs[dir];
            for (int NI = 0; NI < neighbors.size(); NI++) { 
                auto neighbor = neighbors[NI];
                Real_b jac;
                Real_b fac;
                if (!neighbor) { // skip exterior edges
                    continue;
                }

                else if (neighbor->CID < elm->CID || elm->level < neighbor->level) { // case appropriate neighbor exists
                    jac = std::min(elm->width, neighbor->width);
                    neighborNodes = mesh.GetGlobalElemNodes(neighbor->CID);
                    neighborLocals = mesh.GetTrimmedLocalNodes(neighbor->CID, neighborNodes);
                    fac = .5;
                } else { 
                    continue;
                }
                // jump matrix setup
                DD topJump;
                DD bottomJump = -nodalValues(neighborLocals, localNodes[oppdir]);

                // flux matrix setup
                DD topGradX;
                DD topGradY;
                DD bottomGradX = normalX[oppdir] * combinedX(localNodes[oppdir], neighborLocals);
                DD bottomGradY = normalY[oppdir] * combinedY(localNodes[oppdir], neighborLocals);

                if (elm->level == neighbor->level) {
                    topJump = nodalValues(elemLocals, localNodes[dir]);
                    topGradX = normalX[dir] * combinedX(localNodes[dir], elemLocals);
                    topGradY = normalY[dir] * combinedY(localNodes[dir], elemLocals);
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
                DD localElemMat = (DD)(fac * jumpMatrix * quadWeights1D * fluxMatrixX + 
                                    fac * jumpMatrix * quadWeights1D * fluxMatrixY);

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

DvD PMA::AssembleFVec(QTM::QuadTreeMesh& mesh, Real_b f, std::string evalStr) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();
    std::vector<Real_b> gaussPoints = Utils::genGaussPoints(deg);

    // Generate quadrature weight matrices
    DD weightMat = PMA::GenerateQuadWeights(gaussPoints, gaussPoints, numNodes, numNodes, numElemNodes);

    // Turn weight mat int vector and mult. by source since diagonal
    // DvD sourceVec = f * weightMat.diagonal();

    std::vector<std::array<Real_b,2>> allNodesPos = mesh.AllNodePos();
    std::array<Real_b,2>* startpoint = allNodesPos.data(); auto allocSize = allNodesPos.size();
    std::vector<Real_b> fEval = Utils::EvalSymbolicBC(startpoint, allocSize, evalStr);

    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<Real_b>> tripletList;
    tripletList.reserve(nElements * numElemNodes);

    // Integrate over all elements
    DvD localElemMat(numElemNodes);
    auto leaves = mesh.GetAllCells();
    for (auto &elm : leaves) {
        Real_b jac = elm->width; // Jacobian factors
        // calculate local matrix
        auto nodes = elm->nodes;
        std::vector<Real_b> collectSourceVals; collectSourceVals.reserve(numElemNodes);
        collectSourceVals.insert(collectSourceVals.begin(), fEval.begin()+nodes[0], fEval.begin()+nodes[1]+1);
        // for (int i=nodes[0]; i<=nodes[1]; i++) {
        //     collectSourceVals.push_back(fEval[i]);
        // }
        Eigen::Map<DvD> sourceVector(collectSourceVals.data(), numElemNodes);

        localElemMat = weightMat*jac*jac*sourceVector;
        // Get nodes in element
        
        // Generate i,j,v triplets
        for (int i=0; i<numElemNodes; i++) {
            tripletList.emplace_back(nodes[0]+i, 0, localElemMat(i));
            // out[nodes[0]+i] += localElemMat[i];
        }  
    }

    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,1);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    mat.makeCompressed();
    DvD out = (DvD)mat;
    return out;
}

std::vector<Real_b> PMA::ComputeResiduals(QTM::QuadTreeMesh& mesh, DvD& solution, SpD& source) {
    std::vector<Real_b> out;

    // get quadrature weights

    // second derivatives


    for (auto leaf : mesh.leaves) {

    }
    return out;
}

// std::vector<std::shared_ptr<QTM::Cell>> PMA::TestResiduals(DvD& solution, QTM::QuadTreeMesh& mesh, Real_b residualLimit) {
//     std::vector<std::shared_ptr<QTM::Cell>> out; out.reserve(mesh.leaves.size());
//     for (auto leaf : mesh.leaves) {
//         Real_b res = ComputeResidual();
//         if (leaf->level < 5 && ComputeResidual() > residualLimit) {
//              out.push_back(leaf);
//         }
//     }

//     return out;
// }

DvD PMA::EvalSymbolicBoundaryCond(QTM::QuadTreeMesh& inputMesh, std::vector<std::vector<int>>& boundaryNodes, std::vector<int>& allBoundaryNodes, std::vector<std::string>& strs) {
    std::vector<std::array<Real_b,2>> boundaryNodePos;
    for (auto bNodes : boundaryNodes) {
        std::vector<std::array<Real_b,2>> nodePos = inputMesh.GetNodePos(bNodes);
        boundaryNodePos.insert(boundaryNodePos.end(), nodePos.begin(), nodePos.end());
    }
    int numBoundaryNodes = allBoundaryNodes.size();

    // Boundary nodes are given in clockwise order, not column-major order
    std::vector<Real_b> boundaryNodeValues; 
    boundaryNodeValues.reserve(numBoundaryNodes);

    std::vector<Real_b> boundaryCalc;
    std::array<Real_b,2> *currPointer = boundaryNodePos.data();
    int ptrIncr;
    std::string prompt;

    for (int i=0; i<boundaryNodes.size(); i++) {
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

void PMA::GetExtensionMatrices(QTM::QuadTreeMesh& inputMesh,
                                        std::vector<int>& boundaryNodes, 
                                        std::vector<int>& freeNodes,
                                        SpD& nullSpace,
                                        SpD& columnSpace) {
    int nNodes = inputMesh.nNodes();

    std::sort(boundaryNodes.begin(), boundaryNodes.end());
    std::sort(freeNodes.begin(), freeNodes.end());

    std::vector<Eigen::Triplet<Real_b>> tripletListNS;
    std::vector<Eigen::Triplet<Real_b>> tripletListCS;
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

// void PMA::FindConditionNumber(SpD& mat) {
//     using namespace Eigen;
//     Eigen::JacobiSVD<SpD> svd(mat);

//     Real_b cond = svd.singularValues()(0) 
//     / svd.singularValues()(svd.singularValues().size()-1);

//     std::cout << "Condition number: " << cond << std::endl;
// }

void PMA::FindRank(SpD& mat) {
    Eigen::SparseQR<Eigen::SparseMatrix<Real_b>, Eigen::COLAMDOrdering<int>> qr;
    qr.compute(mat);

    // Count the number of non-zero (or significant) diagonal elements in R
    int rank = 0;
    Real_b tol = 1e-10; // Tolerance for considering an element as non-zero
    for(int i = 0; i < qr.matrixR().rows(); ++i) {
        if (std::abs(qr.matrixR().coeff(i, i)) > tol) {
            rank++;
        }
    }

    std::cout << "Rank: " << rank << std::endl;
    std::cout<<"Size: "<<mat.rows()<<std::endl;
}

DvD PMA::ComputeSolutionStationaryLinear(SpD& StiffnessMatrix, DvD& fVec, SpD& columnSpace, SpD& nullSpace, DvD& boundaryVals) {
    // Eliminate rows and columns corr. to boundary nodes
    SpD columnSpaceT = (SpD)(columnSpace.transpose());
    SpD nullSpaceT = (SpD)(nullSpace.transpose());

    // Eliminate boundary rows and boundary columns
    SpD A11 = columnSpaceT * StiffnessMatrix * columnSpace;
    DEBUG_PRINT("Stokes stiffness matrix: ", A11);
    // Eliminate boundary rows and free columns
    SpD A12 = columnSpaceT * StiffnessMatrix * nullSpace;
    // Eliminate boundary rows
    DvD F11 = columnSpaceT * fVec;

    using namespace Eigen;

    // SparseLU<SpD,COLAMDOrdering<int>> solver;
    // solver.analyzePattern(A11);
    // FindRank(A11);
    // solver.factorize(A11);
    // DvD x = solver.solve(F11 - A12 * boundaryVals); 

    ConjugateGradient<SpD, Lower|Upper> cg;
    cg.compute(A11);
    DvD x = cg.solve(F11 - A12 * boundaryVals);

    x = columnSpace * x + nullSpace * boundaryVals;
    return x;
}

DvD PMA::PoissonSolve(QTM::QuadTreeMesh& inputMesh,
                Real_b c,
                Real_b k,
                std::string source,
                std::vector<std::string> bcs,
                Real_b penaltyParam) {
    
    auto boundaryNodes = inputMesh.boundaryNodes;
    std::vector<int> freeNodes = inputMesh.freeNodes;
    int nNodes = inputMesh.nNodes();

    std::vector<int> allBoundaryNodes;
    for (auto bNodes : boundaryNodes) {
        allBoundaryNodes.insert(allBoundaryNodes.end(), bNodes.begin(), bNodes.end());
    }

    
    PMA::quadUtils package = PMA::GenerateAssemblyPackage(inputMesh);

    std::cout<<"Assembling stiffness matrix"<<std::endl;
    SpD KMatrix = PMA::StiffnessMatrix(inputMesh, k);
    std::cout<<"Assembling penalty matrix"<<std::endl;
    SpD PMatrix = PMA::PenaltyMatrix(inputMesh, k, penaltyParam, package);
    std::cout<<"Assembling flux matrix"<<std::endl;
    SpD SMatrix = PMA::FluxMatrix(inputMesh, k, package);
    SpD SMatrixT = (SpD)(SMatrix.transpose());
    std::cout<<"Assembling overall system matrix"<<std::endl;
    SpD StiffnessMatrix = KMatrix - SMatrix - SMatrixT + PMatrix;

    std::cout<<"Assembling RHS vector"<<std::endl;
    DvD FMatrix = PMA::AssembleFVec(inputMesh, 1.0, source);
    std::cout<<"Assembling boundary condition vector"<<std::endl;
    DvD boundaryVals = PMA::EvalSymbolicBoundaryCond(inputMesh, boundaryNodes, allBoundaryNodes, bcs);

    SpD nullSpace(nNodes, allBoundaryNodes.size());
    SpD columnSpace(nNodes, freeNodes.size());
    std::cout<<"Assembling null-orth matrix"<<std::endl;
    PMA::GetExtensionMatrices(inputMesh, allBoundaryNodes, freeNodes, nullSpace, columnSpace);

    SpD columnSpaceT = (SpD)(columnSpace.transpose());
    SpD nullSpaceT = (SpD)(nullSpace.transpose());

    // SpD B11 = columnSpaceT * PMatrix * columnSpace;
    // // Eliminate boundary rows and free columns
    // SpD B12 = columnSpaceT * PMatrix * nullSpace;


    // SpD A11 = columnSpaceT * SMatrix * columnSpace;
    // // Eliminate boundary rows and free columns
    // SpD A12 = columnSpaceT * SMatrix * nullSpace;

    // SpD C11 = columnSpaceT * SMatrixT * columnSpace;
    // // Eliminate boundary rows and free columns
    // SpD C12 = columnSpaceT * SMatrixT * nullSpace;

    // std::cout<<"dropped nodes: \n";
    // for (int i : allBoundaryNodes) {
    //     std::cout<<i<<std::endl;
    // }

    // std::cout<<"pre-reduce (flux): \n"<<SMatrix<<std::endl;
    // std::cout<<"after reduce (flux): \n"<<A11<<std::endl;
    // std::cout<<"right side (flux): \n"<<A12<<std::endl;

    // std::cout<<"--------------------"<<std::endl;

    // std::cout<<"pre-reduce (flux.T): \n"<<SMatrixT<<std::endl;
    // std::cout<<"after reduce (flux.T): \n"<<C11<<std::endl;
    // std::cout<<"right side (flux.T): \n"<<C12<<std::endl;

    // std::cout<<"--------------------"<<std::endl;

    // std::cout<<"pre-reduce (penalty): \n"<<PMatrix<<std::endl;
    // std::cout<<"after reduce (penalty): \n"<<B11<<std::endl;
    // std::cout<<"right side (penalty): \n"<<B12<<std::endl;

    std::cout<<"Solving system with "<<freeNodes.size()<<" nodes"<<std::endl;
    DvD x = PMA::ComputeSolutionStationaryLinear(StiffnessMatrix, FMatrix, columnSpace, nullSpace, boundaryVals);
    std::cout<<"System solved!"<<std::endl;
    return x;
}


PMA::quadUtils PMA::GenerateAssemblyPackage(QTM::QuadTreeMesh& mesh) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = mesh.numElemNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();
    std::vector<Real_b> gaussPoints = Utils::genGaussPoints(deg);

    // get quad weights in 1D
    std::vector<Real_b> integX = Utils::integrateLagrange(gaussPoints);

    // Convert 1D quad weights to diag matrix
    Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    DD quadWeights1D = (DD)integXMat.asDiagonal();

    DD weightMat = PMA::GenerateQuadWeights(gaussPoints, gaussPoints, numNodes, numNodes, numElemNodes);

    // basis func to node mapping
    DD B; B.setIdentity(numNodes, numNodes);
    DD Bs; Bs.setIdentity(numNodes, numNodes);

    // Generate derivative matrix
    std::vector<Real_b> AInitializer; 
    AInitializer.reserve(numNodes * numNodes);

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<Real_b> xPartials = Utils::numDeriv(.00001, k, gaussPoints, gaussPoints);
        AInitializer.insert(AInitializer.end(), xPartials.begin(), xPartials.end());
    }
    
    // map derivative values to matrix
    Eigen::Map<DD> A(AInitializer.data(), numNodes, numNodes); 
    // Gradients of basis function on all element nodes
    DD combinedX(numElemNodes, numElemNodes);
    combinedX << Eigen::kroneckerProduct(B, A);
    DD combinedY(numElemNodes, numElemNodes);
    combinedY << Eigen::kroneckerProduct(A, B);

    // get basis func vals for split cell quad
    std::vector<Real_b> BhInitializer; 
    BhInitializer.reserve((mesh.halfGaussPoints.size()) * numNodes);

    // get unit vecs
    DvD topVec = DvD::Zero(numNodes); topVec(0) = 1;
    DvD bottomVec = DvD::Zero(numNodes); bottomVec(deg) = 1;

    DD topRowVec = DD::Zero(1,numNodes); topRowVec(0,0) = 1;
    DD bottomRowVec = DD::Zero(1,numNodes); bottomRowVec(0,deg) = 1;

    DD fullRowVec = DD::Ones(1,numNodes);

    // Generate values for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<Real_b> xVals = Utils::evalLagrangeInterp(k, mesh.halfGaussPoints, mesh.gaussPoints);
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
    std::vector<Real_b> AhInitializer; 
    AhInitializer.reserve((mesh.halfGaussPoints.size()) * numNodes);

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<Real_b> xPartials = Utils::numDeriv(.00001, k, mesh.halfGaussPoints, gaussPoints);
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

    Real_b normalX[4] = {0,1,0,-1};
    Real_b normalY[4] = {1,0,-1,0};

    std::vector<QTM::Direction> directions = {QTM::Direction::N, QTM::Direction::E, 
                                                QTM::Direction::S, QTM::Direction::W};

                                                
    std::vector<QTM::Direction> oppdirs = {QTM::Direction::S, QTM::Direction::W, 
                                            QTM::Direction::N, QTM::Direction::E};

    std::vector<std::vector<int>> localNodes = {mesh.GetLocalBoundaryNodes(QTM::Direction::N),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::E),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::S),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::W)};

    DD nodalVals; nodalVals.setIdentity(numElemNodes, numElemNodes);
    
    // populate package with arrays
    PMA::quadUtils package;

    package.weightMat = weightMat;
    package.quadWeights1D = quadWeights1D;

    package.combinedX = combinedX;
    package.combinedY = combinedY;

    package.nodalVals = nodalVals;
    package.splitCellVals = splitCellVals;

    package.splitCellGradsX = splitCellGradsX;
    package.splitCellGradsY = splitCellGradsY;

    package.directions = directions;        
    package.oppdirs = oppdirs;
    package.localNodes = localNodes;

    return package;
}