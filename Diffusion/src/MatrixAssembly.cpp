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
    std::cout<<"d1"<<std::endl;
    // map derivative values to matrix
    Eigen::Map<DD> A(AInitializer.data(), numNodes, numNodes); 
std::cout<<"d2"<<std::endl;
    // Generate quadrature weight matrices
    DD weightMat = GenerateQuadWeights(gaussPoints, gaussPoints, numNodes, numNodes, numElemNodes);
std::cout<<"d3"<<std::endl;
    // Generate mass matrices
    DD B; B.setIdentity(numNodes, numNodes);
std::cout<<"d4"<<std::endl;
    DD coeffMat(numElemNodes, numElemNodes);
    coeffMat.setIdentity(); coeffMat *= k;
    
    // Get element-wise matrix intermediates
    DD combinedX(numElemNodes, numElemNodes);
    combinedX << Eigen::kroneckerProduct(B, A);
    DD combinedY(numElemNodes, numElemNodes);
    combinedY << Eigen::kroneckerProduct(A, B);
std::cout<<"d5"<<std::endl;
    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes * numElemNodes);

    auto combineXT = (DD)combinedX.transpose();
    auto combineYT = (DD)combinedY.transpose();
    std::cout<<"d6"<<std::endl;

    // Integrate over all elements
    DD localElemMat(numElemNodes, numElemNodes);
    auto leaves = mesh.GetAllCells();
    std::cout<<"d7"<<std::endl;
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
    std::cout<<"d10"<<std::endl;
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    std::cout<<"d11"<<std::endl;
    return mat;
}

SpD PenaltyMatrix(QTM::QuadTreeMesh& mesh, double k, double alpha) {
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();

    // basis func to node mapping
    DD B; B.setIdentity(numElemNodes, numElemNodes);
    DD Bs; Bs.setIdentity(2*numNodes-1, 2*numNodes-1);
    std::cout<<"db1"<<std::endl;

    // get basis func vals for split cell quad
    std::vector<double> BhInitializer; 
    BhInitializer.reserve(numNodes * numNodes);
    std::cout<<"db1"<<std::endl;

    // get unit vecs
    DvD topVec = DvD::Zero(numNodes); topVec(0) = 1;
    DvD bottomVec = DvD::Zero(numNodes); bottomVec(deg) = 1;
    std::cout<<"db1"<<std::endl;

    // Generate values for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<double> xVals = Utils::evalLagrangeInterp(k, mesh.halfGaussPoints, mesh.gaussPoints);
        BhInitializer.insert(BhInitializer.end(), xVals.begin(), xVals.end());
    }
    std::cout<<"db1"<<std::endl;
    // map basis values to matrix
    Eigen::Map<DD> BhT(BhInitializer.data(), 2*numNodes-1, numNodes); // eval points are traversed first
    DD Bh = (DD)(BhT.transpose());
    std::cout<<"db1"<<std::endl;

    std::array<DD,4> splitCellVals;
    std::cout<<"db1"<<std::endl;

    DD splitCellPlaceholder(numElemNodes, mesh.halfGaussPoints.size());
    splitCellPlaceholder << Eigen::kroneckerProduct(bottomVec, Bh); splitCellVals[0] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(Bh, bottomVec); splitCellVals[1] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(topVec, Bh); splitCellVals[2] = splitCellPlaceholder;
    splitCellPlaceholder << Eigen::kroneckerProduct(Bh, topVec); splitCellVals[3] = splitCellPlaceholder;
    std::cout<<"db1"<<std::endl;
    
    // get quad weights in 1D
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
    std::vector<double> integX = Utils::integrateLagrange(gaussPoints);
    std::cout<<"db1"<<std::endl;

    // Convert 1D quad weights to diag matrix
    Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    DD quadWeights1D = (DD)integXMat.asDiagonal();
    std::cout<<"db12"<<std::endl;

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
            QTM::Direction oppdir = (QTM::Direction)((dir+2)%4);
            cellSign = edgeSign[dir]; neighborSign = -cellSign;
            neighbors = mesh.GetCellNeighbors(dir, elm->CID);
            elemNodes = mesh.GetGlobalElemNodes(dir, elm->CID);
            for (auto neighbor : neighbors) { 
                if ((neighbor!=nullptr) && (neighbor->CID > elm->CID || elm->level > neighbor->level)) { // wrong side of boundary
                    continue;
                }
                if (neighbor == nullptr) { // case curr boundary is exterior
                    neighborNodes = {};
                    neighbor = elm;
                } else { // case appropriate neighbor exists
                    neighborNodes = mesh.GetGlobalElemNodes(oppdir, neighbor->CID);
                }
                // calculate penalty param
                a = 2*alpha/std::min(elm->width, neighbor->width);
                // calculate jump matrix
                DD topJump = cellSign*B(Eigen::all, localNodes[dir]);
                DD bottomJump = neighborSign*B(Eigen::all, localNodes[oppdir]);
                DD jumpMatrix(2*numElemNodes, numNodes);
                jumpMatrix << topJump, bottomJump;
                boundaryNodes.reserve(boundaryNodes.size() + neighborNodes.size());
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
    AInitializer.reserve(numElemNodes * numElemNodes);

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numNodes; k++) { 
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPoints, gaussPoints);
        AInitializer.insert(AInitializer.end(), xPartials.begin(), xPartials.end());
    }
    
    // map derivative values to matrix
    Eigen::Map<DD> A(AInitializer.data(), numNodes, numNodes); 
    A = .5*A;

    // Get element-wise matrix intermediates
    DD combinedX(numElemNodes, numElemNodes);
    combinedX << Eigen::kroneckerProduct(Bs, A);
    DD combinedY(numElemNodes, numElemNodes);
    combinedY << Eigen::kroneckerProduct(A, Bs);

    // Signs for element edges
    int edgeSign[4] = {1,1,-1-1};
    int cellSign; int neighborSign;

    double normalX[4] = {0,1,0,-1};
    double normalY[4] = {1,0,-1,0};

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
            elemNodes = mesh.GetGlobalElemNodes(dir, elm->CID);
            QTM::Direction oppdir = (QTM::Direction)((dir+2)%4);
            for (auto neighbor : neighbors) { 
                if ((neighbor!=nullptr) && (neighbor->CID > elm->CID || elm->level > neighbor->level)) { // wrong side of boundary
                    continue;
                }
                if (neighbor == nullptr) { // case curr boundary is exterior
                    neighborNodes = {};
                    neighbor = elm;
                } else { // case appropriate neighbor exists
                    neighborNodes = mesh.GetGlobalElemNodes(oppdir, neighbor->CID);
                }
                // calculate jump matrix
                DD topJump = cellSign*B(Eigen::all, localNodes[dir]);
                DD bottomJump = neighborSign*B(Eigen::all, localNodes[oppdir]);
                DD jumpMatrix(2*numElemNodes, numNodes);
                jumpMatrix << topJump, bottomJump;

                // calculate flux matrix
                DD topGradX = normalX[dir] * combinedX(localNodes[dir], Eigen::all);
                DD bottomGradX = normalX[oppdir] * combinedX(localNodes[oppdir], Eigen::all );

                DD topGradY = normalY[dir] * combinedY(localNodes[dir], Eigen::all);
                DD bottomGradY = normalY[oppdir] * combinedY(localNodes[oppdir], Eigen::all);

                DD fluxMatrixX(numNodes, 2*numElemNodes);
                DD fluxMatrixY(numNodes, 2*numElemNodes);

                // place partial derivatives in combined mat
                fluxMatrixX << topGradX, bottomGradX;
                fluxMatrixY << topGradY, bottomGradY;

                boundaryNodes.reserve(elemNodes.size() + neighborNodes.size());
                boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                boundaryNodes.insert(boundaryNodes.end(), neighborNodes.begin(), neighborNodes.end());

                // assemble local matrix 
                DD localElemMat = (DD)(jumpMatrix * elm->width*quadWeights1D * fluxMatrixX + 
                                    jumpMatrix * elm->width*quadWeights1D * fluxMatrixY);

                if (neighbor->CID == elm->CID) {
                    localElemMat = (DD)(2*localElemMat);
                }

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
    std::cout<<"db1"<<std::endl;

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
        double Lx = elm->width; // Jacobian factors
        // calculate local matrix
        std::vector<double> collectSourceVals; collectSourceVals.reserve(numElemNodes);
        auto nodes = elm->nodes;
        for (int i=nodes[0]; i<=nodes[1]; i++) {
            collectSourceVals.push_back(fEval[i]);
        }
        Eigen::Map<DvD> sourceVector(collectSourceVals.data(), numElemNodes, 1);

        localElemMat = weightMat*sourceVector*Lx*Lx/4;
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

double ComputeResidual() {

}

std::vector<std::shared_ptr<QTM::Cell>> TestResiduals(DvD& solution, QTM::QuadTreeMesh& mesh, double residualLimit) {
    std::vector<std::shared_ptr<QTM::Cell>> out; out.reserve(mesh.leaves.size());
    for (auto leaf : mesh.leaves) {
        double res = ComputeResidual();
        if (leaf->level < 5 && ComputeResidual() > residualLimit) {
             out.push_back(leaf);
        }
    }

    return out;
}

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
    std::cout<<numBoundaryNodes<<", "<<BcValsSorted.size()<<std::endl;

    Eigen::Map<DvD> boundaryNodeValuesVec(BcValsSorted.data(), numBoundaryNodes, 1);
    return (DvD)boundaryNodeValuesVec;
}

void GetExtensionMatrices(QTM::QuadTreeMesh& inputMesh,
                                        std::vector<int>& boundaryNodes, 
                                        std::vector<int>& freeNodes,
                                        SpD& nullSpace,
                                        SpD& columnSpace) {
    int nNodes = inputMesh.nNodes();

    std::cout<<"db1"<<std::endl;

    std::sort(boundaryNodes.begin(), boundaryNodes.end());

    std::cout<<"db1"<<std::endl;

    std::vector<Eigen::Triplet<double>> tripletListNS;
    std::vector<Eigen::Triplet<double>> tripletListCS;
    std::cout<<"db1"<<std::endl;
    tripletListNS.reserve(boundaryNodes.size());
    tripletListCS.reserve(freeNodes.size());
    std::cout<<"db1"<<std::endl;

    for (int i=0; i<boundaryNodes.size(); i++) {
        tripletListNS.emplace_back(boundaryNodes[i], i, 1.0);
    }
    std::cout<<"db1"<<std::endl;

    for (int i=0; i<freeNodes.size(); i++) {
        tripletListCS.emplace_back(freeNodes[i], i, 1.0);
    }
    std::cout<<"db1"<<std::endl;

    std::cout<<boundaryNodes.size()<<", "<<nullSpace.rows()<<", "<<nullSpace.cols()<<std::endl;
    std::cout<<freeNodes.size()<<", "<<columnSpace.rows()<<", "<<columnSpace.cols()<<std::endl;
    
    nullSpace.setFromTriplets(tripletListNS.begin(), tripletListNS.end());
    columnSpace.setFromTriplets(tripletListCS.begin(), tripletListCS.end());
}

DvD ComputeSolutionStationary(SpD& StiffnessMatrix, SpD& PenaltyMatrix, SpD& FluxMatrix, SpD& fVec, SpD& columnSpace, SpD& nullSpace, DvD& boundaryVals) {
    SpD FluxMatrixT = (SpD)(FluxMatrix.transpose());
    SpD CombinedMatrix = StiffnessMatrix - FluxMatrix + FluxMatrixT + PenaltyMatrix;
    // Eliminate rows and columns corr. to boundary nodes
    SpD columnSpaceT = (SpD)(columnSpace.transpose());
    SpD A11 = columnSpaceT * CombinedMatrix * columnSpace;
    // Eliminate boundary rows and free columns
    SpD A12 = columnSpaceT * CombinedMatrix * nullSpace;
    // Eliminate boundary rows
    SpD F11 = columnSpaceT * fVec;

    Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolver;    
    LuSolver.analyzePattern(A11);
    LuSolver.factorize(A11);

    DvD x = LuSolver.solve(F11 - A12 * boundaryVals);
    x = columnSpace * x + nullSpace * boundaryVals;
    std::cout<<x.rows()<<std::endl;
    return x;
}

DD PoissonSolve(QTM::QuadTreeMesh& inputMesh,
                double c,
                double k,
                std::string source,
                std::vector<std::string> bcs,
                int penaltyParam) {
    
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