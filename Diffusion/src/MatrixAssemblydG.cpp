#include "..\include\MatrixAssembly.hpp"

SpD PMA::BoundaryMatrix(QTM::QuadTreeMesh& mesh, double k, 
                std::vector<bool> isDirichletBC,
                std::vector<std::string> dbcs,
                double alpha, 
                PMA::quadUtils& package) {
    // Integral on Dirichlet boundaries of: (a/h)*[u][v] - {k * grad(u) dot n}[v] - {k * grad(v) dot n}[u]
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = mesh.numElemNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();

    DD quadWeights1D = package.quadWeights1D;
    DD nodalVals = package.nodalVals;

    DD combinedX = package.combinedX;
    DD combinedY = package.combinedY;

    std::cout<<"--------package in boundary------------"<<std::endl;
    std::cout<<combinedX.rows()<<", "<<combinedX.cols()<<std::endl;
    std::cout<<combinedY.rows()<<", "<<combinedY.cols()<<std::endl;
    std::cout<<"--------package in boundary------------"<<std::endl;
    
    std::array<DD, 4> splitCellVals = package.splitCellVals;
    std::array<DD, 4> splitCellGradsX = package.splitCellGradsX;
    std::array<DD, 4> splitCellGradsY = package.splitCellGradsY;

    std::vector<QTM::Direction> directions = package.directions;        
    std::vector<QTM::Direction> oppdirs = package.oppdirs;
    std::vector<std::vector<int>> localNodes = package.localNodes;
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
    std::vector<int> elemNodes;
    std::vector<int> elemLocals;
    std::vector<std::shared_ptr<QTM::Cell>> neighbors;
    std::vector<Eigen::Triplet<double>> tripletList; 
    auto leaves = mesh.GetAllCells();

    std::vector<std::shared_ptr<QTM::Cell>> dirichletCells; 
    std::vector<QTM::Direction> boundaryDirs = {QTM::Direction::N,
                                   QTM::Direction::E,
                                   QTM::Direction::S,
                                   QTM::Direction::W};

    for (int i=0; i<isDirichletBC.size(); i++) {
        if (isDirichletBC[i]) {
            QTM::Direction dir = boundaryDirs[i];
            dirichletCells = mesh.boundaryCells[i];

            for (const auto& cell : dirichletCells) {
                elemNodes = mesh.GetGlobalElemNodes(cell->CID);
                elemLocals = mesh.GetTrimmedLocalNodes(cell->CID, elemNodes);

                auto nodes = cell->nodes;
                double jac = cell->width; // Jacobian factors
                double penalty = alpha / jac / 2;

                // jump matrix setup
                DD topJump;

                // flux matrix setup
                DD topGradX;
                DD topGradY;

                topJump = nodalVals(elemLocals, localNodes[dir]);
                topGradX = normalX[dir] * combinedX(localNodes[dir], elemLocals);
                topGradY = normalY[dir] * combinedY(localNodes[dir], elemLocals);

                // std::cout<<"tj: "<<topJump.rows()<<", "<<topJump.cols()<<std::endl;
                // std::cout<<"tgx: "<<topGradX.rows()<<", "<<topGradX.cols()<<std::endl;
                // std::cout<<"tgy: "<<topGradX.rows()<<", "<<topGradY.cols()<<std::endl;

                // calculate jump matrix
                DD jumpMatrix(elemNodes.size(), numNodes);
                // std::cout<<"jm: "<<jumpMatrix.rows()<<", "<<jumpMatrix.cols()<<std::endl;
                jumpMatrix << topJump;
          
                // Utils::printVec(elemNodes);
                // Utils::printVec(localNodes[dir]);
                // Utils::printVec(elemLocals);

                DD fluxMatrixX(numNodes, elemNodes.size());
                DD fluxMatrixY(numNodes, elemNodes.size());

                // std::cout<<"fmX: "<<fluxMatrixX.rows()<<", "<<fluxMatrixX.cols()<<std::endl;
                // std::cout<<"fmY: "<<fluxMatrixY.rows()<<", "<<fluxMatrixY.cols()<<std::endl;

                // place partial derivatives in combined mat
                fluxMatrixX << topGradX;
                fluxMatrixY << topGradY;

                // assemble local matrix 
                DD jumpMatrixT = (DD)(jumpMatrix.transpose());
                DD fluxMatrixXT = (DD)(fluxMatrixX.transpose());
                DD fluxMatrixYT = (DD)(fluxMatrixY.transpose());

                // std::cout<<"miscellaneous matrix assembly check"<<std::endl;
                // std::cout<<jumpMatrix.rows()<<", "<<jumpMatrix.cols()<<std::endl;
                // std::cout<<fluxMatrixX.rows()<<", "<<fluxMatrixX.cols()<<std::endl;
                // std::cout<<fluxMatrixY.rows()<<", "<<fluxMatrixY.cols()<<std::endl;
                // std::cout<<quadWeights1D.rows()<<std::endl;

                DD localElemMat = penalty*jumpMatrix * jac*quadWeights1D * jumpMatrixT;
                DD localElemMatGrad = (DD)(jumpMatrix * quadWeights1D * fluxMatrixX + 
                                        jumpMatrix * quadWeights1D * fluxMatrixY); // TODO: check that matrix dimensions are good
                DD localElemMatGradT = (DD)(fluxMatrixXT * quadWeights1D * jumpMatrixT + 
                                        fluxMatrixYT * quadWeights1D * jumpMatrixT);

                localElemMat -= localElemMatGrad;
                localElemMat -= localElemMatGradT;

                boundaryNodes.reserve(elemNodes.size());
                boundaryNodes.insert(boundaryNodes.end(), elemNodes.begin(), elemNodes.end());
                
                for (int j=0; j<boundaryNodes.size(); j++) {
                    for (int i=0; i<boundaryNodes.size(); i++) {
                        tripletList.emplace_back(boundaryNodes[i],boundaryNodes[j],localElemMat(i,j));
                    }
                }
                boundaryNodes.clear();
            }
        } else {
            continue;
        }
    }
        
    // Declare and construct sparse matrix from triplets
    SpD mat(nNodes,nNodes);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    mat.makeCompressed();
    return mat;
}

DvD PMA::IntegrateDirichlet(QTM::QuadTreeMesh& mesh,
                        std::vector<bool> isDirichletBC,
                        std::vector<std::string> dbcs,
                        double alpha, 
                        PMA::quadUtils& package) {
    // do integral (v - n dot grad v) * u_D   over dirichlet boundary
    std::cout<<"dbc bool:"<<std::endl;
    for (bool elm : isDirichletBC) {
        std::cout<<elm<<std::endl;
    }
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = mesh.numElemNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);

    std::vector<int> allBoundaryNodes;

    // unpack quad utils
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

    std::vector<std::shared_ptr<QTM::Cell>> dirichletCells; 
    std::vector<QTM::Direction> boundaryDirs = {QTM::Direction::N,
                                   QTM::Direction::E,
                                   QTM::Direction::S,
                                   QTM::Direction::W};

    std::vector<int> elemNodes;
    std::vector<int> elemLocals;

    double normalX[4] = {0,1,0,-1};
    double normalY[4] = {1,0,-1,0};
    int dbcIdx = 0;

    DvD out(nNodes);

    for (int i=0; i<isDirichletBC.size(); i++) {
        if (isDirichletBC[i]) {
            QTM::Direction dir = boundaryDirs[i];
            dirichletCells = mesh.boundaryCells[i];
            std::vector<int> dirichletCellIDs; 
            dirichletCellIDs.reserve(dirichletCells.size());
            std::cout<<"dbc #"<<dbcIdx<<": "<<std::endl;
            for (const auto& cell : dirichletCells) {
                std::cout<<cell->CID<<std::endl;
                dirichletCellIDs.push_back(cell->CID);
            }
            std::cout<<"---------------"<<std::endl;

            std::vector<int> localBoundaryNodes = localNodes[dir];

            std::vector<int> boundaryCellPos = Utils::GetBoundaryNodes(localBoundaryNodes, dirichletCellIDs, numNodes);

            auto nodesPos = mesh.GetNodePos(boundaryCellPos); // cell IDs won't work, rather need node idxs
            auto startpoint = nodesPos.data(); auto allocSize = nodesPos.size();
            auto fEval = Utils::EvalSymbolicBC(startpoint, allocSize, dbcs[dbcIdx++]);
            Utils::printVec(fEval);
            

            int cellIterIdx = 0;
            auto dirichletIdx = fEval.begin();
            for (const auto& cell : dirichletCells) {
                elemNodes = mesh.GetGlobalElemNodes(cell->CID);
                elemLocals = mesh.GetTrimmedLocalNodes(cell->CID, elemNodes);

                auto nodes = cell->nodes;
                double jac = cell->width; // Jacobian factors
                double penalty = alpha / cell->width / 2;

                std::vector<double> collectSourceVals; collectSourceVals.reserve(numNodes);
                collectSourceVals.insert(collectSourceVals.begin(), dirichletIdx, dirichletIdx+numNodes);
                dirichletIdx += numNodes;

                Eigen::Map<DvD> sourceVector(collectSourceVals.data(), numNodes, 1);
                // jump matrix setup
                DD topJump;

                // flux matrix setup
                DD topGradX;
                DD topGradY;

                topJump = nodalVals(elemLocals, localNodes[dir]);
                topGradX = normalX[dir] * combinedX(localNodes[dir], elemLocals);
                topGradY = normalY[dir] * combinedY(localNodes[dir], elemLocals);

                DD topGradXT = topGradX.transpose();
                DD topGradYT = topGradY.transpose();

                std::cout<<"tj: "<<topJump.rows()<<", "<<topJump.cols()<<std::endl;
                std::cout<<"tgx: "<<topGradXT.rows()<<", "<<topGradXT.cols()<<std::endl;
                std::cout<<"tgy: "<<topGradXT.rows()<<", "<<topGradYT.cols()<<std::endl;

                // calculate jump matrix
                DD jumpMatrix(elemNodes.size(), numNodes);

                std::cout<<"jm: "<<jumpMatrix.rows()<<", "<<jumpMatrix.cols()<<std::endl;
                jumpMatrix << topJump;

                DD fluxMatrixX(elemNodes.size(), numNodes);
                DD fluxMatrixY(elemNodes.size(), numNodes);

                // place partial derivatives in combined mat
                fluxMatrixX << topGradXT;
                fluxMatrixY << topGradYT;

                std::cout<<"miscellaneous matrix assembly check"<<std::endl;
                std::cout<<jumpMatrix.rows()<<", "<<jumpMatrix.cols()<<std::endl;
                std::cout<<fluxMatrixX.rows()<<", "<<fluxMatrixX.cols()<<std::endl;
                std::cout<<fluxMatrixY.rows()<<", "<<fluxMatrixY.cols()<<std::endl;
                std::cout<<quadWeights1D.rows()<<", "<<quadWeights1D.cols()<<std::endl;
                std::cout<<sourceVector.rows()<<", "<<sourceVector.cols()<<std::endl;

                // assemble local matrix 
                DvD localElemMat = penalty*topJump * jac*quadWeights1D * sourceVector;
                DvD localElemMatGrad = (DvD)(fluxMatrixX * quadWeights1D * sourceVector + 
                                        fluxMatrixY * quadWeights1D * sourceVector); 

                localElemMat -= localElemMatGrad;
                
                for (int i=0; i<numElemNodes; i++) {
                    out[nodes[0]+i] += localElemMat[i];
                }
            }
        } else {
            continue;
        }
    }
    return out;
} 

DvD PMA::IntegrateNeumann(QTM::QuadTreeMesh& mesh,
                        std::vector<bool> isNeumannBC,
                        std::vector<std::string> nbcs,
                        PMA::quadUtils& package) {
    // do integral v * u_N over neumann boundary
    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = mesh.numElemNodes;
    int nElements = mesh.numLeaves;
    int nNodes = mesh.nNodes();
    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);

    // load package data
    DD quadWeights1D = package.quadWeights1D;
    DD nodalVals = package.nodalVals;

    // basis func to node mapping
    DD B; B.setIdentity(numElemNodes, numElemNodes);
    // Generate quadrature weight matrices
    DD weightMat = GenerateQuadWeights(gaussPoints, gaussPoints, numNodes, numNodes, numElemNodes);

    std::vector<std::shared_ptr<QTM::Cell>> neumannCells;
    std::vector<QTM::Direction> boundaryDirs = {QTM::Direction::N,
                                   QTM::Direction::E,
                                   QTM::Direction::S,
                                   QTM::Direction::W};
    
    // Initialize i,j,v triplet list for sparse matrix
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(nElements * numElemNodes);

    DvD out(nNodes);

    std::vector<std::vector<int>> localNodes = {mesh.GetLocalBoundaryNodes(QTM::Direction::N),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::E),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::S),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::W)};

    double a; // penalty parameters

    std::vector<int> elemNodes;
    std::vector<int> elemLocals;
    
    // Integrate over all elements
    DvD localElemMat(numElemNodes);
    int nbcIdx = 0;
    for (int i=0; i<isNeumannBC.size(); i++) {
        if (isNeumannBC[i]) {
            QTM::Direction dir = boundaryDirs[i];
            neumannCells = mesh.boundaryCells[i];
            std::vector<int> neumannCellIDs; 
            neumannCellIDs.reserve(neumannCells.size());
            for (const auto& cell : neumannCells) {
                neumannCellIDs.push_back(cell->CID);
            }

            auto nodesPos = mesh.GetNodePos(neumannCellIDs);
            auto startpoint = nodesPos.data(); auto allocSize = nodesPos.size();
            auto fEval = Utils::EvalSymbolicBC(startpoint, allocSize, nbcs[nbcIdx++]);
            for (const auto& elm : neumannCells) {
                elemNodes = mesh.GetGlobalElemNodes(elm->CID);
                elemLocals = mesh.GetTrimmedLocalNodes(elm->CID, elemNodes);
                double jac = elm->width; // Jacobian factors
                // calculate local matrix
                auto nodes = elm->nodes;
                std::vector<double> collectSourceVals; collectSourceVals.reserve(numElemNodes);
                collectSourceVals.insert(collectSourceVals.begin(), fEval.begin()+nodes[0], fEval.begin()+nodes[1]+1);

                Eigen::Map<DvD> sourceVector(collectSourceVals.data(), numElemNodes, 1);

                DD topJump = nodalVals(elemLocals, localNodes[dir]);
                localElemMat = topJump * jac*quadWeights1D * sourceVector;
                // Get nodes in element
                
                // Generate i,j,v triplets
                for (int i=0; i<numElemNodes; i++) {
                    // tripletList.emplace_back(nodes[0]+i, 0, localElemMat(i));
                    out[nodes[0]+i] += localElemMat[i];
                }
            }
        } else {
            continue;
        }
    }
    return out;
}

DvD PMA::ComputeSolutionStationaryLinearNoElim(SpD& A, DvD& b) {
    using namespace Eigen;

    SparseLU<SpD,COLAMDOrdering<int>> solver;
    solver.analyzePattern(A);
    PMA::FindRank(A);
    // PMA::FindConditionNumber(A);
    solver.factorize(A);
    DvD x = solver.solve(b); 

    // ConjugateGradient<SpD, Lower|Upper> cg;
    // cg.compute(A11);
    // DvD x = cg.solve(F11 - A12 * boundaryVals);

    return x;
}

DvD PMA::dgPoissonSolve(QTM::QuadTreeMesh& inputMesh,
                double k,
                std::string source,
                std::vector<bool> isDirichletBC,
                std::vector<bool> isNeumannBC,
                std::vector<std::string> dbcs,
                std::vector<std::string> nbcs,
                double penaltyParam,
                double dirichletPenalty) {
    
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
    std::cout<<"Assembling flux matrix"<<std::endl;
    SpD Miscellaneous = PMA::BoundaryMatrix(inputMesh, k, isDirichletBC, dbcs, penaltyParam, package);
    std::cout<<"Assembling overall system matrix"<<std::endl;
    SpD StiffnessMatrix = KMatrix + PMatrix - SMatrix - SMatrixT + Miscellaneous;

    // std::cout<<"All:/n"<<StiffnessMatrix<<std::endl;
    // std::cout<<"Boundary:/n"<<Miscellaneous<<std::endl;

    std::cout<<"Assembling RHS vector"<<std::endl;
    DvD FMatrix = PMA::AssembleFVec(inputMesh, 1.0, source);
    std::cout<<"Assembling boundary condition vectors"<<std::endl;

    DvD dirichletBound = PMA::IntegrateDirichlet(inputMesh, isDirichletBC, dbcs, penaltyParam, package);
    DvD neumannBound = PMA::IntegrateNeumann(inputMesh, isNeumannBC, nbcs, package);

    std::cout<<"BOUND: "<< allBoundaryNodes.size()<<" FREE: "<< freeNodes.size()<<" ALL: "<<inputMesh.nNodes()<<std::endl;
    std::cout<<FMatrix.rows()<<", "<<dirichletBound.rows()<<", "<<neumannBound.rows()<<std::endl;

    DvD b = FMatrix + dirichletBound + neumannBound;

    std::cout<<"Solving system with "<<inputMesh.nNodes()<<" nodes"<<std::endl;
    DvD x = PMA::ComputeSolutionStationaryLinearNoElim(StiffnessMatrix, b);
    std::cout<<"System solved!"<<std::endl;
    return x;
}