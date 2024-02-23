#include "include\MatrixAssembly.hpp"

using namespace QTM;

int main() {
    // int nx = 3; int ny = 3;
    // QuadTreeMesh mesh(2, nx, ny, 3, 2);

    // int deg = mesh.deg;
    // int numNodes = deg+1;
    // int numElemNodes = numNodes * numNodes;
    // int nElements = mesh.numLeaves;
    // int nNodes = nElements * numNodes;

    // std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
    // std::vector<double> integX = Utils::integrateLagrange(gaussPoints);

    // integX = {1,2,3};

    // // Convert 1D quad weights to diag matrix
    // Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    // DD quadWeights1D = (DD)integXMat.asDiagonal();

    // DD B; B.setIdentity(numElemNodes, numElemNodes);
    // std::cout<<"printing B:"<<std::endl;
    // std::cout<<B<<std::endl;

    // std::vector<QTM::Direction> directions = {QTM::Direction::N, QTM::Direction::E, 
    //                                             QTM::Direction::S, QTM::Direction::W};

    // std::vector<std::vector<int>> localNodes = {mesh.GetLocalBoundaryNodes(QTM::Direction::N),
    //                                             mesh.GetLocalBoundaryNodes(QTM::Direction::E),
    //                                             mesh.GetLocalBoundaryNodes(QTM::Direction::S),
    //                                             mesh.GetLocalBoundaryNodes(QTM::Direction::W)};

    // Direction dir = Direction::N;

    // std::vector<int> neighborNodes;
    // std::vector<int> elemNodes;
    // std::vector<std::shared_ptr<QTM::Cell>> neighbors;

    // elemNodes = mesh.GetGlobalBoundaryNodes(dir, 4);
    // neighbors = mesh.GetCellNeighbors(dir, 4);
    // auto neighbor = neighbors[0];
    // neighborNodes = mesh.GetGlobalBoundaryNodes((QTM::Direction)((dir+2)%4), neighbor->CID);

    // DD topJump = B(Eigen::all, localNodes[dir]);
    // DD bottomJump = -B(Eigen::all, localNodes[(QTM::Direction)((dir+2)%4)]);
    // std::cout<<"printing topJump:"<<std::endl;
    // std::cout<<topJump<<std::endl;
    // std::cout<<"printing bottomJump:"<<std::endl;
    // std::cout<<bottomJump<<std::endl;
    // DD jumpMatrix(2*numElemNodes, elemNodes.size());
    // jumpMatrix << topJump, bottomJump;

    // std::cout<<"printing jumpMatrix:"<<std::endl;
    // std::cout<<jumpMatrix<<std::endl;

    // auto jumpMatrixT = jumpMatrix.transpose();

    // DD localElemMat(2*numElemNodes, 2*numElemNodes);
    // localElemMat = (DD)(jumpMatrix * quadWeights1D * jumpMatrixT);

    // std::cout<<"printing quadWeights1D:"<<std::endl;
    // std::cout<<quadWeights1D<<std::endl;

    // std::cout<<"printing localElemMat:"<<std::endl;
    // std::cout<<localElemMat<<std::endl;

    // // Generate derivative matrix
    // std::vector<double> AInitializer; 
    // AInitializer.reserve(numElemNodes * numElemNodes);

    // double *AInitIdx = AInitializer.data(); 

    // // Generate derivatives for each basis function, copy to full array
    // for (int k=0; k<numNodes; k++) { 
    //     std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPoints, gaussPoints);
    //     std::copy(xPartials.begin(), xPartials.end(), AInitIdx);
    //     AInitIdx += numNodes;
    // }

    // AInitializer = {1,2,3,4,5,6,7,8,9};
    
    // // map derivative values to matrix
    // Eigen::Map<DD> A(AInitializer.data(), numNodes, numNodes); 

    // DD Bs; Bs.setIdentity(numNodes, numNodes);

    // // Get element-wise matrix intermediates
    // DD combinedX(numElemNodes, numElemNodes);
    // combinedX << Eigen::kroneckerProduct(Bs, A);
    // DD combinedY(numElemNodes, numElemNodes);
    // combinedY << Eigen::kroneckerProduct(A, Bs);

    // std::cout<<"printing combinedX:"<<std::endl;
    // std::cout<<combinedX<<std::endl;

    // std::cout<<"printing combinedY:"<<std::endl;
    // std::cout<<combinedY<<std::endl;

    // DD topGradX = combinedX(localNodes[dir], Eigen::all);
    // DD bottomGradX = -combinedX(localNodes[(QTM::Direction)((dir+2)%4)], Eigen::all);

    // std::cout<<"printing topGradX:"<<std::endl;
    // std::cout<<topGradX<<std::endl;

    // std::cout<<"printing bottomGradX:"<<std::endl;
    // std::cout<<bottomGradX<<std::endl;

    // DD topGradY = combinedY(localNodes[dir], Eigen::all);
    // DD bottomGradY = -combinedY(localNodes[(QTM::Direction)((dir+2)%4)], Eigen::all);

    // std::cout<<"printing topGradY:"<<std::endl;
    // std::cout<<topGradY<<std::endl;

    // std::cout<<"printing bottomGradY:"<<std::endl;
    // std::cout<<bottomGradY<<std::endl;

    // DD fluxMatrixX(elemNodes.size(), 2*numElemNodes);
    // DD fluxMatrixY(elemNodes.size(), 2*numElemNodes);

    // fluxMatrixX << topGradX, bottomGradX;
    // fluxMatrixY << topGradY, bottomGradY; 

    // std::cout<<"printing fluxMatrixX:"<<std::endl;
    // std::cout<<fluxMatrixX<<std::endl;

    // std::cout<<"printing fluxMatrixY:"<<std::endl;
    // std::cout<<fluxMatrixY<<std::endl;


    // DvD topVec = DvD::Zero(numNodes); topVec(0) = 1;
    // DvD bottomVec = DvD::Zero(numNodes); bottomVec(deg) = 1;

    // for (int i=0; i<numNodes-1; i++) {
    //     topVec(i+1) = 0;
    //     bottomVec(i) = 0;
    // }
    // std::cout<<topVec<<std::endl;
    // std::cout<<"print"<<bottomVec<<std::endl;

    std::array<DD,2> mv;
    DD mat(3,3);
    mat << 1,2,3,4,5,6,7,8,9; mv[0] = mat;
    std::cout<<mv[0]<<std::endl;
    mat<<2,4,6,8,10,12,14,16,18; mv[1] = mat;
    std::cout<<mv[0]<<std::endl<<mv[1]<<std::endl; 
    mat<<9,8,7,6,5,4,3,2,1;
    std::cout<<mv[0]<<std::endl<<mv[1]<<std::endl;
    mat.transposeInPlace();

    auto g=mat.data();
    for (int i=0; i<9; i++) {
        std::cout<<*g<<", ";
        g++;
    }

    std::vector<double> v = {1,2,3,4,5};
    Eigen::Map<DD> topVals(v.data(),1,3);
    Eigen::Map<DD> bottomVals(v.data()+2,1,3);

    std::cout<<bottomVals<<std::endl<<topVals<<std::endl;

}