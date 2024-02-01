#include "include\MatrixAssembly.hpp"

using namespace QTM;

int main() {
    int nx = 3; int ny = 3;
    QuadTreeMesh mesh(2, nx, ny, 3, 2);

    int deg = mesh.deg;
    int numNodes = deg+1;
    int numElemNodes = numNodes * numNodes;
    int nElements = mesh.numLeaves;
    int nNodes = nElements * numNodes;

    std::vector<double> gaussPoints = Utils::genGaussPoints(deg);
    std::vector<double> integX = Utils::integrateLagrange(gaussPoints);

    integX = {1,2,3};

    // Convert 1D quad weights to diag matrix
    Eigen::Map<DvD> integXMat(integX.data(),numNodes);
    DD quadWeights1D = (DD)integXMat.asDiagonal();

    DD B; B.setIdentity(numElemNodes, numElemNodes);
    std::cout<<"printing B:"<<std::endl;
    std::cout<<B<<std::endl;

    std::vector<QTM::Direction> directions = {QTM::Direction::N, QTM::Direction::E, 
                                                QTM::Direction::S, QTM::Direction::W};

    std::vector<std::vector<int>> localNodes = {mesh.GetLocalBoundaryNodes(QTM::Direction::N),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::E),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::S),
                                                mesh.GetLocalBoundaryNodes(QTM::Direction::W)};

    Direction dir = Direction::N;

    std::vector<int> neighborNodes;
    std::vector<int> elemNodes;
    std::vector<std::shared_ptr<QTM::Cell>> neighbors;

    elemNodes = mesh.GetBoundaryNodes(dir, 4);
    neighbors = mesh.GetCellNeighbors(dir, 4);
    auto neighbor = neighbors[0];
    neighborNodes = mesh.GetBoundaryNodes((QTM::Direction)((dir+2)%4), neighbor->CID);

    DD topJump = B(Eigen::all, localNodes[dir]);
    DD bottomJump = -B(Eigen::all, localNodes[(QTM::Direction)((dir+2)%4)]);
    std::cout<<"printing topJump:"<<std::endl;
    std::cout<<topJump<<std::endl;
    std::cout<<"printing bottomJump:"<<std::endl;
    std::cout<<bottomJump<<std::endl;
    DD jumpMatrix(2*numElemNodes, elemNodes.size());
    jumpMatrix << topJump, bottomJump;

    std::cout<<"printing jumpMatrix:"<<std::endl;
    std::cout<<jumpMatrix<<std::endl;

    auto jumpMatrixT = jumpMatrix.transpose();

    DD localElemMat(elemNodes.size() + neighborNodes.size(), elemNodes.size() + neighborNodes.size());
    localElemMat = (DD)(jumpMatrix * quadWeights1D * jumpMatrixT);

    std::cout<<"printing quadWeights1D:"<<std::endl;
    std::cout<<quadWeights1D<<std::endl;

    std::cout<<"printing localElemMat:"<<std::endl;
    std::cout<<localElemMat<<std::endl;
}