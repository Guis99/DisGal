#include "..\..\Dependencies\QTM\include\QTM.hpp"
#include "..\..\Dependencies\Eigen\Core"
#include "..\..\Dependencies\Eigen\Sparse"
#include "..\..\Dependencies\unsupported\Eigen\KroneckerProduct"
#include "..\..\Dependencies\Utils\Utils.hpp"
// #include "..\..\Meshing\Meshing.hpp"

#ifndef MatrixAssembly_diff
#define MatrixAssembly_diff

typedef Eigen::SparseMatrix<double> SpD;
// Dynamically-sized matrix of doubles
typedef Eigen::MatrixXd DD;
// Dynamically-sized vector of doubles
typedef Eigen::VectorXd DvD;

DD GenerateQuadWeights(std::vector<double>& gpX, std::vector<double> &gpY, int numXNodes, int numYNodes, int numElemNodes);
SpD StiffnessMatrix(QTM::QuadTreeMesh& mesh, double k);
SpD PenaltyMatrix(QTM::QuadTreeMesh& mesh, double k, double alpha);
SpD FluxMatrix(QTM::QuadTreeMesh& mesh, double k);
DvD AssembleFVec(QTM::QuadTreeMesh& mesh, double f, std::string evalStr);;
std::vector<double> ComputeResiduals(QTM::QuadTreeMesh&, DvD& solution, SpD& source);
std::vector<std::shared_ptr<QTM::Cell>> TestResiduals(DvD& solution, QTM::QuadTreeMesh& mesh, double residualLimit);
DvD EvalDirichletBoundaryCond(QTM::QuadTreeMesh& inputMesh, std::vector<std::vector<int>>& boundaryNodes, std::vector<int>& allBoundaryNodes, std::vector<std::string>& strs);
void GetExtensionMatrices(QTM::QuadTreeMesh& inputMesh,
                                        std::vector<int>& boundaryNodes, 
                                        std::vector<int>& freeNodes,
                                        SpD& nullSpace,
                                        SpD& columnSpace);
DvD ComputeSolutionStationaryLinear(SpD& KMatrix, DvD& FMatrix, SpD& columnSpace, SpD& nullSpace, DvD& dirichletBoundaryVals);
DD PoissonSolve(QTM::QuadTreeMesh& inputMesh,
                double c,
                double k,
                std::string source,
                std::vector<std::string> bcs,
                double penaltyParam);
#endif