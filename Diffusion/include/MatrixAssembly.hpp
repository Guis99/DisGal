#include "../../Dependencies/QTM/include/QTM.hpp"
#include "../../Dependencies/Eigen/Core"
#include "../../Dependencies/Eigen/Sparse"
#include "../../Dependencies/unsupported/Eigen/KroneckerProduct"
#include "../../Dependencies/Utils/Utils.hpp"
// #include "..\..\Meshing\Meshing.hpp"

#ifndef MatrixAssembly_diff
#define MatrixAssembly_diff

typedef Eigen::SparseMatrix<double> SpD;
// Dynamically-sized matrix of doubles
typedef Eigen::MatrixXd DD;
// Dynamically-sized vector of doubles
typedef Eigen::VectorXd DvD;

namespace PMA { // PMA = Poisson Matrix Assembly

struct quadUtils {
    // Quadrature weights
    DD quadWeights1D;
    DD weightMat;

    // full cell nodal vals
    DD nodalVals;

    // full cell nodal grads on all points
    DD combinedX;
    DD combinedY;

    // split cell nodal vals
    std::array<DD,4> splitCellVals;

    // split cell nodal grads
    std::array<DD,4> splitCellGradsX;
    std::array<DD,4> splitCellGradsY;

    // neighbor-finding
    std::vector<QTM::Direction> directions;        
    std::vector<QTM::Direction> oppdirs;
    std::vector<std::vector<int>> localNodes;
};

quadUtils GenerateAssemblyPackage(QTM::QuadTreeMesh& mesh);

DD GenerateQuadWeights(std::vector<double>& gpX, std::vector<double> &gpY, int numXNodes, int numYNodes, int numElemNodes);

SpD StiffnessMatrix(QTM::QuadTreeMesh& mesh, double k);
SpD PenaltyMatrix(QTM::QuadTreeMesh& mesh, double k, double alpha, quadUtils& package);
SpD FluxMatrix(QTM::QuadTreeMesh& mesh, double k, quadUtils& package);
SpD BoundaryMatrix(QTM::QuadTreeMesh& mesh, double k, 
                std::vector<bool> isDirichletBC,
                std::vector<std::string> dbcs,
                double alpha, 
                quadUtils& package);

DvD IntegrateDirichlet(QTM::QuadTreeMesh& mesh,
                        std::vector<bool> isDirichletBC,
                        std::vector<std::string> dbcs,
                        double alpha, 
                        quadUtils& package);

DvD IntegrateNeumann(QTM::QuadTreeMesh& mesh,
                        std::vector<bool> isNeumannBC,
                        std::vector<std::string> nbcs,
                        quadUtils& package);                        

DvD AssembleFVec(QTM::QuadTreeMesh& mesh, double f, std::string evalStr);;
std::vector<double> ComputeResiduals(QTM::QuadTreeMesh&, DvD& solution, SpD& source);
std::vector<std::shared_ptr<QTM::Cell>> TestResiduals(DvD& solution, QTM::QuadTreeMesh& mesh, double residualLimit);
DvD EvalSymbolicBoundaryCond(QTM::QuadTreeMesh& inputMesh, std::vector<std::vector<int>>& boundaryNodes, std::vector<int>& allBoundaryNodes, std::vector<std::string>& strs);
void GetExtensionMatrices(QTM::QuadTreeMesh& inputMesh,
                                        std::vector<int>& boundaryNodes, 
                                        std::vector<int>& freeNodes,
                                        SpD& nullSpace,
                                        SpD& columnSpace);

DvD ComputeSolutionStationaryLinear(SpD& KMatrix, DvD& FMatrix, SpD& columnSpace, SpD& nullSpace, DvD& dirichletBoundaryVals);
DvD ComputeSolutionStationaryLinearNoElim(SpD& A, DvD& b);
void FindRank(SpD& mat);
void FindConditionNumber(SpD& mat);

DvD PoissonSolve(QTM::QuadTreeMesh& inputMesh,
                double c,
                double k,
                std::string source,
                std::vector<std::string> bcs,
                double penaltyParam);
DvD dgPoissonSolve(QTM::QuadTreeMesh& inputMesh,
                double k,
                std::string source,
                std::vector<bool> isDirichletBC,
                std::vector<bool> isNeumannBC,
                std::vector<std::string> dbcs,
                std::vector<std::string> nbcs,
                double penaltyParam,
                double dirichletPenalty);
}
#endif