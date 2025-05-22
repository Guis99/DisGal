#include "../../Dependencies/QTM/include/QTM.hpp"
#include "../../Dependencies/Eigen/Core"
#include "../../Dependencies/Eigen/Sparse"
#include "../../Dependencies/unsupported/Eigen/KroneckerProduct"
#include "../../Utils/Utils.hpp"
// #include "..\..\Meshing\Meshing.hpp"

#ifndef MatrixAssembly_diff
#define MatrixAssembly_diff

#ifndef VARPRINT
#define VARPRINT
template<typename... Args>
void var_print(Args&&... args) {
    (std::cout << ... << args) << std::endl; // Fold expression to print all arguments
}
#endif

#ifndef DEBUG_PRINT
#ifdef VERBOSE
#define DEBUG_PRINT(...) var_print(__VA_ARGS__)
#else
#define DEBUG_PRINT(...)
#endif
#endif

#ifdef MULTITHREAD
    #define CONFIRM_MT(numThreads) var_print("Multithreading active. Using ", numThreads, " threads")
    #define MT_ACTIVE 1
#else
    #define CONFIRM_MT(numThreads)
    #define MT_ACTIVE 0
#endif

// static_assert(sizeof(Real_b) == 8);

typedef Eigen::SparseMatrix<Real_b> SpD;
// Dynamically-sized matrix of Real_bs
typedef Eigen::Matrix<Real_b, Eigen::Dynamic, Eigen::Dynamic> DD;
// Dynamically-sized vector of Real_bs
typedef Eigen::Matrix<Real_b, Eigen::Dynamic, 1> DvD;

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

DD GenerateQuadWeights(std::vector<Real_b>& gpX, std::vector<Real_b> &gpY, int numXNodes, int numYNodes, int numElemNodes);

SpD StiffnessMatrix(QTM::QuadTreeMesh& mesh, Real_b k);
SpD PenaltyMatrix(QTM::QuadTreeMesh& mesh, Real_b k, Real_b alpha, quadUtils& package);
SpD FluxMatrix(QTM::QuadTreeMesh& mesh, Real_b k, quadUtils& package);
SpD BoundaryMatrix(QTM::QuadTreeMesh& mesh, Real_b k, 
                std::vector<bool> isDirichletBC,
                std::vector<std::string> dbcs,
                Real_b alpha, 
                quadUtils& package);

DvD IntegrateDirichlet(QTM::QuadTreeMesh& mesh,
                        std::vector<bool> isDirichletBC,
                        std::vector<std::string> dbcs,
                        Real_b alpha, 
                        quadUtils& package);

DvD IntegrateNeumann(QTM::QuadTreeMesh& mesh,
                        std::vector<bool> isNeumannBC,
                        std::vector<std::string> nbcs,
                        quadUtils& package);                        

DvD AssembleFVec(QTM::QuadTreeMesh& mesh, Real_b f, std::string evalStr);;
std::vector<Real_b> ComputeResiduals(QTM::QuadTreeMesh&, DvD& solution, SpD& source);
std::vector<std::shared_ptr<QTM::Cell>> TestResiduals(DvD& solution, QTM::QuadTreeMesh& mesh, Real_b residualLimit);
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
                Real_b c,
                Real_b k,
                std::string source,
                std::vector<std::string> bcs,
                Real_b penaltyParam);
DvD dgPoissonSolve(QTM::QuadTreeMesh& inputMesh,
                Real_b k,
                std::string source,
                std::vector<bool> isDirichletBC,
                std::vector<bool> isNeumannBC,
                std::vector<std::string> dbcs,
                std::vector<std::string> nbcs,
                Real_b penaltyParam,
                Real_b dirichletPenalty,
                uint64_t numThreads);

Real_b dgComputeResidual(QTM::QuadTreeMesh& inputMesh,
                Real_b k,
                std::string source,
                std::vector<bool> isDirichletBC,
                std::vector<bool> isNeumannBC,
                std::vector<std::string> dbcs,
                std::vector<std::string> nbcs,
                Real_b penaltyParam,
                Real_b dirichletPenalty);
}
#endif
