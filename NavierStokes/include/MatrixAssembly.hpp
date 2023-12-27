#include "..\..\Dependencies\Eigen\Core"
#include "..\..\Dependencies\Eigen\Sparse"
#include "..\..\Dependencies\unsupported\Eigen\KroneckerProduct"
#include "..\..\Meshing\Meshing.hpp"

#ifndef SOLVERS
#define SOLVERS

typedef Eigen::SparseMatrix<double> SpD;
// Dynamically-sized matrix of doubles
typedef Eigen::MatrixXd DD;
// Dynamically-sized vector of doubles
typedef Eigen::VectorXd DvD;

namespace Solvers {
    namespace MatrixAssembly {
        DD GenerateQuadWeights(std::vector<double> &gpX, std::vector<double> &gpY, int numXNodes, int numYNodes, int numElemNodes);
        void AssembleMatrices(Meshing::BasicMesh::BasicMesh2D &inputMesh,
                    SpD &MassMatrix,
                    SpD &StiffnessMatrix,
                    SpD &SourceVector,
                    double c,
                    double k,
                    double f);
        SpD MassMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double c);
        SpD StiffnessMatrix(Meshing::BasicMesh::BasicMesh2D &inputMesh, double k);
        SpD ConvectionMatrix(Meshing::BasicMesh::BasicMesh2D & inputMesh, DvD state);
        SpD AssembleFVec(Meshing::BasicMesh::BasicMesh2D &inputMesh, double f);
        void GetExtensionMatrices(Meshing::BasicMesh::BasicMesh2D &inputMesh,
                                        std::vector<int> &boundaryNodes, 
                                        std::vector<int> &freeNodes,
                                        SpD &nullSpace,
                                        SpD &columnSpace);
        DvD EvalDirichletBoundaryCond(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &boundaryNodes);
        DvD EvalNeumannBoundaryCond(Meshing::BasicMesh::BasicMesh2D &inputMesh, std::vector<int> &boundaryNodes);
        // Solve discretized PDE
        DvD ComputeSolutionStationary(SpD &StiffnessMatrix, SpD &fVec, SpD &columnSpace, SpD &nullSpace, DvD &boundaryVals);
        std::vector<DvD> ComputeSolutionTimeDependent1stOrder(SpD &StiffnessMatrix, 
                                        SpD &MassMatrix, 
                                        SpD &fVec, SpD &columnSpace, 
                                        SpD &nullSpace, DvD &boundaryVals, 
                                        DvD &initialCondition,
                                        double timeStep,
                                        int numTimeSteps);

        // Set up problem
        DD PoissonSolve(Meshing::BasicMesh::BasicMesh2D &inputMesh,
                double c,
                double k,
                double f);

    }
}

#endif