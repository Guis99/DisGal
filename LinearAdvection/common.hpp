#include "..\Dependencies\Utils\Utils.hpp"
#include "..\Dependencies\Eigen\Core"
#include "..\Dependencies\Eigen\Sparse"

#include <iostream>
#include <fstream>
#include <functional>

typedef Eigen::SparseMatrix<double> SpD;
// Dynamically-sized matrix of doubles
typedef Eigen::MatrixXd DD;
// Dynamically-sized vector of doubles
typedef Eigen::VectorXd DvD;

struct Node {
    int ID;
    double position;
};

struct Element {
    int ID;
    std::vector<int> nodes;
    std::vector<int> dofs;
};

struct BasicMesh1D {
    int interpolantOrder;
    int numElems; 
    double elemWidth;
    std::vector<double> spacing;
    std::vector<Element> elements;
    std::vector<Node> nodes;
};

BasicMesh1D CreateMesh(int interpolantOrder, int numElems, 
                        double elemWidth);

std::vector<double> GetPosOfElemNOdes(BasicMesh1D& mesh, int elem);
DD GenerateQuadWeights(std::vector<double> &gpX, int numXNodes);
void GetExtensionMatrices(BasicMesh1D &inputMesh,
                                        std::vector<int> &boundaryNodes, 
                                        std::vector<int> &freeNodes,
                                        SpD &nullSpace,
                                        SpD &columnSpace);

namespace TimeStep {
    
    std::vector<DvD> solver_RK4(SpD &A, 
                                SpD &columnSpace, SpD &nullSpace, 
                                DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps);

    std::vector<DvD> solver_RK4_NonLinear_Burger(SpD &A, 
                                SpD &columnSpace, SpD &nullSpace, 
                                DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps);

    std::vector<DvD> solver_GL1(SpD &A, 
                                SpD &columnSpace, SpD &nullSpace, 
                                DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps);

    std::vector<DvD> solver_GL2(SpD &A, 
                                SpD &columnSpace, SpD &nullSpace, 
                                DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps);

    std::vector<DvD> solver_GL2_nonlinear(SpD &A, 
                                SpD &columnSpace, SpD &nullSpace, 
                                DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps);
}

void AssembleSystemMatrices(BasicMesh1D& mesh, SpD &MassMatrix, SpD &StiffnessMatrix);

std::vector<DvD> ComputeTransientSolution(SpD &StiffnessMatrix, 
                                SpD &MassMatrix, SpD &columnSpace, 
                                SpD &nullSpace, DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps,
                                int integrator);

std::vector<DvD> ComputeTransientSolutionBurger(SpD &StiffnessMatrix, 
                                SpD &MassMatrix, SpD SIPMatrix,
                                SpD &columnSpace, 
                                SpD &nullSpace, DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps,
                                int integrator);