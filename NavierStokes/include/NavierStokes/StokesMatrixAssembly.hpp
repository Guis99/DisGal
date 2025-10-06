#include <Utils/Utils.hpp>

#ifndef MatrixAssembly_Stokes
#define MatrixAssembly_Stokes

enum FunctionType {
    TRIAL,
    TEST
};

struct FunctionInfo {
    std::string identifier;
    uint8_t variableIndex;
    uint8_t order;
    FunctionType functionType;
};

struct EquationInfo {
    std::string equation;
    std::unordered_map<std::string, Real_b> parameters;
};

// Incompressible
namespace NSMA {
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
    
    TrDV InteriorBodyIntegral_Stokes(QTM::QuadTreeMesh& mesh, std::vector<EquationInfo>& equations, std::unordered_map<std::string, FunctionInfo>& variables);
}

#endif