#include "..\..\Diffusion\include\MatrixAssembly.hpp"

#ifndef MatrixAssembly_NS
#define MatrixAssembly_NS

// Incompressible
namespace NSMA {
    SpD PVMatrix(QTM::QuadTreeMesh& mesh, int diffDir);
    // "proper" formulation
    SpD PressureJumpMatrix(QTM::QuadTreeMesh& mesh, int diffDir); // c3i
    SpD PressureAvgMatrix(QTM::QuadTreeMesh& mesh, int diffDir); // bi3

    SpD PenaltyMatrix(QTM::QuadTreeMesh& mesh, double k, double alpha, PMA::quadUtils& package);
    SpD FluxMatrix(QTM::QuadTreeMesh& mesh, double k, PMA::quadUtils& package);

    // "bad" formulation
    SpD PressureFluxMatrix(QTM::QuadTreeMesh& mesh, int diffDir);

    // 3,3 block
    SpD PressurePenaltyMatrix(QTM::QuadTreeMesh& mesh, double mu);

    SpD ConvectionMatrix(QTM::QuadTreeMesh& mesh, double rho, DvD&& U, DvD&& V, DvD&& P);
    SpD ConvectionMatrix(QTM::QuadTreeMesh& mesh, double rho, DvD& state);

    SpD GeneralizedNeumannBC(QTM::QuadTreeMesh& mesh, double mu);
    DvD EvalPartialSymbolicBoundaryCond(QTM::QuadTreeMesh& inputMesh, 
                                        std::vector<std::vector<int>>& boundaryNodes, 
                                        std::vector<std::string>& strs, 
                                        std::vector<int>& allBoundaryNodes, 
                                        int offset);

    SpD AssembleStokesSystem(QTM::QuadTreeMesh& mesh, 
                                SpD& dgPoissonMat,
                                SpD& mat1, 
                                SpD& mat2,
                                SpD& mat3,
                                SpD& mat4);

    SpD AssembleStokesSystem(QTM::QuadTreeMesh& mesh, 
                                SpD& dgPoissonMat,
                                SpD& QVMatrixX, 
                                SpD& QVMatrixY);

    DvD AssembleStokesSource(QTM::QuadTreeMesh& mesh, 
                                DvD& XSource,
                                DvD& YSource);

    DvD IncompressibleStokesSolve(QTM::QuadTreeMesh& inputMesh,
                    double rho,
                    double mu,
                    double penaltyParam,
                    std::vector<std::string> source,
                    std::vector<std::string> Ubcs,
                    std::vector<std::string> Vbcs,
                    std::vector<std::string> Pbcs,
                    std::vector<std::string> Nbcs,
                    std::vector<std::vector<int>>& boundaryNodes,
                    std::vector<bool> velEss,
                    std::vector<bool> pressEss,
                    std::vector<bool> natBC);

    DD IncompressibleNavierStokesSolve(QTM::QuadTreeMesh& inputMesh,
                    double rho,
                    double mu,
                    std::vector<std::string> source,
                    std::vector<std::string> Ubcs,
                    std::vector<std::string> Vbcs,
                    std::vector<std::string> Pbcs,
                    double penaltyParam);
                    
    
// Compressible
#endif
}
