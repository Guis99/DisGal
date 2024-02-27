#include "..\..\Diffusion\include\MatrixAssembly.hpp"

#ifndef MatrixAssembly_NS
#define MatrixAssembly_NS

// Incompressible
SpD QVMatrix(QTM::QuadTreeMesh& mesh, int diffDir);
SpD PressureFluxMatrix(QTM::QuadTreeMesh& mesh, int diffDir);
SpD ConvectionMatrix(QTM::QuadTreeMesh& mesh, double rho, DvD&& U, DvD&& V, DvD&& P);
SpD ConvectionMatrix(QTM::QuadTreeMesh& mesh, double rho, DvD& state);

SpD GeneralizedNeumannBC(QTM::QuadTreeMesh& mesh, double mu);
DvD EvalPartialDirichletBoundaryCond(QTM::QuadTreeMesh& inputMesh, std::vector<std::vector<int>>& boundaryNodes, std::vector<std::string>& strs);

DD IncompressibleStokesSolve(QTM::QuadTreeMesh& inputMesh,
                double rho,
                double mu,
                std::vector<std::string> source,
                std::vector<std::string> Ubcs,
                std::vector<std::string> Vbcs,
                std::vector<std::string> Pbcs,
                double penaltyParam);

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