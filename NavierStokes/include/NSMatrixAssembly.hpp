#include "..\Diffusion\include\MatrixAssembly.hpp"

#ifndef MatrixAssembly_NS
#define MatrixAssembly_NS

// Incompressible
SpD QVMatrix(QTM::QuadTreeMesh& mesh, int diffDir);
SpD PressureFluxMatrix(QTM::QuadTreeMesh& mesh, int diffDir);
SpD ConvectionMatrix(QTM::QuadTreeMesh& mesh, double rho, DvD&& U, DvD&& V, DvD&& P);
SpD ConvectionMatrix(QTM::QuadTreeMesh& mesh, double rho, DvD& state);

DD IncompressibleStokesSolve(QTM::QuadTreeMesh& inputMesh,
                double rho,
                double mu,
                std::vector<std::string> source,
                std::vector<std::string> bcs,
                double penaltyParam);

DD IncompressibleNavierStokesSolve(QTM::QuadTreeMesh& inputMesh,
                double rho,
                double mu,
                std::vector<std::string> source,
                std::vector<std::string> bcs,
                double penaltyParam);

// Compressible









#endif