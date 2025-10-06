#include <NavierStokes/NSMatrixAssembly.hpp>
#include <NavierStokes/StokesMatrixAssembly.hpp>

TrDV InteriorBodyIntegral_Stokes(QTM::QuadTreeMesh& mesh, std::vector<EquationInfo>& equations, std::unordered_map<std::string, FunctionInfo>& variables) {
    // grad(u1) dot grad(v1) - partial(v1, x) * p -> partial(u1, x) * partial(v1, x) + partial(u1, y) * partial(v1, y) - partial(v1, x) * p
    // grad(u2) dot grad(v2) - partial(v2, y) * p -> partial(u2, x) * partial(v2, x) + partial(u2, y) * partial(v2, y) - partial(v1, x) * p
    // u1 * partial(q, x) + u2 * partial(q, y))
}

// boundary integrals:
// alpha/h * jump(u1)*jump(v1)
// avg(p) * n1 * jump(u1)
// q * (n1 * avg(u1) + n2 * avg(u2))