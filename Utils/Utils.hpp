#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

#include "../Dependencies/MathParser/MathParser.hpp"

#define PI 3.1415926535897932384

#ifndef UTILS
#define UTILS
 
// #define USE_SINGLE_PRECISION

#ifdef USE_SINGLE_PRECISION
using Real_b = float;
static_assert(sizeof(Real_b) == 4);
#else
using Real_b = double;
static_assert(sizeof(Real_b) == 8);
#endif
 
namespace Utils {
    std::vector<Real_b> evalLagrangeInterp(int k, std::vector<Real_b> evalPoints, std::vector<Real_b> &gaussPoints);
    std::vector<Real_b> genGaussPoints(int degree);
    std::vector<Real_b> numDeriv(Real_b h, int k, std::vector<Real_b>& evalPoints, std::vector<Real_b>& gaussPoints);
    std::vector<Real_b> numDeriv2(Real_b h, int k, std::vector<Real_b>& evalPoints, std::vector<Real_b>& gaussPoints);
    std::vector<Real_b> integrateLagrange(std::vector<Real_b> &gaussPoints);
    std::vector<Real_b> ReshuffleNodeVals(std::vector<int> RmOrder, std::vector<int> CwOrder, std::vector<Real_b> shuffleArray);
    std::vector<Real_b> EvalSymbolicBC(std::array<Real_b,2>* startpoint, int allocSize, std::string evalString);
    std::vector<Real_b> EvalSymbolicBC1D(Real_b* startpoint, int allocSize, std::string evalString);
    Real_b TransformPoint(Real_b x, Real_b a, Real_b b);
    std::array<Real_b, 2> GetNormalVector(std::array<Real_b, 2> &pos1, std::array<Real_b, 2> &pos2); 
    std::vector<int> GetBoundaryNodes(std::vector<int>& localElemNodes, std::vector<int>& offsets, int boundarySize);
}

namespace Utils {
    template<typename T>
    void printVec(std::vector<T>& in) {
        std::cout<<"--------BEGIN VECTOR---------"<<std::endl;
        std::cout<<"vector size: "<< in.size()<<std::endl;
        for (int i=0; i < in.size(); i++) {
            std::cout<<i<<": "<<in[i]<<std::endl;
        }
        std::cout<<"---------END VECTOR--------"<<std::endl;
    }
}

#endif