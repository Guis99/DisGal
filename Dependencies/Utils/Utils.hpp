#define PI 3.1415926535897932384

#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

#ifndef UTILS
#define UTILS
 
namespace Utils {
    std::vector<double> evalLagrangeInterp(int k, std::vector<double> evalPoints, std::vector<double> &gaussPoints);
    std::vector<double> genGaussPoints(int degree);
    std::vector<double> numDeriv(double h, int k, std::vector<double>& evalPoints, std::vector<double>& gaussPoints);
    std::vector<double> numDeriv2(double h, int k, std::vector<double>& evalPoints, std::vector<double>& gaussPoints);
    std::vector<double> integrateLagrange(std::vector<double> &gaussPoints);
    std::vector<double> ReshuffleNodeVals(std::vector<int> RmOrder, std::vector<int> CwOrder, std::vector<double> shuffleArray);
    std::vector<double> EvalSymbolicBC(double* startpoint, int allocSize, std::string prompt);
    std::vector<double> EvalSymbolicBC(std::array<double,2>* startpoint, int allocSize, std::string evalString);
    std::vector<double> EvalSymbolicBC1D(double* startpoint, int allocSize, std::string evalString);
    double TransformPoint(double x, double a, double b);
    std::array<double, 2> GetNormalVector(std::array<double, 2> &pos1, std::array<double, 2> &pos2); 
    void printVec(std::vector<double> in);
    void printVec(std::vector<int> in);
    std::vector<int> GetBoundaryNodes(std::vector<int>& localElemNodes, std::vector<int>& offsets, int boundarySize);
}

#endif