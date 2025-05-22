#include "Utils.hpp"

std::vector<Real_b> Utils::genGaussPoints(int degree) {
    std::vector<Real_b> gaussPoints;
    gaussPoints.reserve(degree+1);

    for (int i=0; i < degree+1; i++) {
        gaussPoints.push_back(-cos(i*PI/degree));
    }

    return gaussPoints;
}

std::vector<Real_b> Utils::evalLagrangeInterp(int k, std::vector<Real_b> evalPoints, std::vector<Real_b> &gaussPoints) {
    // Evaluate the kth Lagrange interpolant at the given locations
    int nNodes = gaussPoints.size();
    int nPoints = evalPoints.size();

    std::vector<Real_b> out;
    out.reserve(nPoints);  

    Real_b curr;
    for (int i = 0; i < nPoints; i++) {
        curr = 1.0;
        Real_b point = evalPoints[i];
        for (int j = 0; j < nNodes; j++) {
            if (j != k) {
            curr *= (point-gaussPoints[j])/(gaussPoints[k]-gaussPoints[j]);
            }
        }
        out.push_back(curr);
    }

    return out;
}

std::vector<Real_b> Utils::numDeriv(Real_b h, int k, std::vector<Real_b>& evalPoints, std::vector<Real_b>& gaussPoints) {
    // Fourth-order four-point derivative approximation given by
    // (-f(x + 2h) + 8f(x + h) − 8f(x -h) + f(x 2h))/12h

    // use step size h <= .001 for best results

    std::vector<Real_b> out;
    out.reserve(evalPoints.size());
    Real_b denom = 12*h;

    for (int i = 0; i < evalPoints.size(); i++) {
        Real_b currPoint = evalPoints[i];
        std::vector<Real_b> centeredPoints = {currPoint+2*h,currPoint+h,currPoint-h,currPoint-2*h};
        std::vector<Real_b> funcVals = Utils::evalLagrangeInterp(k, centeredPoints, gaussPoints);

        out.push_back((-funcVals[0]+8*funcVals[1]-8*funcVals[2]+funcVals[3])/denom);
    }
    
    return out;
}

std::vector<Real_b> Utils::numDeriv2(Real_b h, int k, std::vector<Real_b>& evalPoints, std::vector<Real_b>& gaussPoints) {
    // 4th-order centered second derivative given by 
    // (-f(x + 2h) + 16f(x + h) - 30f(x) + 16f(x - h) - f(x - 2h))/(12h^2)

    // use step size h <= .001 for best results
    if (gaussPoints.size() == 2) {
        std::vector<Real_b> out(2,0);
        return out;            
    }
    std::vector<Real_b> out;
    out.reserve(evalPoints.size());
    Real_b denom = 12*h*h;

    for (int i = 0; i < evalPoints.size(); i++) {
        Real_b currPoint = evalPoints[i];
        std::vector<Real_b> centeredPoints = {currPoint+2*h,currPoint+h,currPoint,currPoint-h,currPoint-2*h};
        std::vector<Real_b> funcVals = Utils::evalLagrangeInterp(k, centeredPoints, gaussPoints);

        out.push_back((-funcVals[0]+16*funcVals[1]-30*funcVals[2]+16*funcVals[3]-funcVals[4])/denom);
    }

    return out;
}

std::vector<Real_b> Utils::integrateLagrange(std::vector<Real_b> &gaussPoints) {
    // Integrates given basis polynomials with Simpson's rule
    int nInt = gaussPoints.size();

    std::vector<Real_b> evalPoints;
    std::vector<Real_b> out;
    int numDiv = 1000;
    Real_b h = 2/Real_b(numDiv);

    evalPoints.resize(numDiv);
    out.reserve(nInt);

    Real_b currLoc = -1;
    for (int i = 0; i < numDiv; i++) {
        evalPoints[i] = currLoc;
        currLoc += h;
    }
    evalPoints.push_back(1);

    Real_b currInt;
    for (int k = 0; k < nInt; k++) {
        currInt = 0.0;
        std::vector<Real_b> vals = Utils::evalLagrangeInterp(k, evalPoints, gaussPoints);
        for (int i = 0; i < evalPoints.size()-2; i+=2) {
            currInt += h*(vals[i]+4*vals[i+1]+vals[i+2]);
        }
        out.push_back(currInt/3);
    }

    return out;
}

std::vector<Real_b> Utils::ReshuffleNodeVals(std::vector<int> RmOrder, std::vector<int> CwOrder, std::vector<Real_b> shuffleArray) {
    int numObjs = CwOrder.size();
    std::vector<int> shuffledIdxs;
    shuffledIdxs.reserve(numObjs);
    for (int i=0; i<numObjs; i++) {
        int elmToMove = CwOrder[i];
        int lb = 0; int ub = numObjs - 1;
        int mid = (lb+ub)/2;
        while (elmToMove != RmOrder[mid]) {
            mid = (lb+ub)/2;
            if (RmOrder[mid] > elmToMove) {
                ub = mid;
            }
            else {
                lb = mid+1;
            }
        }
        shuffledIdxs.push_back(mid);
    }

    std::vector<Real_b> out;
    out.resize(numObjs);
    for (int i=0; i<numObjs; i++) {
        out[shuffledIdxs[i]] = shuffleArray[i];
    }

    return out;
}

std::vector<Real_b> Utils::EvalSymbolicBC1D(Real_b* startpoint, int allocSize, std::string evalString) {

    MathParser::QA x;
    x.reserve(allocSize);

    for (int i=0; i<allocSize; i++) {
        auto coord = *(startpoint+i);
        x.push_back(*(startpoint+i));
    }

    MathParser::InitMaps();
    MathParser::SetVariable("x", x);

    std::cout<<std::endl;
    auto result = MathParser::ParseText(evalString);
    
    std::vector<Real_b> out;

    if (result.size() == x.size()) {
        out = result.release();
    }
    else {
        out = std::vector<Real_b>(x.size(), result.release()[0]);
    }

    return out;
}

std::vector<Real_b> Utils::EvalSymbolicBC(std::array<Real_b,2>* startpoint, int allocSize, std::string evalString) {
    MathParser::QA x;
    MathParser::QA y;

    x.reserve(allocSize); y.reserve(allocSize);

    for (int i=0; i<allocSize; i++) {
        auto coord = *(startpoint+i);
        x.push_back(coord[0]);
        y.push_back(coord[1]);
    }

    MathParser::InitMaps();
    MathParser::SetVariable("x", x);
    MathParser::SetVariable("y", y);

    auto result = MathParser::ParseText(evalString);
    
    std::vector<Real_b> out;

    if (result.size() == x.size()) {
        out = result.release();
    }
    else {
        out = std::vector<Real_b>(x.size(), result.release()[0]);
    }

    return out;
}

Real_b Utils::TransformPoint(Real_b x, Real_b a, Real_b b) {
    Real_b result = a + ((b - a) / 2.0) * (x + 1.0);
    return result;
}

std::array<Real_b, 2> Utils::GetNormalVector(std::array<Real_b, 2> &pos1, std::array<Real_b, 2> &pos2) {
    auto xdiff = pos2[0]-pos1[0];
    auto ydiff = pos2[1]-pos1[1];
    Real_b vecNorm = std::sqrt(std::pow(xdiff,2) + std::pow(ydiff,2));

    std::array<Real_b,2> out;
    out[0] = -ydiff / vecNorm;
    out[1] = xdiff / vecNorm;

    return out;
}

std::vector<int> Utils::GetBoundaryNodes(std::vector<int>& localElemNodes, std::vector<int>& offsets, int boundarySize) {
    int numElemNodes = boundarySize * boundarySize;
    std::vector<int> out; out.reserve(offsets.size() * numElemNodes);
    for (int cellID : offsets) {
        int offset = numElemNodes * cellID;
        for (int localNodeID : localElemNodes) {
            out.push_back(localNodeID + offset);
        }
    }

    return out;
}