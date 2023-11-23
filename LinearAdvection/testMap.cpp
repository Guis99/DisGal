#include "..\Dependencies\Utils\Utils.hpp"
#include "..\Dependencies\Eigen\Core"
#include "..\Dependencies\Eigen\Sparse"

#include <iostream>
#include <fstream>

int main() {
    typedef Eigen::MatrixXd DD;

    int numElemNodes = 5;


    std::vector<double> AInitializer;
    AInitializer.reserve(numElemNodes * numElemNodes);

    double *AxInitIdx = AInitializer.data(); 

    for (int k=0; k<numElemNodes; k++) { 
        std::vector<double> xPartials = {0,1,2,3,4};
        std::copy(xPartials.begin(), xPartials.end(), AxInitIdx);
        AxInitIdx += numElemNodes;
    }
    Eigen::Map<DD> Ax(AInitializer.data(), numElemNodes, numElemNodes); 
    std::cout<<Ax<<std::endl;
    return 0;
}