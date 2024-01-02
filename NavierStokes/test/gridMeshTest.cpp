#include "../../Dependencies/Meshing/include/BasicMesh2D.hpp"

#include <iostream>
#include <vector>

int main() {
    std::cout<<__cplusplus<<std::endl;
    std::vector<double> xdivs;
    int nxElem = 4;
    xdivs.reserve(nxElem);
    for (int i = 0; i < nxElem; i++) {
        xdivs.push_back(2);
    }

    std::vector<double> ydivs;
    int nyElem = 3;
    ydivs.reserve(nyElem);

    for (int i = 0; i < nyElem; i++) {
        ydivs.push_back(2);
    }

    int xdeg = 4; int ydeg = 4;

    int widthX = nxElem*(xdeg+1);
    int widthY = nyElem*(ydeg+1);

    Meshing::BasicMesh::BasicMesh2D mesh(xdeg, ydeg, xdivs, ydivs, 0, 0);

    std::cout<<mesh.nNodes()<<std::endl<<mesh.nElements()<<std::endl;

    int nodeIdx; std::cin>>nodeIdx;
    Meshing::BasicMesh::Element elem = mesh.Elements[nodeIdx];

    std::cout<<elem.Nodes.size()<<std::endl;
    
    std::vector<std::array<double, 2>> nodePos = mesh.posInElem(nodeIdx);
    int ind = 0;
    for (int y = 0; y < ydeg+1; y++) {
        for (int x = 0; x < xdeg+1; x++) {
            // std::cout<<nodePos[ind][0]<<", "<<nodePos[ind][1]<<" | ";
            std::cout<<elem.Nodes[ind]<<" | ";
            ind++;
        }
        std::cout<<std::endl;
    }

    // for (int y = j*ydeg; y < (j+1)*ydeg+1; y++) {
    //     for (int x = i*xdeg; x < (i+1)*xdeg+1; x++) {
    //         std::cout<<elem.Nodes[y*widthX+x]<<" | ";
    //     }
    //     std::cout<<std::endl;
    // }


}