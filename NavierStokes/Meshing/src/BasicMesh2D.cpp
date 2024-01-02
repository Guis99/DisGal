#include "..\include\BasicMesh2D.hpp"

#include <iostream>

Meshing::BasicMesh::Element::Element(int EID, std::vector<int> Nodes, std::vector<double> Boundaries) {
    this->EID = EID;
    this->Nodes = Nodes;
    this->Boundaries = Boundaries;
}

double Meshing::BasicMesh::Element::getWidth() {
    return Boundaries[1] - Boundaries[0];
}

double Meshing::BasicMesh::Element::getHeight() {
    return Boundaries[3] - Boundaries[2];
}

Meshing::BasicMesh::Node::Node(int NID, std::array<double, 2> Position, int nClass) {
    this->NID = NID;
    this->Position = Position;
    this->nClass = nClass;
}

Meshing::BasicMesh::BasicMesh2D::BasicMesh2D(int xdeg, int ydeg, 
                    std::vector<double> xdiv, std::vector<double> ydiv, 
                    double xstart, double ystart) {
    this->xdeg = xdeg;
    this->ydeg = ydeg;
    this->ydiv = ydiv;
    this->xdiv = xdiv;

    int nElemX = xdiv.size();
    int nElemY = ydiv.size();

    xOffsets.resize(nElemX+1);
    yOffsets.resize(nElemY+1);

    xOffsets[0] = xstart;
    yOffsets[0] = ystart;

    for (int i = 1; i < nElemX+1; i++) {
        xOffsets[i] = xOffsets[i-1]+xdiv[i-1];
    }

    for (int i = 1; i < nElemY+1; i++) {
        yOffsets[i] = yOffsets[i-1]+ydiv[i-1];
    }

    std::vector<double> xSpacing = Utils::genGaussPoints(xdeg);
    std::vector<double> ySpacing = Utils::genGaussPoints(ydeg);

    int widthX = nElemX*(xdeg+1);
    int widthY = nElemY*(ydeg+1);

    Elements.reserve(nElemX*nElemY);
    Nodes.reserve(widthX*widthY);

    int NID = 0; int nClass;
    int EID = 0;

    for (int j=0; j<nElemX; j++) {
        for (int i=0; i<nElemY; i++) { // loop over elements
            std::vector<int> currNodes;
            currNodes.reserve((xdeg+1)*(ydeg+1));

            double ax = xOffsets[i]; double bx = xOffsets[i+1];
            double ay = yOffsets[j]; double by = yOffsets[j+1];
            for (int jj=0; jj<ydeg+1; jj++) {
                for (int ii=0; ii<xdeg+1; ii++) { // loop over nodes/dofs to add
                    currNodes.push_back(NID);
                    std::array<double,2> pos;
                    pos[0] = transformPoint(xSpacing[ii], ax, bx);
                    pos[1] = transformPoint(ySpacing[jj], ay, by);

                    if (j==0 || j==ydeg + 1 || i==0 || i==xdeg) {
                        nClass = 1;
                    }
                    else {
                        nClass = 0;
                    }

                    Nodes.emplace_back(NID, pos, nClass);
                    NID++;
                }
            }

            std::vector<double> bounds = {ax, bx, ay, by};
            Elements.emplace_back(EID, currNodes, bounds);
            EID++;
        }
    }
}


std::vector<std::array<double, 2>> Meshing::BasicMesh::BasicMesh2D::allNodePos() {
    std::vector<std::array<double, 2>> out;
    out.reserve(Nodes.size());

    for (const auto &node : Nodes) {
        out.push_back(node.Position);
    }

    return out;
}

std::vector<std::array<double, 2>> Meshing::BasicMesh::BasicMesh2D::posInElem(int ElemID) {
    std::vector<std::array<double, 2>> out;
    out.reserve((xdeg+1)*(ydeg+1));

    for (const auto &node : Elements[ElemID].Nodes) {
        out.push_back(Nodes[node].Position);
    }

    return out;
}

std::vector<std::array<double, 2>> Meshing::BasicMesh::BasicMesh2D::posOfNodes(std::vector<int> NodeIds) {
    std::vector<std::array<double, 2>> out;
    out.reserve(NodeIds.size());

    for (const auto nodeID : NodeIds) {
        out.push_back(Nodes[nodeID].Position);
    }

    return out;
}

std::vector<int> Meshing::BasicMesh::BasicMesh2D::getBoundaryNodes() {
    int numXElems = xOffsets.size()-1;
    int numYElems = yOffsets.size()-1;
    int xWidth = xdeg*numXElems; int yWidth = ydeg*numYElems;

    int numBoundaryNodes = 2*xWidth + 2*yWidth;
    std::vector<int> boundaryNodes(numBoundaryNodes,0);

    for (int i=0; i<xWidth; i++) {
        boundaryNodes[i] = i;
        boundaryNodes[i+xWidth+yWidth] = this->nNodes() - xWidth + i;
    }

    for (int i=0; i<yWidth; i++) {
        boundaryNodes[i+xWidth] = i*(xWidth+1)+xWidth;
        boundaryNodes[i+2*xWidth+yWidth] = (i+1)*(xWidth+1);
    }

    return boundaryNodes;
}

std::vector<int> Meshing::BasicMesh::BasicMesh2D::getFreeNodes() {
    int numXElems = xOffsets.size()-1;
    int numYElems = yOffsets.size()-1;
    int xWidth = xdeg*numXElems; int yWidth = ydeg*numYElems;

    int numBoundaryNodes = 2*xWidth + 2*yWidth;
    std::vector<int> freeNodes; freeNodes.reserve(this->nNodes() - numBoundaryNodes);

    for (auto &node : Nodes) {
        if (node.nClass == 0) {
            freeNodes.push_back(node.NID);
        }
    }

    return freeNodes;
}

int Meshing::BasicMesh::BasicMesh2D::nNodes() {
    return Nodes.size();
}

int Meshing::BasicMesh::BasicMesh2D::nElements() {
    return Elements.size();
}

double Meshing::BasicMesh::BasicMesh2D::transformPoint(double x, double a, double b) {
    double result = a + ((b - a) / 2.0) * (x + 1.0);
    return result;
}