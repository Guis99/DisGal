#include "..\..\Utils\Utils.hpp"

#include <array>

#ifndef MESH_RECT_2D
#define MESH_RECT_2D

namespace Meshing {
    namespace BasicMesh {
        class Element {
            public:
                int EID;
                std::vector<int> Nodes;
                std::vector<double> Boundaries;
                std::vector<int> Neighbors;

                Element(int EID, std::vector<int> Nodes, std::vector<double> Boundaries);
                double getWidth();
                double getHeight();
        };

        class Node {
            public:
                int NID;
                std::array<double, 2> Position;
                int nClass;

                Node(int NID, std::array<double, 2> Position, int nClass);
        };

        class Face {
            public:
                int FID;
                std::vector<int> plusElm;
                std::vector<int> minusElm;
                std::array<int, 2> endPoints; // Node label
                std::array<double, 2> normal;

                Face(int FID, std::array<int, 2> endPoints);
        };

        class BasicMesh2D {
            public:
                int xdeg;
                int ydeg;
                std::vector<double> xdiv;
                std::vector<double> ydiv;
                std::vector<Node> Nodes;
                std::vector<Element> Elements;
                std::vector<double> xOffsets;
                std::vector<double> yOffsets;

                BasicMesh2D(int xdeg, int ydeg, std::vector<double> xdiv, std::vector<double> ydiv, double xstart, double ystart);

                int numNodes() { return Nodes.size(); }
                int numElems() { return Elements.size(); }

                std::vector<std::array<double, 2>> allNodePos();
                std::vector<std::array<double, 2>> posInElem(int ElemID);
                std::vector<std::array<double, 2>> posOfNodes(std::vector<int> NodeIds);
                std::vector<int> getBoundaryNodes();
                std::vector<int> getFreeNodes();
                int nNodes();
                int nElements();
                static double transformPoint(double x, double a, double b);
        };
    }
}

#endif 