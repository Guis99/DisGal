#include "..\..\Utils\Utils.hpp"

#include <array>
#include <unordered_map>

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

                Node();
                Node(int NID, std::array<double, 2> Position, int nClass);
        };

        class Face {
            public:
                int FID;
                std::vector<int> PlusElm;
                std::vector<int> MinusElm;
                
                std::unordered_map<int, std::vector<int>> PlusNodes;
                std::unordered_map<int, std::vector<int>> MinusNodes;
                
                Meshing::BasicMesh::Node EndPoint1;
                Meshing::BasicMesh::Node EndPoint2;

                std::array<double, 2> Normal; // Calculated from endpoints

                Face(int FID, Meshing::BasicMesh::Node &endPoint1, Meshing::BasicMesh::Node &endPoint2);
        };

        class BasicMesh2D {
            public:
                int xdeg;
                int ydeg;
                std::vector<double> xdiv;
                std::vector<double> ydiv;
                std::vector<Node> Nodes;
                std::vector<Element> Elements;
                std::vector<Face> Faces;
                
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