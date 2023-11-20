#include "..\Dependencies\Utils\Utils.hpp"
#include "..\Dependencies\Eigen\Core"
#include "..\Dependencies\Eigen\Sparse"

typedef Eigen::SparseMatrix<double> SpD;
// Dynamically-sized matrix of doubles
typedef Eigen::MatrixXd DD;
// Dynamically-sized vector of doubles
typedef Eigen::VectorXd DvD;

struct Node {
    int ID;
    double position;
};

struct Element {
    int ID;
    std::vector<int> nodes;
    std::vector<int> dofs;
};

struct BasicMesh1D {
    int interpolantOrder;
    int numElems; 
    double elemWidth;
    std::vector<double> spacing;
    std::vector<Element> elements;
    std::vector<Node> nodes;
};

BasicMesh1D CreateMesh(int interpolantOrder, int numElems, 
                        double elemWidth) {
    BasicMesh1D mesh;

    // Populate basic mesh props
    mesh.interpolantOrder = interpolantOrder;
    mesh.numElems = numElems;
    mesh.elemWidth = elemWidth;

    // Generate nodes for basis func
    auto elemIntervals = Utils::genGaussPoints(interpolantOrder);
    mesh.spacing = elemIntervals;

    std::vector<Node> nodes;
    std::vector<Element> elements;

    // Populate node list
    Node node;
    int NID = 0;
    nodes.reserve(numElems * interpolantOrder + 1);
    for (int i=0; i<numElems; i++) { // iter over elements
        double offset = elemWidth * i;
        for (int j=0; j<interpolantOrder; j++) { // iter from element left bound
            // Calculate pos and push
            node.ID = NID++;
            node.position = elemIntervals[j] + offset;
            nodes.push_back(node);
        }
    }
    // Add final node
    node.ID = NID;
    node.position = numElems * elemWidth;
    nodes.push_back(node);

    std::cout<<"verifying element construction"<<std::endl;

    std::vector<int> nodeIds; std::vector<int> dofs;
    nodeIds.reserve(interpolantOrder + 1);
    dofs.reserve(interpolantOrder + 1);
    int dofId = 0;
    Element elem; int EID = 0; elements.reserve(numElems);
    for (int i=0; i<numElems; i++) {
        elem.ID = EID++;
        for (int j=0; j<interpolantOrder + 1; j++) {
            std::cout<<i*interpolantOrder + j<<std::endl;
            nodeIds.push_back(i*interpolantOrder + j);
            dofs.push_back(dofId++);
        }
        elem.nodes = nodeIds;
        elem.dofs = dofs;
        elements.push_back(elem);
        nodeIds.clear(); dofs.clear();
    }

    std::cout<<"-----------"<<std::endl;
    for (auto node : nodes) {
        std::cout<<node.ID<<std::endl;
    }

    std::cout<<mesh.numElems<<std::endl;

    for (auto elem : elements) {
        std::cout<<elem.ID<<std::endl;
        for (int i=0; i<interpolantOrder + 1; i++) {
            std::cout<<elem.nodes[i]<<": "<<elem.dofs[i]<<", ";
        }
        std::cout<<std::endl;
    }

    mesh.nodes = nodes;
    mesh.elements = elements;

    return mesh;
}

std::vector<double> GetPosOfElemNOdes(BasicMesh1D& mesh, int elem) {
    std::vector<double> out;
    auto elemNodes = mesh.elements[elem].nodes;
    out.reserve(mesh.interpolantOrder + 1);
    for (int i=0; i<mesh.interpolantOrder + 1; i++) {
        out.push_back(mesh.nodes[elemNodes[i]].position);
    }
    return out;
}

DD GenerateQuadWeights(std::vector<double> &gpX, int numXNodes) {
    // Generate quadrature weight matrices
    // Get integrals of basis functions
    std::vector<double> integX = Utils::integrateLagrange(gpX);

    Eigen::Map<DvD> integXMat(integX.data(),numXNodes);

    DD weightMat = integXMat.asDiagonal();

    return weightMat;
}

void AssembleSystemMatrices(BasicMesh1D& mesh, SpD &MassMatrix, SpD &StiffnessMatrix) {
    int xdeg = mesh.interpolantOrder;
    int numElemNodes = xdeg + 1;
    int nNodes = mesh.nodes.size(); int nElements = mesh.elements.size();
    std::vector<double> gaussPointsX = Utils::genGaussPoints(xdeg);

    // Generate derivative matrices
    std::vector<double> AInitializer;
    AInitializer.reserve(numElemNodes * numElemNodes);

    double *AxInitIdx = AInitializer.data(); 

    // Generate derivatives for each basis function, copy to full array
    for (int k=0; k<numElemNodes; k++) { 
        std::vector<double> xPartials = Utils::numDeriv(.00001, k, gaussPointsX, gaussPointsX);
        std::copy(xPartials.begin(), xPartials.end(), AxInitIdx);
        AxInitIdx += numElemNodes;
    }
    Eigen::Map<DD> Ax(AInitializer.data(), numElemNodes, numElemNodes); 

    // Generate quadrature weight matrices
    DD weightMat = GenerateQuadWeights(gaussPointsX, numElemNodes);

    std::vector<Eigen::Triplet<double>> tripletListM;
    std::vector<Eigen::Triplet<double>> tripletListK;
    tripletListM.reserve(nElements * numElemNodes);
    tripletListK.reserve(nElements * numElemNodes * numElemNodes);

    DvD massVec = weightMat.diagonal();
    DvD localElemVecMass(numElemNodes);
    DD localElemMatK(numElemNodes, numElemNodes);

    // Add contributions from basis functions (stiffness matrix is block-diagonal)
    for (auto &elm : mesh.elements) {
        double Lx = mesh.elemWidth; // Jacobian factors
        // calculate local matrix
        localElemVecMass = massVec*Lx/2;
        localElemMatK = Ax*weightMat*Ax.transpose();
        
        // Get nodes in element
        std::vector<int> dofsInElem = elm.dofs;
        std::vector<int> nodesInElem = elm.nodes;
        // Generate i,j,v triplets
        for (int j=0; j<numElemNodes; j++) {
            tripletListM.emplace_back(dofsInElem[j], dofsInElem[j], localElemVecMass(j));
            for (int i=0; i<numElemNodes; i++) {
                tripletListK.emplace_back(dofsInElem[i],dofsInElem[j],localElemMatK(i,j));
            }
        }
    }

    // Add contributions from flux terms
    for (auto &elm : mesh.elements) {

    }

    // Declare and construct sparse matrix from triplets
    MassMatrix.setFromTriplets(tripletListM.begin(), tripletListM.end());
    StiffnessMatrix.setFromTriplets(tripletListK.begin(), tripletListK.end());
}

int main() {
    int deg;
    int numElem;
    double width;
    std::cout<<"Specify domain length: ";
    std::cin>>width;
    std::cout<<std::endl<<"Specify num elems: ";
    std::cin>>numElem;
    std::cout<<std::endl<<"Specify shape func degree: ";
    std::cin>>deg;
    std::cout<<std::endl;
    double spacing = width / numElem;
    BasicMesh1D mesh = CreateMesh(deg, numElem, spacing);

    int numDofs = numElem * (deg + 1);
    SpD MassMatrix(numDofs, numDofs); SpD StiffnessMatrix(numDofs, numDofs);
    AssembleSystemMatrices(mesh, MassMatrix, StiffnessMatrix);
    return 0;
}

