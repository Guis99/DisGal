#include "..\Dependencies\Utils\Utils.hpp"
#include "..\Dependencies\Eigen\Core"
#include "..\Dependencies\Eigen\Sparse"

#include <iostream>
#include <fstream>

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
        for (int j=0; j<interpolantOrder + 1; j++) { // iter from element left bound
            // Calculate pos and push
            node.ID = NID++;
            node.position = Utils::TransformPoint(elemIntervals[j], offset, offset+elemWidth);
            nodes.push_back(node);
        }
    }

    std::vector<int> nodeIds; std::vector<int> dofs;
    nodeIds.reserve(interpolantOrder + 1);
    dofs.reserve(interpolantOrder + 1);
    int dofId = 0;
    Element elem; int EID = 0; elements.reserve(numElems);
    for (int i=0; i<numElems; i++) {
        elem.ID = EID++;
        for (int j=0; j<interpolantOrder + 1; j++) {
            nodeIds.push_back(i*(interpolantOrder+1) + j);
            dofs.push_back(dofId++);
        }
        elem.dofs = dofs;
        elements.push_back(elem);
        dofs.clear();
    }

    // std::cout<<"-----------"<<std::endl;
    // for (auto node : nodes) {
    //     std::cout<<node.ID<<std::endl;
    // }

    // std::cout<<mesh.numElems<<std::endl;

    mesh.nodes = nodes;
    mesh.elements = elements;

    // for (auto elem : elements) {
    //     std::cout<<elem.ID<<std::endl;
    //     for (int i=0; i<interpolantOrder + 1; i++) {
    //         std::cout<<elem.dofs[i]<<": "<<mesh.nodes[elem.dofs[i]].position<<", ";
    //     }
    //     std::cout<<std::endl;
    // }

    return mesh;
}

std::vector<double> GetPosOfElemNOdes(BasicMesh1D& mesh, int elem) {
    std::vector<double> out;
    auto elemNodes = mesh.elements[elem].dofs;
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
        localElemMatK = -Ax*weightMat;
        
        // Get nodes in element
        std::vector<int> dofsInElem = elm.dofs;
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
        int elmID = elm.ID;
        int leftVal = elm.dofs[0];
        auto rightVal = (elm.dofs).back(); 
        if (elmID == 0) {
            tripletListK.emplace_back(rightVal, rightVal, 1.0);
        }
        else {
            tripletListK.emplace_back(leftVal, leftVal - 1, -1.0); // left boundary
            tripletListK.emplace_back(rightVal, rightVal, 1.0); // right boundary
        }
    }

    // Declare and construct sparse matrix from triplets
    MassMatrix.setFromTriplets(tripletListM.begin(), tripletListM.end());
    StiffnessMatrix.setFromTriplets(tripletListK.begin(), tripletListK.end());
}

void GetExtensionMatrices(BasicMesh1D &inputMesh,
                                        std::vector<int> &boundaryNodes, 
                                        std::vector<int> &freeNodes,
                                        SpD &nullSpace,
                                        SpD &columnSpace) {
    int nNodes = (inputMesh.nodes).size();

    std::sort(boundaryNodes.begin(), boundaryNodes.end());

    std::vector<Eigen::Triplet<double>> tripletListNS;
    std::vector<Eigen::Triplet<double>> tripletListCS;
    tripletListNS.reserve(boundaryNodes.size());
    tripletListCS.reserve(freeNodes.size());

    for (int i=0; i<boundaryNodes.size(); i++) {
        tripletListNS.emplace_back(boundaryNodes[i], i, 1.0);
    }

    for (int i=0; i<freeNodes.size(); i++) {
        tripletListCS.emplace_back(freeNodes[i], i, 1.0);
    }
    
    nullSpace.setFromTriplets(tripletListNS.begin(), tripletListNS.end());
    columnSpace.setFromTriplets(tripletListCS.begin(), tripletListCS.end());
}

std::vector<DvD> ComputeTransientSolution(SpD &StiffnessMatrix, 
                                SpD &MassMatrix, SpD &columnSpace, 
                                SpD &nullSpace, DvD &boundaryVals, 
                                DvD &initialCondition,
                                double timeStep,
                                int numTimeSteps) {
    std::cout<<"here1"<<std::endl;
    SpD K11 = columnSpace.transpose() * StiffnessMatrix * columnSpace;
    std::cout<<"here1"<<std::endl;
    SpD M11 = columnSpace.transpose() * MassMatrix * columnSpace;
    std::cout<<"here1"<<std::endl;
    // Eliminate boundary rows and free columns
    SpD K12 = columnSpace.transpose() * StiffnessMatrix * nullSpace;
    std::cout<<"here1"<<std::endl;
    SpD combinedMats = timeStep * K11 + M11;
    std::cout<<"here1"<<std::endl;
    int system_size = combinedMats.rows();
    std::cout<<"here1"<<std::endl;
    SpD dummyId(system_size, system_size); 
    std::cout<<"here1"<<std::endl;
    std::vector<Eigen::Triplet<double>> tripletListID;
    std::cout<<"here1"<<std::endl;

    tripletListID.reserve(combinedMats.rows());
    std::cout<<"here1"<<std::endl;

    for (int i=0; i<combinedMats.rows(); i++) {
        tripletListID.emplace_back(i, i, 1.0);
    }
    
    std::cout<<"here2"<<std::endl;
    dummyId.setFromTriplets(tripletListID.begin(), tripletListID.end());
    std::cout<<"here1"<<std::endl;
    Eigen::SparseLU<SpD, Eigen::COLAMDOrdering<int>> LuSolver;    
    LuSolver.analyzePattern(combinedMats);
    LuSolver.factorize(combinedMats);
    std::cout<<"here1"<<std::endl;
    SpD combinedMatsInv = LuSolver.solve(dummyId);

    std::vector<DvD> out;
    out.reserve(numTimeSteps);
    std::cout<<"here1"<<std::endl;
    std::cout<<"col: "<<columnSpace<<std::endl;
    std::cout<<"ic: "<<initialCondition;
    std::cout<<"null"<<nullSpace<<std::endl<<"boundary"<<boundaryVals<<std::endl;
    out.push_back(columnSpace * initialCondition + nullSpace * boundaryVals);
    DvD prevState = initialCondition;
    DvD x;

    // Time-stepping
    std::cout<<"here1"<<std::endl;
    for (int i=1; i<numTimeSteps; i++) {
        x = combinedMatsInv * (M11 * prevState - timeStep * K12 * boundaryVals);
        prevState = x;
        out.push_back(columnSpace * x + nullSpace * boundaryVals);
    }
    return out;
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

    double timeStep;
    int numTimeSteps;

    std::cout<<"Specify time-step: ";
    std::cin>>timeStep;
    std::cout<<std::endl<<"Specify number of time-steps: ";
    std::cin>>numTimeSteps;

    std::cout<<std::endl<<"Generating mesh"<<std::endl;
    BasicMesh1D mesh = CreateMesh(deg, numElem, spacing);
    auto nodes = mesh.nodes;
    int nNodes = nodes.size();

    std::vector<double> ICVec;
    ICVec.reserve(nNodes - 1);

    double currPos;

    for (int i=1; i<nNodes; i++) {
        currPos = nodes[i].position;
        if (currPos<.5) {
            ICVec.push_back(-16 * currPos * (currPos - .5));
        }
        else {
            ICVec.push_back(0);
        }
    }

    Eigen::Map<DvD> InitialCondition(ICVec.data(), nNodes - 1, 1);
    DvD initialCondition = (DvD)InitialCondition;

    std::vector<int> boundaryNodes = {0};
    std::vector<int> freeNodes; freeNodes.reserve(nNodes - 1);

    for (int i=1; i<nNodes; i++) {
        freeNodes.push_back(i);
    }

    int numDofs = numElem * (deg + 1);
    std::cout<<"Generating system matrices"<<std::endl;
    SpD MassMatrix(numDofs, numDofs); SpD StiffnessMatrix(numDofs, numDofs);
    AssembleSystemMatrices(mesh, MassMatrix, StiffnessMatrix);

    SpD nullSpace(nNodes, boundaryNodes.size());
    SpD columnSpace(nNodes, freeNodes.size());
    std::cout<<"Assembling null-orth matrix"<<std::endl;
    GetExtensionMatrices(mesh, boundaryNodes, freeNodes, nullSpace, columnSpace);

    std::vector<double> oneVal = {0};
    Eigen::Map<DvD> boundaryValsMap(oneVal.data(), 1,  1);
    DvD boundaryVals = (DvD) boundaryValsMap;

    std::cout<<"Computing solution"<<std::endl;
    std::vector<DvD> solns = ComputeTransientSolution(StiffnessMatrix, MassMatrix, 
                                                        columnSpace, 
                                                        nullSpace, boundaryVals, 
                                                        initialCondition, timeStep, numTimeSteps);

    
    
    std::vector<double> xPositions;
    xPositions.reserve(nNodes);
    for (int i=0; i<nNodes; i++) {
        xPositions.push_back(nodes[i].position);
    }

    Eigen::Map<DD> xOffsets(xPositions.data(), nNodes, 1);
    std::ofstream fileX("xt.txt");
    std::ofstream fileZ("zt.txt");

    if (fileZ.is_open())
    {
        for (auto x : solns) {
            fileZ << x;
            fileZ << "\n";
        }
    }
    if (fileX.is_open())
    {
        fileX << xOffsets;
    }

    return 0;
}