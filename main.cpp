#include <nlohmann/json.hpp>
#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

struct Node
{
    double x;
};

struct Element
{
    int node1, node2;
    double E, A;
};

Eigen::Matrix2d elementStiffness(const Element& e, const std::vector<Node>& nodes)
{
    double x1 = nodes[e.node1].x;
    double x2 = nodes[e.node2].x;
    double L = x2 - x1;
    double k = e.E * e.A / L;

    Eigen::Matrix2d K;
    K << k, -k,
        -k, k;

    return K;
}

int main()
{
    std::ifstream inputFile("structure.json");
    nlohmann::json data = nlohmann::json::parse(inputFile);

    int numElements = data["numElements"];
    double barLength = data["barLength"];
    double E = data["E"];
    double A = data["A"];
    double force = data["force"];

    std::vector<Node> nodes;
    for (int i = 0; i <= numElements; i++)
    {
        nodes.push_back(Node{i * barLength});
    }

    std::vector<Element> elements;
    for (int i = 0; i < numElements; i++)
    {
        elements.push_back(Element{i, i + 1, E, A});
    }

    int numNodes = numElements + 1;

    Eigen::MatrixXd K_global(numNodes, numNodes);
    K_global.setZero();

    for (int e = 0; e < elements.size(); e++)
    {
        Eigen::Matrix2d Klocal = elementStiffness(elements[e], nodes);
        int n1 = elements[e].node1;

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                K_global(n1 + i, n1 + j) += Klocal(i, j);
            }
        }
    }
    std::cout << "K_global =\n" << K_global << "\n";

    Eigen::VectorXd F(numNodes);
    F.setZero();
    F(numNodes - 1) = force;

    double penalty = K_global(0, 0) * 1e8;
    K_global(0, 0) += penalty;

    Eigen::VectorXd u = K_global.ldlt().solve(F);
    std::cout << "u =\n" << u << "\n";

    for (int e = 0; e < elements.size(); e++)
    {
        int n1 = elements[e].node1;
        int n2 = elements[e].node2;

        double x1 = nodes[n1].x;
        double x2 = nodes[n2].x;
        double L = x2 - x1;

        double strain = (u(n2) - u(n1)) / L;
        double stress = elements[e].E * strain;

        std::cout << "Element " << e << " stress = " << stress << " Pa\n";
    }

    return 0;
}