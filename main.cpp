#include <nlohmann/json.hpp>
#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>

struct Node
{
    double x, y;
};

struct Element
{
    int node1, node2;
    double E, A;
};

Eigen::Matrix4d elementStiffness(const Element& e, const std::vector<Node>& nodes)
{
    double x1 = nodes[e.node1].x;
    double y1 = nodes[e.node1].y;
    double x2 = nodes[e.node2].x;
    double y2 = nodes[e.node2].y;

    double dx = x2 - x1;
    double dy = y2 - y1;
    double L  = std::sqrt(dx*dx + dy*dy);
    double c  = dx / L;
    double s  = dy / L;
    double k  = e.E * e.A / L;

    Eigen::Matrix4d K;
    K << c*c,  c*s, -c*c, -c*s,
         c*s,  s*s, -c*s, -s*s,
        -c*c, -c*s,  c*c,  c*s,
        -c*s, -s*s,  c*s,  s*s;

    return k * K;
}

int main()
{
    std::ifstream inputFile("structure.json");
    nlohmann::json data = nlohmann::json::parse(inputFile);

    std::vector<Node> nodes;
    for (auto& jn : data["nodes"])
    {
        nodes.push_back(Node{jn["x"], jn["y"]});
    }

    std::vector<Element> elements;
    for (auto& je : data["elements"])
    {
        elements.push_back(Element{je["node1"], je["node2"], je["E"], je["A"]});
    }

    double force   = data["force"];
    int numNodes   = nodes.size();
    int numDOFs    = 2 * numNodes;

    Eigen::MatrixXd K_global(numDOFs, numDOFs);
    K_global.setZero();

    for (int e = 0; e < elements.size(); e++)
    {
        Eigen::Matrix4d Klocal = elementStiffness(elements[e], nodes);
        int n1 = elements[e].node1;
        int n2 = elements[e].node2;

        int dofs[4] = {2*n1, 2*n1+1, 2*n2, 2*n2+1};

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                K_global(dofs[i], dofs[j]) += Klocal(i, j);
            }
        }
    }

    Eigen::VectorXd F(numDOFs);
    F.setZero();
    F(numDOFs - 1) = force;

    

double penalty = K_global(0, 0) * 1e8;
K_global(0, 0) += penalty;   // node 0 u — pinned
K_global(1, 1) += penalty;   // node 0 v — pinned
K_global(3, 3) += penalty;   // node 1 v — roller

    Eigen::VectorXd u = K_global.ldlt().solve(F);
    std::cout << "u =\n" << u << "\n";

    for (int e = 0; e < elements.size(); e++)
    {
        int n1 = elements[e].node1;
        int n2 = elements[e].node2;

        double dx = nodes[n2].x - nodes[n1].x;
        double dy = nodes[n2].y - nodes[n1].y;
        double L  = std::sqrt(dx*dx + dy*dy);

        double u1 = u(2*n1),   v1 = u(2*n1+1);
        double u2 = u(2*n2),   v2 = u(2*n2+1);

        double elongation = (u2-u1)*dx/L + (v2-v1)*dy/L;
        double strain     = elongation / L;
        double stress     = elements[e].E * strain;

        std::cout << "Element " << e << " stress = " << stress << " Pa\n";
    }

    return 0;
}