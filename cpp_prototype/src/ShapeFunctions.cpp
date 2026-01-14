#include "ShapeFunctions.h"
#include <stdexcept>
#include <cmath>

using namespace elmer;

std::vector<double> ShapeFunctions::linear1D(double xi, int nNodes) {
    std::vector<double> N(nNodes, 0.0);
    
    if (nNodes == 2) {
        // Linear shape functions for 1D bar element
        N[0] = 0.5 * (1.0 - xi);
        N[1] = 0.5 * (1.0 + xi);
    } else if (nNodes == 3) {
        // Quadratic shape functions for 1D bar element
        N[0] = 0.5 * xi * (xi - 1.0);
        N[1] = 1.0 - xi * xi;
        N[2] = 0.5 * xi * (xi + 1.0);
    } else {
        throw std::invalid_argument("Unsupported number of nodes for 1D element");
    }
    
    return N;
}

std::vector<double> ShapeFunctions::linear1DDerivatives(double xi, int nNodes) {
    std::vector<double> dNdxi(nNodes, 0.0);
    
    if (nNodes == 2) {
        // Derivatives of linear shape functions
        dNdxi[0] = -0.5;
        dNdxi[1] = 0.5;
    } else if (nNodes == 3) {
        // Derivatives of quadratic shape functions
        dNdxi[0] = xi - 0.5;
        dNdxi[1] = -2.0 * xi;
        dNdxi[2] = xi + 0.5;
    } else {
        throw std::invalid_argument("Unsupported number of nodes for 1D element");
    }
    
    return dNdxi;
}

std::vector<double> ShapeFunctions::linearQuadrilateral(double xi, double eta, int nNodes) {
    std::vector<double> N(nNodes, 0.0);
    
    if (nNodes == 4) {
        // Bilinear shape functions for quadrilateral element
        N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
        N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
        N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
        N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
    } else if (nNodes == 8) {
        // Quadratic shape functions for quadrilateral element
        double xi2 = xi * xi;
        double eta2 = eta * eta;
        
        N[0] = 0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0);
        N[1] = 0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0);
        N[2] = 0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0);
        N[3] = 0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0);
        N[4] = 0.5 * (1.0 - xi2) * (1.0 - eta);
        N[5] = 0.5 * (1.0 + xi) * (1.0 - eta2);
        N[6] = 0.5 * (1.0 - xi2) * (1.0 + eta);
        N[7] = 0.5 * (1.0 - xi) * (1.0 - eta2);
    } else {
        throw std::invalid_argument("Unsupported number of nodes for quadrilateral element");
    }
    
    return N;
}

std::pair<std::vector<double>, std::vector<double>> 
ShapeFunctions::linearQuadrilateralDerivatives(double xi, double eta, int nNodes) {
    std::vector<double> dNdxi(nNodes, 0.0);
    std::vector<double> dNdeta(nNodes, 0.0);
    
    if (nNodes == 4) {
        // Derivatives of bilinear shape functions
        dNdxi[0] = -0.25 * (1.0 - eta);
        dNdxi[1] = 0.25 * (1.0 - eta);
        dNdxi[2] = 0.25 * (1.0 + eta);
        dNdxi[3] = -0.25 * (1.0 + eta);
        
        dNdeta[0] = -0.25 * (1.0 - xi);
        dNdeta[1] = -0.25 * (1.0 + xi);
        dNdeta[2] = 0.25 * (1.0 + xi);
        dNdeta[3] = 0.25 * (1.0 - xi);
    } else if (nNodes == 8) {
        // Derivatives of quadratic shape functions
        double xi2 = xi * xi;
        double eta2 = eta * eta;
        
        dNdxi[0] = 0.25 * (1.0 - eta) * (2.0 * xi + eta);
        dNdxi[1] = 0.25 * (1.0 - eta) * (2.0 * xi - eta);
        dNdxi[2] = 0.25 * (1.0 + eta) * (2.0 * xi + eta);
        dNdxi[3] = 0.25 * (1.0 + eta) * (2.0 * xi - eta);
        dNdxi[4] = -xi * (1.0 - eta);
        dNdxi[5] = 0.5 * (1.0 - eta2);
        dNdxi[6] = -xi * (1.0 + eta);
        dNdxi[7] = -0.5 * (1.0 - eta2);
        
        dNdeta[0] = 0.25 * (1.0 - xi) * (xi + 2.0 * eta);
        dNdeta[1] = 0.25 * (1.0 + xi) * (-xi + 2.0 * eta);
        dNdeta[2] = 0.25 * (1.0 + xi) * (xi + 2.0 * eta);
        dNdeta[3] = 0.25 * (1.0 - xi) * (-xi + 2.0 * eta);
        dNdeta[4] = -0.5 * (1.0 - xi2);
        dNdeta[5] = -eta * (1.0 + xi);
        dNdeta[6] = 0.5 * (1.0 - xi2);
        dNdeta[7] = -eta * (1.0 - xi);
    } else {
        throw std::invalid_argument("Unsupported number of nodes for quadrilateral element");
    }
    
    return {dNdxi, dNdeta};
}

std::vector<double> ShapeFunctions::linearHexahedron(double xi, double eta, double zeta, int nNodes) {
    std::vector<double> N(nNodes, 0.0);
    
    if (nNodes == 8) {
        // Trilinear shape functions for hexahedral element
        N[0] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta);
        N[1] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta);
        N[2] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta);
        N[3] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta);
        N[4] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta);
        N[5] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta);
        N[6] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta);
        N[7] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta);
    } else if (nNodes == 20) {
        // Quadratic shape functions for hexahedral element
        double xi2 = xi * xi;
        double eta2 = eta * eta;
        double zeta2 = zeta * zeta;
        
        // Corner nodes
        N[0] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta) * (-xi - eta - zeta - 2.0);
        N[1] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta) * (xi - eta - zeta - 2.0);
        N[2] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta) * (xi + eta - zeta - 2.0);
        N[3] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta) * (-xi + eta - zeta - 2.0);
        N[4] = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta) * (-xi - eta + zeta - 2.0);
        N[5] = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta) * (xi - eta + zeta - 2.0);
        N[6] = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta) * (xi + eta + zeta - 2.0);
        N[7] = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta) * (-xi + eta + zeta - 2.0);
        
        // Edge nodes (midpoints)
        N[8] = 0.25 * (1.0 - xi2) * (1.0 - eta) * (1.0 - zeta);
        N[9] = 0.25 * (1.0 + xi) * (1.0 - eta2) * (1.0 - zeta);
        N[10] = 0.25 * (1.0 - xi2) * (1.0 + eta) * (1.0 - zeta);
        N[11] = 0.25 * (1.0 - xi) * (1.0 - eta2) * (1.0 - zeta);
        N[12] = 0.25 * (1.0 - xi2) * (1.0 - eta) * (1.0 + zeta);
        N[13] = 0.25 * (1.0 + xi) * (1.0 - eta2) * (1.0 + zeta);
        N[14] = 0.25 * (1.0 - xi2) * (1.0 + eta) * (1.0 + zeta);
        N[15] = 0.25 * (1.0 - xi) * (1.0 - eta2) * (1.0 + zeta);
        N[16] = 0.25 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta2);
        N[17] = 0.25 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta2);
        N[18] = 0.25 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta2);
        N[19] = 0.25 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta2);
    } else {
        throw std::invalid_argument("Unsupported number of nodes for hexahedral element");
    }
    
    return N;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
ShapeFunctions::linearHexahedronDerivatives(double xi, double eta, double zeta, int nNodes) {
    std::vector<double> dNdxi(nNodes, 0.0);
    std::vector<double> dNdeta(nNodes, 0.0);
    std::vector<double> dNdzeta(nNodes, 0.0);
    
    if (nNodes == 8) {
        // Derivatives of trilinear shape functions
        dNdxi[0] = -0.125 * (1.0 - eta) * (1.0 - zeta);
        dNdxi[1] = 0.125 * (1.0 - eta) * (1.0 - zeta);
        dNdxi[2] = 0.125 * (1.0 + eta) * (1.0 - zeta);
        dNdxi[3] = -0.125 * (1.0 + eta) * (1.0 - zeta);
        dNdxi[4] = -0.125 * (1.0 - eta) * (1.0 + zeta);
        dNdxi[5] = 0.125 * (1.0 - eta) * (1.0 + zeta);
        dNdxi[6] = 0.125 * (1.0 + eta) * (1.0 + zeta);
        dNdxi[7] = -0.125 * (1.0 + eta) * (1.0 + zeta);
        
        dNdeta[0] = -0.125 * (1.0 - xi) * (1.0 - zeta);
        dNdeta[1] = -0.125 * (1.0 + xi) * (1.0 - zeta);
        dNdeta[2] = 0.125 * (1.0 + xi) * (1.0 - zeta);
        dNdeta[3] = 0.125 * (1.0 - xi) * (1.0 - zeta);
        dNdeta[4] = -0.125 * (1.0 - xi) * (1.0 + zeta);
        dNdeta[5] = -0.125 * (1.0 + xi) * (1.0 + zeta);
        dNdeta[6] = 0.125 * (1.0 + xi) * (1.0 + zeta);
        dNdeta[7] = 0.125 * (1.0 - xi) * (1.0 + zeta);
        
        dNdzeta[0] = -0.125 * (1.0 - xi) * (1.0 - eta);
        dNdzeta[1] = -0.125 * (1.0 + xi) * (1.0 - eta);
        dNdzeta[2] = -0.125 * (1.0 + xi) * (1.0 + eta);
        dNdzeta[3] = -0.125 * (1.0 - xi) * (1.0 + eta);
        dNdzeta[4] = 0.125 * (1.0 - xi) * (1.0 - eta);
        dNdzeta[5] = 0.125 * (1.0 + xi) * (1.0 - eta);
        dNdzeta[6] = 0.125 * (1.0 + xi) * (1.0 + eta);
        dNdzeta[7] = 0.125 * (1.0 - xi) * (1.0 + eta);
    } else {
        throw std::invalid_argument("Quadratic hexahedral derivatives not implemented");
    }
    
    return {dNdxi, dNdeta, dNdzeta};
}

std::vector<double> ShapeFunctions::linearTriangle(double xi, double eta, int nNodes) {
    std::vector<double> N(nNodes, 0.0);
    
    if (nNodes == 3) {
        // Linear shape functions for triangular element (area coordinates)
        N[0] = 1.0 - xi - eta;
        N[1] = xi;
        N[2] = eta;
    } else if (nNodes == 6) {
        // Quadratic shape functions for triangular element
        double L1 = 1.0 - xi - eta;
        double L2 = xi;
        double L3 = eta;
        
        N[0] = L1 * (2.0 * L1 - 1.0);
        N[1] = L2 * (2.0 * L2 - 1.0);
        N[2] = L3 * (2.0 * L3 - 1.0);
        N[3] = 4.0 * L1 * L2;
        N[4] = 4.0 * L2 * L3;
        N[5] = 4.0 * L3 * L1;
    } else {
        throw std::invalid_argument("Unsupported number of nodes for triangular element");
    }
    
    return N;
}

std::pair<std::vector<double>, std::vector<double>> 
ShapeFunctions::linearTriangleDerivatives(double xi, double eta, int nNodes) {
    std::vector<double> dNdxi(nNodes, 0.0);
    std::vector<double> dNdeta(nNodes, 0.0);
    
    if (nNodes == 3) {
        // Derivatives of linear shape functions
        dNdxi[0] = -1.0;
        dNdxi[1] = 1.0;
        dNdxi[2] = 0.0;
        
        dNdeta[0] = -1.0;
        dNdeta[1] = 0.0;
        dNdeta[2] = 1.0;
    } else if (nNodes == 6) {
        // Derivatives of quadratic shape functions
        double L1 = 1.0 - xi - eta;
        double L2 = xi;
        double L3 = eta;
        
        dNdxi[0] = -4.0 * L1 + 1.0;
        dNdxi[1] = 4.0 * L2 - 1.0;
        dNdxi[2] = 0.0;
        dNdxi[3] = 4.0 * (L1 - L2);
        dNdxi[4] = 4.0 * L3;
        dNdxi[5] = -4.0 * L3;
        
        dNdeta[0] = -4.0 * L1 + 1.0;
        dNdeta[1] = 0.0;
        dNdeta[2] = 4.0 * L3 - 1.0;
        dNdeta[3] = -4.0 * L2;
        dNdeta[4] = 4.0 * L2;
        dNdeta[5] = 4.0 * (L1 - L3);
    } else {
        throw std::invalid_argument("Unsupported number of nodes for triangular element");
    }
    
    return {dNdxi, dNdeta};
}

std::vector<double> ShapeFunctions::linearTetrahedron(double xi, double eta, double zeta, int nNodes) {
    std::vector<double> N(nNodes, 0.0);
    
    if (nNodes == 4) {
        // Linear shape functions for tetrahedral element (volume coordinates)
        N[0] = 1.0 - xi - eta - zeta;
        N[1] = xi;
        N[2] = eta;
        N[3] = zeta;
    } else if (nNodes == 10) {
        // Quadratic shape functions for tetrahedral element
        double L1 = 1.0 - xi - eta - zeta;
        double L2 = xi;
        double L3 = eta;
        double L4 = zeta;
        
        N[0] = L1 * (2.0 * L1 - 1.0);
        N[1] = L2 * (2.0 * L2 - 1.0);
        N[2] = L3 * (2.0 * L3 - 1.0);
        N[3] = L4 * (2.0 * L4 - 1.0);
        N[4] = 4.0 * L1 * L2;
        N[5] = 4.0 * L2 * L3;
        N[6] = 4.0 * L3 * L1;
        N[7] = 4.0 * L1 * L4;
        N[8] = 4.0 * L2 * L4;
        N[9] = 4.0 * L3 * L4;
    } else {
        throw std::invalid_argument("Unsupported number of nodes for tetrahedral element");
    }
    
    return N;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
ShapeFunctions::linearTetrahedronDerivatives(double xi, double eta, double zeta, int nNodes) {
    std::vector<double> dNdxi(nNodes, 0.0);
    std::vector<double> dNdeta(nNodes, 0.0);
    std::vector<double> dNdzeta(nNodes, 0.0);
    
    if (nNodes == 4) {
        // Derivatives of linear shape functions
        dNdxi[0] = -1.0;
        dNdxi[1] = 1.0;
        dNdxi[2] = 0.0;
        dNdxi[3] = 0.0;
        
        dNdeta[0] = -1.0;
        dNdeta[1] = 0.0;
        dNdeta[2] = 1.0;
        dNdeta[3] = 0.0;
        
        dNdzeta[0] = -1.0;
        dNdzeta[1] = 0.0;
        dNdzeta[2] = 0.0;
        dNdzeta[3] = 1.0;
    } else {
        throw std::invalid_argument("Quadratic tetrahedral derivatives not implemented");
    }
    
    return {dNdxi, dNdeta, dNdzeta};
}

ShapeFunctions::ShapeResult ShapeFunctions::computeShapeFunctions(ElementType elementType, 
                                                                 const std::vector<Node>& nodes,
                                                                 double xi, double eta, double zeta) {
    int nNodes = nodes.size();
    ShapeResult result(nNodes);
    
    // Compute natural derivatives based on element type
    switch (elementType) {
        case ElementType::LINEAR:
            // For LINEAR type, determine dimensionality based on node count
            if (nNodes == 2) {
                // 1D linear element
                result = compute1DShapeFunctions(nodes, xi);
            } else if (nNodes == 3) {
                // Triangle element
                result = computeTriangleShapeFunctions(nodes, xi, eta);
            } else if (nNodes == 4) {
                // Quadrilateral element
                result = compute2DQuadrilateralShapeFunctions(nodes, xi, eta);
            } else if (nNodes == 8) {
                // Hexahedron element
                result = compute3DHexahedronShapeFunctions(nodes, xi, eta, zeta);
            } else {
                throw std::invalid_argument("Unsupported number of nodes for LINEAR element type");
            }
            break;
        case ElementType::HEXAHEDRON:
            result = compute3DHexahedronShapeFunctions(nodes, xi, eta, zeta);
            break;
        case ElementType::TETRAHEDRON:
            result = computeTetrahedronShapeFunctions(nodes, xi, eta, zeta);
            break;
        default:
            throw std::invalid_argument("Unsupported element type");
    }
    
    return result;
}

std::array<std::array<double, 3>, 3> 
ShapeFunctions::computeJacobianMatrix(const std::vector<Node>& nodes,
                                     const std::vector<double>& dNdxi,
                                     const std::vector<double>& dNdeta,
                                     const std::vector<double>& dNdzeta) {
    std::array<std::array<double, 3>, 3> jac = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    
    int nNodes = nodes.size();
    for (int i = 0; i < nNodes; ++i) {
        jac[0][0] += nodes[i].x * dNdxi[i];    // dx/dxi
        jac[0][1] += nodes[i].y * dNdxi[i];    // dy/dxi
        jac[0][2] += nodes[i].z * dNdxi[i];    // dz/dxi
        
        jac[1][0] += nodes[i].x * dNdeta[i];   // dx/deta
        jac[1][1] += nodes[i].y * dNdeta[i];   // dy/deta
        jac[1][2] += nodes[i].z * dNdeta[i];   // dz/deta
        
        jac[2][0] += nodes[i].x * dNdzeta[i];  // dx/dzeta
        jac[2][1] += nodes[i].y * dNdzeta[i];  // dy/dzeta
        jac[2][2] += nodes[i].z * dNdzeta[i];  // dz/dzeta
    }
    
    return jac;
}

std::array<std::array<double, 3>, 3> 
ShapeFunctions::computeInverseJacobian(const std::array<std::array<double, 3>, 3>& jac) {
    std::array<std::array<double, 3>, 3> invJac;
    
    double detJ = computeJacobianDeterminant(jac);
    
    if (std::abs(detJ) < 1e-15) {
        throw std::runtime_error("Jacobian determinant is zero or too small");
    }
    
    // Compute inverse using formula for 3x3 matrix
    invJac[0][0] = (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) / detJ;
    invJac[0][1] = (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) / detJ;
    invJac[0][2] = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) / detJ;
    
    invJac[1][0] = (jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2]) / detJ;
    invJac[1][1] = (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) / detJ;
    invJac[1][2] = (jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2]) / detJ;
    
    invJac[2][0] = (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]) / detJ;
    invJac[2][1] = (jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1]) / detJ;
    invJac[2][2] = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) / detJ;
    
    return invJac;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
ShapeFunctions::transformDerivatives(const std::array<std::array<double, 3>, 3>& invJac,
                                    const std::vector<double>& dNdxi,
                                    const std::vector<double>& dNdeta,
                                    const std::vector<double>& dNdzeta) {
    int nNodes = dNdxi.size();
    std::vector<double> dNdx(nNodes, 0.0);
    std::vector<double> dNdy(nNodes, 0.0);
    std::vector<double> dNdz(nNodes, 0.0);
    
    for (int i = 0; i < nNodes; ++i) {
        dNdx[i] = invJac[0][0] * dNdxi[i] + invJac[0][1] * dNdeta[i] + invJac[0][2] * dNdzeta[i];
        dNdy[i] = invJac[1][0] * dNdxi[i] + invJac[1][1] * dNdeta[i] + invJac[1][2] * dNdzeta[i];
        dNdz[i] = invJac[2][0] * dNdxi[i] + invJac[2][1] * dNdeta[i] + invJac[2][2] * dNdzeta[i];
    }
    
    return {dNdx, dNdy, dNdz};
}

double ShapeFunctions::computeJacobianDeterminant(const std::array<std::array<double, 3>, 3>& jac) {
    return jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) -
           jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) +
           jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
}

int ShapeFunctions::getNumberOfNodes(ElementType elementType, int order) {
    switch (elementType) {
        case ElementType::LINEAR:
            return order + 1;  // 2 for linear, 3 for quadratic
        case ElementType::HEXAHEDRON:
            return (order == 1) ? 8 : 20; // 8 for linear, 20 for quadratic
        case ElementType::TETRAHEDRON:
            return (order == 1) ? 4 : 10; // 4 for linear, 10 for quadratic
        default:
            return 0;
    }
}

bool ShapeFunctions::isSupported(ElementType elementType) {
    return (elementType == ElementType::LINEAR ||
            elementType == ElementType::HEXAHEDRON ||
            elementType == ElementType::TETRAHEDRON);
}

int ShapeFunctions::getDimension(ElementType elementType) {
    switch (elementType) {
        case ElementType::LINEAR:
            return 1;
        case ElementType::HEXAHEDRON:
        case ElementType::TETRAHEDRON:
            return 3;
        default:
            return 0;
    }
}

// Private implementation methods

ShapeFunctions::ShapeResult ShapeFunctions::compute1DShapeFunctions(const std::vector<Node>& nodes,
                                                                   double xi) {
    int nNodes = nodes.size();
    ShapeResult result(nNodes);
    
    result.values = linear1D(xi, nNodes);
    result.dNdxi = linear1DDerivatives(xi, nNodes);
    
    // For 1D elements, Jacobian is simply the length derivative
    double dxdxi = 0.0;
    for (int i = 0; i < nNodes; ++i) {
        dxdxi += nodes[i].x * result.dNdxi[i];
    }
    result.detJ = std::abs(dxdxi);
    
    // Global derivatives are just scaled by Jacobian
    for (int i = 0; i < nNodes; ++i) {
        result.dNdx[i] = result.dNdxi[i] / result.detJ;
    }
    
    return result;
}

ShapeFunctions::ShapeResult ShapeFunctions::compute2DQuadrilateralShapeFunctions(const std::vector<Node>& nodes,
                                                                               double xi, double eta) {
    int nNodes = nodes.size();
    ShapeResult result(nNodes);
    
    result.values = linearQuadrilateral(xi, eta, nNodes);
    auto [dNdxi, dNdeta] = linearQuadrilateralDerivatives(xi, eta, nNodes);
    result.dNdxi = dNdxi;
    result.dNdeta = dNdeta;
    
    // Compute 2D Jacobian matrix
    std::array<std::array<double, 2>, 2> jac2D = {{{0.0, 0.0}, {0.0, 0.0}}};
    for (int i = 0; i < nNodes; ++i) {
        jac2D[0][0] += nodes[i].x * dNdxi[i];    // dx/dxi
        jac2D[0][1] += nodes[i].y * dNdxi[i];    // dy/dxi
        jac2D[1][0] += nodes[i].x * dNdeta[i];   // dx/deta
        jac2D[1][1] += nodes[i].y * dNdeta[i];   // dy/deta
    }
    
    // Determinant for 2D elements
    result.detJ = jac2D[0][0] * jac2D[1][1] - jac2D[0][1] * jac2D[1][0];
    
    // Compute inverse Jacobian for 2D
    std::array<std::array<double, 2>, 2> invJac2D = {{{0.0, 0.0}, {0.0, 0.0}}};
    if (std::abs(result.detJ) > 1e-15) {
        invJac2D[0][0] = jac2D[1][1] / result.detJ;
        invJac2D[0][1] = -jac2D[0][1] / result.detJ;
        invJac2D[1][0] = -jac2D[1][0] / result.detJ;
        invJac2D[1][1] = jac2D[0][0] / result.detJ;
    } else {
        throw std::runtime_error("Jacobian determinant is zero or too small for 2D quadrilateral");
    }
    
    // Transform derivatives to global coordinates for 2D
    for (int i = 0; i < nNodes; ++i) {
        result.dNdx[i] = invJac2D[0][0] * dNdxi[i] + invJac2D[0][1] * dNdeta[i];
        result.dNdy[i] = invJac2D[1][0] * dNdxi[i] + invJac2D[1][1] * dNdeta[i];
    }
    
    return result;
}

ShapeFunctions::ShapeResult ShapeFunctions::compute3DHexahedronShapeFunctions(const std::vector<Node>& nodes,
                                                                             double xi, double eta, double zeta) {
    int nNodes = nodes.size();
    ShapeResult result(nNodes);
    
    result.values = linearHexahedron(xi, eta, zeta, nNodes);
    auto [dNdxi, dNdeta, dNdzeta] = linearHexahedronDerivatives(xi, eta, zeta, nNodes);
    result.dNdxi = dNdxi;
    result.dNdeta = dNdeta;
    result.dNdzeta = dNdzeta;
    
    // Compute Jacobian matrix
    auto jac = computeJacobianMatrix(nodes, dNdxi, dNdeta, dNdzeta);
    result.detJ = computeJacobianDeterminant(jac);
    
    // Compute inverse Jacobian
    auto invJac = computeInverseJacobian(jac);
    
    // Transform derivatives to global coordinates
    auto [dNdx, dNdy, dNdz] = transformDerivatives(invJac, dNdxi, dNdeta, dNdzeta);
    result.dNdx = dNdx;
    result.dNdy = dNdy;
    result.dNdz = dNdz;
    
    return result;
}

ShapeFunctions::ShapeResult ShapeFunctions::computeTriangleShapeFunctions(const std::vector<Node>& nodes,
                                                                         double xi, double eta) {
    int nNodes = nodes.size();
    ShapeResult result(nNodes);
    
    result.values = linearTriangle(xi, eta, nNodes);
    auto [dNdxi, dNdeta] = linearTriangleDerivatives(xi, eta, nNodes);
    result.dNdxi = dNdxi;
    result.dNdeta = dNdeta;
    
    // Compute Jacobian matrix (2D triangle)
    std::array<std::array<double, 3>, 3> jac = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    for (int i = 0; i < nNodes; ++i) {
        jac[0][0] += nodes[i].x * dNdxi[i];    // dx/dxi
        jac[0][1] += nodes[i].y * dNdxi[i];    // dy/dxi
        jac[1][0] += nodes[i].x * dNdeta[i];   // dx/deta
        jac[1][1] += nodes[i].y * dNdeta[i];   // dy/deta
    }
    jac[2][2] = 1.0;  // For 2D elements, dz/dzeta = 1
    
    // Determinant for 2D elements (area element)
    result.detJ = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
    
    // Inverse Jacobian for 2D
    std::array<std::array<double, 3>, 3> invJac = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    if (std::abs(result.detJ) > 1e-15) {
        invJac[0][0] = jac[1][1] / result.detJ;
        invJac[0][1] = -jac[0][1] / result.detJ;
        invJac[1][0] = -jac[1][0] / result.detJ;
        invJac[1][1] = jac[0][0] / result.detJ;
        invJac[2][2] = 1.0;
    }
    
    // Transform derivatives to global coordinates
    for (int i = 0; i < nNodes; ++i) {
        result.dNdx[i] = invJac[0][0] * dNdxi[i] + invJac[0][1] * dNdeta[i];
        result.dNdy[i] = invJac[1][0] * dNdxi[i] + invJac[1][1] * dNdeta[i];
    }
    
    return result;
}

ShapeFunctions::ShapeResult ShapeFunctions::computeTetrahedronShapeFunctions(const std::vector<Node>& nodes,
                                                                             double xi, double eta, double zeta) {
    int nNodes = nodes.size();
    ShapeResult result(nNodes);
    
    result.values = linearTetrahedron(xi, eta, zeta, nNodes);
    auto [dNdxi, dNdeta, dNdzeta] = linearTetrahedronDerivatives(xi, eta, zeta, nNodes);
    result.dNdxi = dNdxi;
    result.dNdeta = dNdeta;
    result.dNdzeta = dNdzeta;
    
    // Compute Jacobian matrix
    auto jac = computeJacobianMatrix(nodes, dNdxi, dNdeta, dNdzeta);
    result.detJ = computeJacobianDeterminant(jac);
    
    // Compute inverse Jacobian
    auto invJac = computeInverseJacobian(jac);
    
    // Transform derivatives to global coordinates
    auto [dNdx, dNdy, dNdz] = transformDerivatives(invJac, dNdxi, dNdeta, dNdzeta);
    result.dNdx = dNdx;
    result.dNdy = dNdy;
    result.dNdz = dNdz;
    
    return result;
}

ShapeFunctions::ShapeResult ShapeFunctions::computePrismShapeFunctions(const std::vector<Node>& nodes,
                                                                      double xi, double eta, double zeta) {
    int nNodes = nodes.size();
    ShapeResult result(nNodes);
    
    // Use wedge nodal basis functions for prism elements
    if (nNodes == 6) {
        result.values = wedgeNodalPBasisAll(xi, eta, zeta);
        auto gradients = dWedgeNodalPBasisAll(xi, eta, zeta);
        
        // Extract derivatives from gradient matrix
        for (int i = 0; i < nNodes; ++i) {
            result.dNdxi[i] = gradients[i][0];
            result.dNdeta[i] = gradients[i][1];
            result.dNdzeta[i] = gradients[i][2];
        }
    } else {
        throw std::invalid_argument("Unsupported number of nodes for prism element");
    }
    
    // Compute Jacobian matrix
    auto jac = computeJacobianMatrix(nodes, result.dNdxi, result.dNdeta, result.dNdzeta);
    result.detJ = computeJacobianDeterminant(jac);
    
    // Compute inverse Jacobian
    auto invJac = computeInverseJacobian(jac);
    
    // Transform derivatives to global coordinates
    auto [dNdx, dNdy, dNdz] = transformDerivatives(invJac, result.dNdxi, result.dNdeta, result.dNdzeta);
    result.dNdx = dNdx;
    result.dNdy = dNdy;
    result.dNdz = dNdz;
    
    return result;
}

// Wedge/Prism shape function implementations

double ShapeFunctions::wedgeNodalPBasis(int node, double u, double v, double w) {
    double value = 0.0;
    
    switch (node) {
        case 1:
            value = wedgeL(1, u, v) * (1 - w);
            break;
        case 2:
            value = wedgeL(2, u, v) * (1 - w);
            break;
        case 3:
            value = wedgeL(3, u, v) * (1 - w);
            break;
        case 4:
            value = wedgeL(1, u, v) * (1 + w);
            break;
        case 5:
            value = wedgeL(2, u, v) * (1 + w);
            break;
        case 6:
            value = wedgeL(3, u, v) * (1 + w);
            break;
        default:
            throw std::invalid_argument("Unknown node for wedge element");
    }
    
    return value / 2.0;
}

std::vector<double> ShapeFunctions::wedgeNodalPBasisAll(double u, double v, double w) {
    std::vector<double> phi(6, 0.0);
    std::vector<double> tri(3, 0.0);
    std::vector<double> line(2, 0.0);
    
    const double half = 0.5;
    const double c3 = 1.0 / std::sqrt(3.0);
    
    tri[0] = half * (1.0 - u - c3 * v);
    tri[1] = half * (1.0 + u - c3 * v);
    tri[2] = c3 * v;
    
    line[0] = half * (1 - w);
    line[1] = half * (1 + w);
    
    phi[0] = line[0] * tri[0];
    phi[1] = line[0] * tri[1];
    phi[2] = line[0] * tri[2];
    phi[3] = line[1] * tri[0];
    phi[4] = line[1] * tri[1];
    phi[5] = line[1] * tri[2];
    
    return phi;
}

std::vector<double> ShapeFunctions::wedgeNodalLBasisAll(double u, double v, double w) {
    std::vector<double> phi(6, 0.0);
    std::vector<double> tri(3, 0.0);
    std::vector<double> line(2, 0.0);
    
    const double half = 0.5;
    
    tri[0] = 1.0 - u - v;
    tri[1] = u;
    tri[2] = v;
    
    line[0] = half * (1 - w);
    line[1] = half * (1 + w);
    
    phi[0] = line[0] * tri[0];
    phi[1] = line[0] * tri[1];
    phi[2] = line[0] * tri[2];
    phi[3] = line[1] * tri[0];
    phi[4] = line[1] * tri[1];
    phi[5] = line[1] * tri[2];
    
    return phi;
}

std::vector<double> ShapeFunctions::dWedgeNodalPBasis(int node, double u, double v, double w) {
    std::vector<double> grad(3, 0.0);
    double signW = 0.0;
    
    switch (node) {
        case 1: case 2: case 3:
            signW = -1.0;
            break;
        case 4: case 5: case 6:
            signW = 1.0;
            break;
        default:
            throw std::invalid_argument("Unknown node for wedge element");
    }
    
    // Calculate gradient from the general form
    auto dL = dWedgeL(node, u, v);
    double L = wedgeL(node, u, v);
    
    grad[0] = 0.5 * dL[0] * (1 + signW * w);
    grad[1] = 0.5 * dL[1] * (1 + signW * w);
    grad[2] = signW * 0.5 * L;
    
    return grad;
}

std::vector<std::vector<double>> ShapeFunctions::dWedgeNodalPBasisAll(double u, double v, double w) {
    std::vector<std::vector<double>> gradPhi(6, std::vector<double>(3, 0.0));
    std::vector<double> tri(3, 0.0);
    std::vector<double> line(2, 0.0);
    std::vector<std::vector<double>> gradTri(3, std::vector<double>(2, 0.0));
    std::vector<double> gradLine(2, 0.0);
    
    const double half = 0.5;
    const double c3 = 1.0 / std::sqrt(3.0);
    
    tri[0] = half * (1.0 - u - c3 * v);
    tri[1] = half * (1.0 + u - c3 * v);
    tri[2] = c3 * v;
    
    line[0] = half * (1 - w);
    line[1] = half * (1 + w);
    
    gradTri[0][0] = -half;
    gradTri[0][1] = -half * c3;
    gradTri[1][0] = half;
    gradTri[1][1] = -half * c3;
    gradTri[2][0] = 0.0;
    gradTri[2][1] = c3;
    
    gradLine[0] = -half;
    gradLine[1] = half;
    
    // Calculate gradients for nodes 1-3
    for (int i = 0; i < 3; ++i) {
        gradPhi[i][0] = gradTri[i][0] * line[0];
        gradPhi[i][1] = gradTri[i][1] * line[0];
        gradPhi[i][2] = tri[i] * gradLine[0];
    }
    
    // Calculate gradients for nodes 4-6
    for (int i = 0; i < 3; ++i) {
        gradPhi[i+3][0] = gradTri[i][0] * line[1];
        gradPhi[i+3][1] = gradTri[i][1] * line[1];
        gradPhi[i+3][2] = tri[i] * gradLine[1];
    }
    
    return gradPhi;
}

std::vector<std::vector<double>> ShapeFunctions::dWedgeNodalLBasisAll(double u, double v, double w) {
    std::vector<std::vector<double>> gradPhi(6, std::vector<double>(3, 0.0));
    std::vector<double> tri(3, 0.0);
    std::vector<double> line(2, 0.0);
    std::vector<std::vector<double>> gradTri(3, std::vector<double>(2, 0.0));
    std::vector<double> gradLine(2, 0.0);
    
    const double half = 0.5;
    
    tri[0] = 1.0 - u - v;
    tri[1] = u;
    tri[2] = v;
    
    line[0] = half * (1 - w);
    line[1] = half * (1 + w);
    
    gradTri[0][0] = -1.0;
    gradTri[0][1] = -1.0;
    gradTri[1][0] = 1.0;
    gradTri[1][1] = 0.0;
    gradTri[2][0] = 0.0;
    gradTri[2][1] = 1.0;
    
    gradLine[0] = -half;
    gradLine[1] = half;
    
    // Calculate gradients for nodes 1-3
    for (int i = 0; i < 3; ++i) {
        gradPhi[i][0] = gradTri[i][0] * line[0];
        gradPhi[i][1] = gradTri[i][1] * line[0];
        gradPhi[i][2] = tri[i] * gradLine[0];
    }
    
    // Calculate gradients for nodes 4-6
    for (int i = 0; i < 3; ++i) {
        gradPhi[i+3][0] = gradTri[i][0] * line[1];
        gradPhi[i+3][1] = gradTri[i][1] * line[1];
        gradPhi[i+3][2] = tri[i] * gradLine[1];
    }
    
    return gradPhi;
}

std::vector<std::vector<double>> ShapeFunctions::ddWedgeNodalPBasis(int node, double u, double v, double w) {
    std::vector<std::vector<double>> grad(3, std::vector<double>(3, 0.0));
    double signW = 0.0;
    
    switch (node) {
        case 1: case 2: case 3:
            signW = -1.0;
            break;
        case 4: case 5: case 6:
            signW = 1.0;
            break;
        default:
            throw std::invalid_argument("Unknown node for wedge element");
    }
    
    // Calculate second derivatives
    auto dL = dWedgeL(node, u, v);
    
    grad[0][2] = dL[0] * signW / 2.0;
    grad[1][2] = dL[1] * signW / 2.0;
    grad[2][0] = grad[0][2];
    grad[2][1] = grad[1][2];
    
    return grad;
}

// Helper functions for wedge shape functions

double ShapeFunctions::wedgeL(int node, double u, double v) {
    switch (node) {
        case 1:
            return 1.0 - u - v;
        case 2:
            return u;
        case 3:
            return v;
        default:
            throw std::invalid_argument("Unknown node for wedge linear basis");
    }
}

std::vector<double> ShapeFunctions::dWedgeL(int node, double u, double v) {
    std::vector<double> grad(2, 0.0);
    
    switch (node) {
        case 1:
            grad[0] = -1.0;
            grad[1] = -1.0;
            break;
        case 2:
            grad[0] = 1.0;
            grad[1] = 0.0;
            break;
        case 3:
            grad[0] = 0.0;
            grad[1] = 1.0;
            break;
        default:
            throw std::invalid_argument("Unknown node for wedge linear basis");
    }
    
    return grad;
}
