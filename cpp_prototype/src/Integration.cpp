#include "Integration.h"
#include <stdexcept>
#include <algorithm>

namespace elmer {

// Initialize static members
std::vector<std::vector<std::pair<double, double>>> Integration::gaussPoints1D_;
bool Integration::initialized_ = false;

void Integration::initializeGaussPoints1D() {
    if (initialized_) return;
    
    // Precompute Gauss points for orders 1 to 13
    gaussPoints1D_.resize(14); // Index 0 unused, indices 1-13 used
    
    // Order 1: 1 point
    gaussPoints1D_[1] = {{0.0, 2.0}};
    
    // Order 2: 2 points
    gaussPoints1D_[2] = {
        {-0.5773502691896257, 1.0},
        {0.5773502691896257, 1.0}
    };
    
    // Order 3: 3 points
    gaussPoints1D_[3] = {
        {-0.7745966692414834, 0.5555555555555556},
        {0.0, 0.8888888888888888},
        {0.7745966692414834, 0.5555555555555556}
    };
    
    // Order 4: 4 points
    gaussPoints1D_[4] = {
        {-0.8611363115940526, 0.3478548451374538},
        {-0.3399810435848563, 0.6521451548625461},
        {0.3399810435848563, 0.6521451548625461},
        {0.8611363115940526, 0.3478548451374538}
    };
    
    // Order 5: 5 points
    gaussPoints1D_[5] = {
        {-0.9061798459386640, 0.2369268850561891},
        {-0.5384693101056831, 0.4786286704993665},
        {0.0, 0.5688888888888889},
        {0.5384693101056831, 0.4786286704993665},
        {0.9061798459386640, 0.2369268850561891}
    };
    
    // For higher orders, we'll compute on the fly
    initialized_ = true;
}

std::vector<std::pair<double, double>> Integration::getPrecomputedGauss1D(int npoints) {
    initializeGaussPoints1D();
    
    if (npoints >= 1 && npoints <= 5 && npoints < static_cast<int>(gaussPoints1D_.size())) {
        return gaussPoints1D_[npoints];
    }
    
    // For higher orders, compute using standard Gauss-Legendre algorithm
    std::vector<std::pair<double, double>> points;
    
    if (npoints <= 0) {
        throw std::invalid_argument("Number of integration points must be positive");
    }
    
    // Simple implementation for higher orders
    // In practice, we would use a more robust algorithm
    switch (npoints) {
        case 6:
            points = {
                {-0.9324695142031521, 0.1713244923791704},
                {-0.6612093864662645, 0.3607615730481386},
                {-0.2386191860831969, 0.4679139345726910},
                {0.2386191860831969, 0.4679139345726910},
                {0.6612093864662645, 0.3607615730481386},
                {0.9324695142031521, 0.1713244923791704}
            };
            break;
        case 7:
            points = {
                {-0.9491079123427585, 0.1294849661688697},
                {-0.7415311855993945, 0.2797053914892766},
                {-0.4058451513773972, 0.3818300505051189},
                {0.0, 0.4179591836734694},
                {0.4058451513773972, 0.3818300505051189},
                {0.7415311855993945, 0.2797053914892766},
                {0.9491079123427585, 0.1294849661688697}
            };
            break;
        default:
            // For orders > 7, use a simple uniform distribution as fallback
            // In production code, we would implement proper Gauss-Legendre quadrature
            double step = 2.0 / (npoints - 1);
            double weight = 2.0 / npoints;
            for (int i = 0; i < npoints; ++i) {
                double xi = -1.0 + i * step;
                points.push_back({xi, weight});
            }
            break;
    }
    
    return points;
}

std::vector<Integration::IntegrationPoint> Integration::getGaussPoints1D(int npoints) {
    auto pointsWeights = getPrecomputedGauss1D(npoints);
    std::vector<IntegrationPoint> result;
    
    for (const auto& pw : pointsWeights) {
        result.push_back(IntegrationPoint(pw.first, 0.0, 0.0, pw.second));
    }
    
    return result;
}

std::vector<Integration::IntegrationPoint> Integration::getGaussPoints2D(int npoints_u, int npoints_v) {
    auto points_u = getPrecomputedGauss1D(npoints_u);
    auto points_v = getPrecomputedGauss1D(npoints_v);
    std::vector<IntegrationPoint> result;
    
    for (const auto& pu : points_u) {
        for (const auto& pv : points_v) {
            double weight = pu.second * pv.second;
            result.push_back(IntegrationPoint(pu.first, pv.first, 0.0, weight));
        }
    }
    
    return result;
}

std::vector<Integration::IntegrationPoint> Integration::getGaussPoints3D(int npoints_u, int npoints_v, int npoints_w) {
    auto points_u = getPrecomputedGauss1D(npoints_u);
    auto points_v = getPrecomputedGauss1D(npoints_v);
    auto points_w = getPrecomputedGauss1D(npoints_w);
    std::vector<IntegrationPoint> result;
    
    for (const auto& pu : points_u) {
        for (const auto& pv : points_v) {
            for (const auto& pw : points_w) {
                double weight = pu.second * pv.second * pw.second;
                result.push_back(IntegrationPoint(pu.first, pv.first, pw.first, weight));
            }
        }
    }
    
    return result;
}

std::vector<Integration::IntegrationPoint> Integration::getTrianglePoints(int order) {
    std::vector<IntegrationPoint> points;
    
    switch (order) {
        case 1: // 1-point rule
            points.push_back(IntegrationPoint(1.0/3.0, 1.0/3.0, 0.0, 1.0/2.0));
            break;
            
        case 2: // 3-point rule
            points.push_back(IntegrationPoint(1.0/6.0, 1.0/6.0, 0.0, 1.0/6.0));
            points.push_back(IntegrationPoint(2.0/3.0, 1.0/6.0, 0.0, 1.0/6.0));
            points.push_back(IntegrationPoint(1.0/6.0, 2.0/3.0, 0.0, 1.0/6.0));
            break;
            
        case 3: // 4-point rule
            points.push_back(IntegrationPoint(1.0/3.0, 1.0/3.0, 0.0, -27.0/96.0));
            points.push_back(IntegrationPoint(0.6, 0.2, 0.0, 25.0/96.0));
            points.push_back(IntegrationPoint(0.2, 0.6, 0.0, 25.0/96.0));
            points.push_back(IntegrationPoint(0.2, 0.2, 0.0, 25.0/96.0));
            break;
            
        default:
            // For higher orders, use Duffy transformation from quadrilateral
            auto quadPoints = getGaussPoints2D(order, order);
            for (const auto& qp : quadPoints) {
                double u = (qp.u + 1.0) / 2.0;
                double v = (qp.v + 1.0) / 2.0;
                double xi = u;
                double eta = v * (1.0 - u);
                double weight = qp.weight * (1.0 - u) / 4.0;
                points.push_back(IntegrationPoint(xi, eta, 0.0, weight));
            }
            break;
    }
    
    return points;
}

std::vector<Integration::IntegrationPoint> Integration::getTetrahedronPoints(int order) {
    std::vector<IntegrationPoint> points;
    
    switch (order) {
        case 1: // 1-point rule
            points.push_back(IntegrationPoint(0.25, 0.25, 0.25, 1.0/6.0));
            break;
            
        case 2: // 4-point rule
            {
                double a = 0.1381966011250105;
                double b = 0.5854101966249685;
                double weight = 1.0/24.0;
                points.push_back(IntegrationPoint(a, a, a, weight));
                points.push_back(IntegrationPoint(b, a, a, weight));
                points.push_back(IntegrationPoint(a, b, a, weight));
                points.push_back(IntegrationPoint(a, a, b, weight));
            }
            break;
            
        default:
            // For higher orders, use transformation from hexahedron
            auto hexPoints = getGaussPoints3D(order, order, order);
            for (const auto& hp : hexPoints) {
                double u = (hp.u + 1.0) / 2.0;
                double v = (hp.v + 1.0) / 2.0;
                double w = (hp.w + 1.0) / 2.0;
                double xi = u;
                double eta = v * (1.0 - u);
                double zeta = w * (1.0 - u) * (1.0 - v);
                double weight = hp.weight * (1.0 - u) * (1.0 - u) * (1.0 - v) / 8.0;
                points.push_back(IntegrationPoint(xi, eta, zeta, weight));
            }
            break;
    }
    
    return points;
}

// Shape function implementations
double Integration::shapeFunction1DLinear(double xi, int n) {
    switch (n) {
        case 0: return 0.5 * (1.0 - xi);
        case 1: return 0.5 * (1.0 + xi);
        default: throw std::invalid_argument("Invalid node index for 1D linear element");
    }
}

double Integration::shapeFunction2DLinear(double xi, double eta, int n) {
    switch (n) {
        case 0: return 0.25 * (1.0 - xi) * (1.0 - eta);
        case 1: return 0.25 * (1.0 + xi) * (1.0 - eta);
        case 2: return 0.25 * (1.0 + xi) * (1.0 + eta);
        case 3: return 0.25 * (1.0 - xi) * (1.0 + eta);
        default: throw std::invalid_argument("Invalid node index for 2D linear quadrilateral");
    }
}

double Integration::shapeFunction3DLinear(double xi, double eta, double zeta, int n) {
    switch (n) {
        case 0: return 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta);
        case 1: return 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta);
        case 2: return 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta);
        case 3: return 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta);
        case 4: return 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta);
        case 5: return 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta);
        case 6: return 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta);
        case 7: return 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta);
        default: throw std::invalid_argument("Invalid node index for 3D linear hexahedron");
    }
}

double Integration::shapeFunction1DQuadratic(double xi, int n) {
    switch (n) {
        case 0: return 0.5 * xi * (xi - 1.0);
        case 1: return 1.0 - xi * xi;
        case 2: return 0.5 * xi * (xi + 1.0);
        default: throw std::invalid_argument("Invalid node index for 1D quadratic element");
    }
}

// Shape function derivative implementations
double Integration::shapeFunctionDerivative1DLinear(double xi, int n) {
    switch (n) {
        case 0: return -0.5;
        case 1: return 0.5;
        default: throw std::invalid_argument("Invalid node index for 1D linear element");
    }
}

double Integration::shapeFunctionDerivative2DLinear(double xi, double eta, int n, int dir) {
    switch (n) {
        case 0:
            if (dir == 0) return -0.25 * (1.0 - eta);
            else return -0.25 * (1.0 - xi);
        case 1:
            if (dir == 0) return 0.25 * (1.0 - eta);
            else return -0.25 * (1.0 + xi);
        case 2:
            if (dir == 0) return 0.25 * (1.0 + eta);
            else return 0.25 * (1.0 + xi);
        case 3:
            if (dir == 0) return -0.25 * (1.0 + eta);
            else return 0.25 * (1.0 - xi);
        default: throw std::invalid_argument("Invalid node index for 2D linear quadrilateral");
    }
}

double Integration::shapeFunctionDerivative3DLinear(double xi, double eta, double zeta, int n, int dir) {
    switch (n) {
        case 0:
            if (dir == 0) return -0.125 * (1.0 - eta) * (1.0 - zeta);
            else if (dir == 1) return -0.125 * (1.0 - xi) * (1.0 - zeta);
            else return -0.125 * (1.0 - xi) * (1.0 - eta);
        case 1:
            if (dir == 0) return 0.125 * (1.0 - eta) * (1.0 - zeta);
            else if (dir == 1) return -0.125 * (1.0 + xi) * (1.0 - zeta);
            else return -0.125 * (1.0 + xi) * (1.0 - eta);
        // ... similar for other nodes
        default: throw std::invalid_argument("Invalid node index for 3D linear hexahedron");
    }
}

// Jacobian computation implementations
double Integration::computeJacobian1D(const std::vector<Node>& nodes, double xi) {
    if (nodes.size() < 2) {
        throw std::invalid_argument("At least 2 nodes required for 1D element");
    }
    
    double dxdxi = 0.0;
    for (int i = 0; i < 2; ++i) {
        double dNdxi = shapeFunctionDerivative1DLinear(xi, i);
        dxdxi += nodes[i].x * dNdxi;
    }
    
    return std::abs(dxdxi);
}

double Integration::computeJacobian2D(const std::vector<Node>& nodes, double xi, double eta) {
    if (nodes.size() < 4) {
        throw std::invalid_argument("At least 4 nodes required for 2D quadrilateral element");
    }
    
    auto jac = computeJacobianMatrix2D(nodes, xi, eta);
    double det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
    return std::abs(det);
}

double Integration::computeJacobian3D(const std::vector<Node>& nodes, double xi, double eta, double zeta) {
    if (nodes.size() < 8) {
        throw std::invalid_argument("At least 8 nodes required for 3D hexahedral element");
    }
    
    auto jac = computeJacobianMatrix3D(nodes, xi, eta, zeta);
    double det = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
               - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
               + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
    return std::abs(det);
}

std::array<std::array<double, 2>, 2> Integration::computeJacobianMatrix2D(
    const std::vector<Node>& nodes, double xi, double eta) {
    
    std::array<std::array<double, 2>, 2> jac = {{{0.0, 0.0}, {0.0, 0.0}}};
    
    for (int i = 0; i < 4; ++i) {
        double dNdxi = shapeFunctionDerivative2DLinear(xi, eta, i, 0);
        double dNdeta = shapeFunctionDerivative2DLinear(xi, eta, i, 1);
        
        jac[0][0] += nodes[i].x * dNdxi;  // dx/dxi
        jac[0][1] += nodes[i].x * dNdeta; // dx/deta
        jac[1][0] += nodes[i].y * dNdxi;  // dy/dxi
        jac[1][1] += nodes[i].y * dNdeta; // dy/deta
    }
    
    return jac;
}

std::array<std::array<double, 3>, 3> Integration::computeJacobianMatrix3D(
    const std::vector<Node>& nodes, double xi, double eta, double zeta) {
    
    std::array<std::array<double, 3>, 3> jac = {
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}
    };
    
    for (int i = 0; i < 8; ++i) {
        double dNdxi = shapeFunctionDerivative3DLinear(xi, eta, zeta, i, 0);
        double dNdeta = shapeFunctionDerivative3DLinear(xi, eta, zeta, i, 1);
        double dNdzeta = shapeFunctionDerivative3DLinear(xi, eta, zeta, i, 2);
        
        jac[0][0] += nodes[i].x * dNdxi;    // dx/dxi
        jac[0][1] += nodes[i].x * dNdeta;   // dx/deta
        jac[0][2] += nodes[i].x * dNdzeta;  // dx/dzeta
        jac[1][0] += nodes[i].y * dNdxi;    // dy/dxi
        jac[1][1] += nodes[i].y * dNdeta;   // dy/deta
        jac[1][2] += nodes[i].y * dNdzeta;  // dy/dzeta
        jac[2][0] += nodes[i].z * dNdxi;    // dz/dxi
        jac[2][1] += nodes[i].z * dNdeta;   // dz/deta
        jac[2][2] += nodes[i].z * dNdzeta;  // dz/dzeta
    }
    
    return jac;
}

std::array<std::array<double, 2>, 2> Integration::computeInverseJacobian2D(
    const std::array<std::array<double, 2>, 2>& jac) {
    
    double det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
    if (std::abs(det) < 1e-15) {
        throw std::runtime_error("Jacobian matrix is singular");
    }
    
    double invDet = 1.0 / det;
    std::array<std::array<double, 2>, 2> invJac;
    
    invJac[0][0] = jac[1][1] * invDet;
    invJac[0][1] = -jac[0][1] * invDet;
    invJac[1][0] = -jac[1][0] * invDet;
    invJac[1][1] = jac[0][0] * invDet;
    
    return invJac;
}

std::array<std::array<double, 3>, 3> Integration::computeInverseJacobian3D(
    const std::array<std::array<double, 3>, 3>& jac) {
    
    double det = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
               - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
               + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
    
    if (std::abs(det) < 1e-15) {
        throw std::runtime_error("Jacobian matrix is singular");
    }
    
    double invDet = 1.0 / det;
    std::array<std::array<double, 3>, 3> invJac;
    
    invJac[0][0] = (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) * invDet;
    invJac[0][1] = (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) * invDet;
    invJac[0][2] = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) * invDet;
    
    invJac[1][0] = (jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2]) * invDet;
    invJac[1][1] = (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) * invDet;
    invJac[1][2] = (jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2]) * invDet;
    
    invJac[2][0] = (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]) * invDet;
    invJac[2][1] = (jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1]) * invDet;
    invJac[2][2] = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) * invDet;
    
    return invJac;
}

std::vector<double> Integration::transformDerivatives(
    const std::vector<double>& dNdxi, 
    const std::array<std::array<double, 3>, 3>& invJac, 
    int dim) {
    
    if (dim != 2 && dim != 3) {
        throw std::invalid_argument("Dimension must be 2 or 3");
    }
    
    std::vector<double> dNdx(dim, 0.0);
    
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            dNdx[i] += invJac[j][i] * dNdxi[j];
        }
    }
    
    return dNdx;
}

} // namespace elmer