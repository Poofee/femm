// GaussIntegration.cpp - Elmer FEM C++ Gaussian Integration Module Implementation
// Corresponds to Fortran module: Integration.F90

#include "GaussIntegration.h"
#include <cmath>
#include <stdexcept>

using namespace elmer;

// 1D Gaussian quadrature data (points and weights)
const std::array<std::vector<double>, 6> GaussIntegration::gauss1DData = {{
    {0.0, 2.0},                           // 1 point
    {-0.5773502691896257, 1.0, 0.5773502691896257, 1.0}, // 2 points
    {-0.7745966692414834, 0.5555555555555556, 0.0, 0.8888888888888888, 0.7745966692414834, 0.5555555555555556}, // 3 points
    {-0.8611363115940526, 0.3478548451374538, -0.3399810435848563, 0.6521451548625461, 
     0.3399810435848563, 0.6521451548625461, 0.8611363115940526, 0.3478548451374538}, // 4 points
    {-0.9061798459386640, 0.2369268850561891, -0.5384693101056831, 0.4786286704993665, 
     0.0, 0.5688888888888889, 0.5384693101056831, 0.4786286704993665, 
     0.9061798459386640, 0.2369268850561891}, // 5 points
    {-0.9324695142031521, 0.1713244923791704, -0.6612093864662645, 0.3607615730481386, 
     -0.2386191860831969, 0.4679139345726910, 0.2386191860831969, 0.4679139345726910, 
     0.6612093864662645, 0.3607615730481386, 0.9324695142031521, 0.1713244923791704} // 6 points
}};

// Triangular integration data (xi, eta, weight)
const std::array<std::vector<double>, 7> GaussIntegration::gaussTriangularData = {{
    {1.0/3.0, 1.0/3.0, 0.0, 0.5}, // 1 point, order 1
    {0.1666666666666667, 0.1666666666666667, 0.0, 0.1666666666666667, 
     0.6666666666666667, 0.1666666666666667, 0.0, 0.1666666666666667, 
     0.1666666666666667, 0.6666666666666667, 0.0, 0.1666666666666667}, // 3 points, order 2
    {0.3333333333333333, 0.3333333333333333, 0.0, -0.28125, 
     0.2, 0.2, 0.0, 0.2604166666666667, 
     0.6, 0.2, 0.0, 0.2604166666666667, 
     0.2, 0.6, 0.0, 0.2604166666666667}, // 4 points, order 3
    {0.4459484909159650, 0.4459484909159650, 0.0, 0.1116907948390057, 
     0.0915762135097707, 0.0915762135097707, 0.0, 0.0549758718276609, 
     0.1081030181680700, 0.4459484909159650, 0.0, 0.1116907948390057, 
     0.8168475729804585, 0.0915762135097707, 0.0, 0.0549758718276609, 
     0.4459484909159650, 0.1081030181680700, 0.0, 0.1116907948390057, 
     0.0915762135097707, 0.8168475729804585, 0.0, 0.0549758718276609}, // 6 points, order 4
    {0.3333333333333333, 0.3333333333333333, 0.0, 0.1125, 
     0.4701420641051151, 0.4701420641051151, 0.0, 0.0661970763942531, 
     0.0597158717897698, 0.4701420641051151, 0.0, 0.0661970763942531, 
     0.4701420641051151, 0.0597158717897698, 0.0, 0.0661970763942531, 
     0.1012865073234563, 0.1012865073234563, 0.0, 0.0629695902724136, 
     0.7974269853530873, 0.1012865073234563, 0.0, 0.0629695902724136, 
     0.1012865073234563, 0.7974269853530873, 0.0, 0.0629695902724136} // 7 points, order 5
}};

// Tetrahedral integration data (xi, eta, zeta, weight)
const std::array<std::vector<double>, 5> GaussIntegration::gaussTetrahedralData = {{
    {0.25, 0.25, 0.25, 0.25, 1.0/6.0}, // 1 point, order 1
    {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.0416666666666667, 
     0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.0416666666666667, 
     0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.0416666666666667, 
     0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.0416666666666667}, // 4 points, order 2
    {0.25, 0.25, 0.25, 0.25, -0.1333333333333333, 
     0.5, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.075, 
     0.1666666666666667, 0.5, 0.1666666666666667, 0.1666666666666667, 0.075, 
     0.1666666666666667, 0.1666666666666667, 0.5, 0.1666666666666667, 0.075, 
     0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.5, 0.075}, // 5 points, order 3
    {0.25, 0.25, 0.25, 0.25, 0.0197530864197531, 
     0.7857142857142857, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0119895139631698, 
     0.0714285714285714, 0.7857142857142857, 0.0714285714285714, 0.0714285714285714, 0.0119895139631698, 
     0.0714285714285714, 0.0714285714285714, 0.7857142857142857, 0.0714285714285714, 0.0119895139631698, 
     0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857, 0.0119895139631698, 
     0.3994035761667992, 0.3994035761667992, 0.1005964238332008, 0.1005964238332008, 0.0115113678710454, 
     0.3994035761667992, 0.1005964238332008, 0.3994035761667992, 0.1005964238332008, 0.0115113678710454, 
     0.3994035761667992, 0.1005964238332008, 0.1005964238332008, 0.3994035761667992, 0.0115113678710454, 
     0.1005964238332008, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008, 0.0115113678710454, 
     0.1005964238332008, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992, 0.0115113678710454, 
     0.1005964238332008, 0.1005964238332008, 0.3994035761667992, 0.3994035761667992, 0.0115113678710454} // 11 points, order 4
}};

std::vector<GaussIntegration::IntegrationPoint> GaussIntegration::get1DPoints(int nPoints) {
    if (nPoints < 1 || nPoints > 6) {
        throw std::invalid_argument("Number of 1D integration points must be between 1 and 6");
    }
    
    std::vector<IntegrationPoint> points;
    const auto& data = gauss1DData[nPoints - 1];
    
    for (size_t i = 0; i < nPoints; ++i) {
        points.emplace_back(data[i * 2], 0.0, 0.0, data[i * 2 + 1]);
    }
    
    return points;
}

std::vector<GaussIntegration::IntegrationPoint> GaussIntegration::getQuadrilateralPoints(int nPoints) {
    if (nPoints < 1 || nPoints > 6) {
        throw std::invalid_argument("Number of quadrilateral integration points must be between 1 and 6");
    }
    
    std::vector<IntegrationPoint> points;
    const auto& xiData = gauss1DData[nPoints - 1];
    const auto& etaData = gauss1DData[nPoints - 1];
    
    for (size_t i = 0; i < nPoints; ++i) {
        for (size_t j = 0; j < nPoints; ++j) {
            double weight = xiData[i * 2 + 1] * etaData[j * 2 + 1];
            points.emplace_back(xiData[i * 2], etaData[j * 2], 0.0, weight);
        }
    }
    
    return points;
}

std::vector<GaussIntegration::IntegrationPoint> GaussIntegration::getHexahedralPoints(int nPoints) {
    if (nPoints < 1 || nPoints > 6) {
        throw std::invalid_argument("Number of hexahedral integration points must be between 1 and 6");
    }
    
    std::vector<IntegrationPoint> points;
    const auto& xiData = gauss1DData[nPoints - 1];
    const auto& etaData = gauss1DData[nPoints - 1];
    const auto& zetaData = gauss1DData[nPoints - 1];
    
    for (size_t i = 0; i < nPoints; ++i) {
        for (size_t j = 0; j < nPoints; ++j) {
            for (size_t k = 0; k < nPoints; ++k) {
                double weight = xiData[i * 2 + 1] * etaData[j * 2 + 1] * zetaData[k * 2 + 1];
                points.emplace_back(xiData[i * 2], etaData[j * 2], zetaData[k * 2], weight);
            }
        }
    }
    
    return points;
}

std::vector<GaussIntegration::IntegrationPoint> GaussIntegration::getTriangularPoints(int order) {
    if (order < 1 || order > 5) {
        throw std::invalid_argument("Triangular integration order must be between 1 and 5");
    }
    
    std::vector<IntegrationPoint> points;
    
    // Map order to number of points
    int nPoints = 0;
    switch (order) {
        case 1: nPoints = 1; break;
        case 2: nPoints = 3; break;
        case 3: nPoints = 4; break;
        case 4: nPoints = 6; break;
        case 5: nPoints = 7; break;
    }
    
    const auto& data = gaussTriangularData[order - 1];
    
    for (int i = 0; i < nPoints; ++i) {
        points.emplace_back(data[i * 4], data[i * 4 + 1], 0.0, data[i * 4 + 3]);
    }
    
    return points;
}

std::vector<GaussIntegration::IntegrationPoint> GaussIntegration::getTetrahedralPoints(int order) {
    if (order < 1 || order > 4) {
        throw std::invalid_argument("Tetrahedral integration order must be between 1 and 4");
    }
    
    std::vector<IntegrationPoint> points;
    
    // Map order to number of points
    int nPoints = 0;
    switch (order) {
        case 1: nPoints = 1; break;
        case 2: nPoints = 4; break;
        case 3: nPoints = 5; break;
        case 4: nPoints = 11; break;
    }
    
    const auto& data = gaussTetrahedralData[order - 1];
    
    for (int i = 0; i < nPoints; ++i) {
        points.emplace_back(data[i * 5], data[i * 5 + 1], data[i * 5 + 2], data[i * 5 + 4]);
    }
    
    return points;
}