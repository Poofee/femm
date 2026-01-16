// GaussIntegration.h - Elmer FEM C++ Gaussian Integration Module
// Corresponds to Fortran module: Integration.F90

#pragma once

#include "Types.h"
#include <vector>
#include <array>

namespace elmer {

/**
 * @brief Gaussian quadrature integration system
 * 
 * This module provides Gaussian quadrature points and weights
 * for numerical integration in finite element analysis.
 */
class GaussIntegration {
public:
    /**
     * @brief Integration point with coordinates and weight
     */
    struct IntegrationPoint {
        double xi;      ///< Natural coordinate in xi direction
        double eta;     ///< Natural coordinate in eta direction
        double zeta;    ///< Natural coordinate in zeta direction
        double weight;  ///< Integration weight
        
        IntegrationPoint(double xi_val, double eta_val, double zeta_val, double w)
            : xi(xi_val), eta(eta_val), zeta(zeta_val), weight(w) {}
    };

    /**
     * @brief Get Gaussian quadrature points for 1D integration
     * @param nPoints Number of integration points (1-6)
     * @return Vector of integration points
     */
    static std::vector<IntegrationPoint> get1DPoints(int nPoints);
    
    /**
     * @brief Get Gaussian quadrature points for 2D quadrilateral integration
     * @param nPoints Number of integration points per dimension (1-6)
     * @return Vector of integration points
     */
    static std::vector<IntegrationPoint> getQuadrilateralPoints(int nPoints);
    
    /**
     * @brief Get Gaussian quadrature points for 3D hexahedral integration
     * @param nPoints Number of integration points per dimension (1-6)
     * @return Vector of integration points
     */
    static std::vector<IntegrationPoint> getHexahedralPoints(int nPoints);
    
    /**
     * @brief Get Gaussian quadrature points for triangular integration
     * @param order Integration order (1-5)
     * @return Vector of integration points
     */
    static std::vector<IntegrationPoint> getTriangularPoints(int order);
    
    /**
     * @brief Get Gaussian quadrature points for tetrahedral integration
     * @param order Integration order (1-5)
     * @return Vector of integration points
     */
    static std::vector<IntegrationPoint> getTetrahedralPoints(int order);
    
    /**
     * @brief Perform numerical integration using Gaussian quadrature
     * @tparam Func Type of function to integrate
     * @param points Integration points
     * @param func Function to integrate (takes xi, eta, zeta as parameters)
     * @return Integral value
     */
    template<typename Func>
    static double integrate(const std::vector<IntegrationPoint>& points, Func func) {
        double result = 0.0;
        for (const auto& point : points) {
            result += point.weight * func(point.xi, point.eta, point.zeta);
        }
        return result;
    }
    
private:
    // 1D Gaussian quadrature data
    static const std::array<std::vector<double>, 6> gauss1DData;
    
    // Triangular integration data
    static const std::array<std::vector<double>, 7> gaussTriangularData;
    
    // Tetrahedral integration data  
    static const std::array<std::vector<double>, 5> gaussTetrahedralData;
};

} // namespace elmer