#pragma once

#include "ElmerCpp.h"
#include <vector>
#include <array>
#include <cmath>
#include <memory>

namespace elmer {

/**
 * @brief Numerical integration utilities for finite element analysis
 * 
 * This module provides Gauss integration points and weights for various
 * element types, as well as shape function evaluation routines.
 * Corresponds to Elmer FEM's Integration.F90 module.
 */
class Integration {
public:
    /**
     * @brief Integration point structure for 1D, 2D, and 3D elements
     */
    struct IntegrationPoint {
        double u;  ///< Natural coordinate in first direction
        double v;  ///< Natural coordinate in second direction
        double w;  ///< Natural coordinate in third direction
        double weight; ///< Integration weight
        
        IntegrationPoint(double u_val = 0.0, double v_val = 0.0, 
                        double w_val = 0.0, double weight_val = 0.0)
            : u(u_val), v(v_val), w(w_val), weight(weight_val) {}
    };

    /**
     * @brief Get Gauss integration points for 1D elements
     * @param npoints Number of integration points
     * @return Vector of integration points with weights
     */
    static std::vector<IntegrationPoint> getGaussPoints1D(int npoints);
    
    /**
     * @brief Get Gauss integration points for 2D quadrilateral elements
     * @param npoints_u Number of points in u direction
     * @param npoints_v Number of points in v direction
     * @return Vector of integration points with weights
     */
    static std::vector<IntegrationPoint> getGaussPoints2D(int npoints_u, int npoints_v);
    
    /**
     * @brief Get Gauss integration points for 3D hexahedral elements
     * @param npoints_u Number of points in u direction
     * @param npoints_v Number of points in v direction
     * @param npoints_w Number of points in w direction
     * @return Vector of integration points with weights
     */
    static std::vector<IntegrationPoint> getGaussPoints3D(int npoints_u, int npoints_v, int npoints_w);
    
    /**
     * @brief Get integration points for triangular elements
     * @param order Integration order (1-7)
     * @return Vector of integration points with weights
     */
    static std::vector<IntegrationPoint> getTrianglePoints(int order);
    
    /**
     * @brief Get integration points for tetrahedral elements
     * @param order Integration order (1-5)
     * @return Vector of integration points with weights
     */
    static std::vector<IntegrationPoint> getTetrahedronPoints(int order);
    
    /**
     * @brief Evaluate linear shape functions for 1D elements
     * @param xi Natural coordinate
     * @param n Node index (0 or 1)
     * @return Shape function value
     */
    static double shapeFunction1DLinear(double xi, int n);
    
    /**
     * @brief Evaluate linear shape functions for 2D quadrilateral elements
     * @param xi Natural coordinate in u direction
     * @param eta Natural coordinate in v direction
     * @param n Node index (0-3)
     * @return Shape function value
     */
    static double shapeFunction2DLinear(double xi, double eta, int n);
    
    /**
     * @brief Evaluate linear shape functions for 3D hexahedral elements
     * @param xi Natural coordinate in u direction
     * @param eta Natural coordinate in v direction
     * @param zeta Natural coordinate in w direction
     * @param n Node index (0-7)
     * @return Shape function value
     */
    static double shapeFunction3DLinear(double xi, double eta, double zeta, int n);
    
    /**
     * @brief Evaluate quadratic shape functions for 1D elements
     * @param xi Natural coordinate
     * @param n Node index (0-2)
     * @return Shape function value
     */
    static double shapeFunction1DQuadratic(double xi, int n);
    
    /**
     * @brief Evaluate shape function derivatives for 1D linear elements
     * @param xi Natural coordinate
     * @param n Node index (0 or 1)
     * @return Shape function derivative
     */
    static double shapeFunctionDerivative1DLinear(double xi, int n);
    
    /**
     * @brief Evaluate shape function derivatives for 2D linear quadrilateral elements
     * @param xi Natural coordinate in u direction
     * @param eta Natural coordinate in v direction
     * @param n Node index (0-3)
     * @param dir Derivative direction (0 for d/dxi, 1 for d/deta)
     * @return Shape function derivative
     */
    static double shapeFunctionDerivative2DLinear(double xi, double eta, int n, int dir);
    
    /**
     * @brief Evaluate shape function derivatives for 3D linear hexahedral elements
     * @param xi Natural coordinate in u direction
     * @param eta Natural coordinate in v direction
     * @param zeta Natural coordinate in w direction
     * @param n Node index (0-7)
     * @param dir Derivative direction (0 for d/dxi, 1 for d/deta, 2 for d/dzeta)
     * @return Shape function derivative
     */
    static double shapeFunctionDerivative3DLinear(double xi, double eta, double zeta, int n, int dir);
    
    /**
     * @brief Compute Jacobian determinant for 1D elements
     * @param nodes Node coordinates
     * @param xi Natural coordinate
     * @return Jacobian determinant
     */
    static double computeJacobian1D(const std::vector<Node>& nodes, double xi);
    
    /**
     * @brief Compute Jacobian determinant for 2D quadrilateral elements
     * @param nodes Node coordinates
     * @param xi Natural coordinate in u direction
     * @param eta Natural coordinate in v direction
     * @return Jacobian determinant
     */
    static double computeJacobian2D(const std::vector<Node>& nodes, double xi, double eta);
    
    /**
     * @brief Compute Jacobian determinant for 3D hexahedral elements
     * @param nodes Node coordinates
     * @param xi Natural coordinate in u direction
     * @param eta Natural coordinate in v direction
     * @param zeta Natural coordinate in w direction
     * @return Jacobian determinant
     */
    static double computeJacobian3D(const std::vector<Node>& nodes, 
                                   double xi, double eta, double zeta);
    
    /**
     * @brief Compute Jacobian matrix for 2D elements
     * @param nodes Node coordinates
     * @param xi Natural coordinate in u direction
     * @param eta Natural coordinate in v direction
     * @return 2x2 Jacobian matrix
     */
    static std::array<std::array<double, 2>, 2> computeJacobianMatrix2D(
        const std::vector<Node>& nodes, double xi, double eta);
    
    /**
     * @brief Compute Jacobian matrix for 3D elements
     * @param nodes Node coordinates
     * @param xi Natural coordinate in u direction
     * @param eta Natural coordinate in v direction
     * @param zeta Natural coordinate in w direction
     * @return 3x3 Jacobian matrix
     */
    static std::array<std::array<double, 3>, 3> computeJacobianMatrix3D(
        const std::vector<Node>& nodes, double xi, double eta, double zeta);
    
    /**
     * @brief Compute inverse Jacobian matrix for 2D elements
     * @param jac Jacobian matrix
     * @return Inverse Jacobian matrix
     */
    static std::array<std::array<double, 2>, 2> computeInverseJacobian2D(
        const std::array<std::array<double, 2>, 2>& jac);
    
    /**
     * @brief Compute inverse Jacobian matrix for 3D elements
     * @param jac Jacobian matrix
     * @return Inverse Jacobian matrix
     */
    static std::array<std::array<double, 3>, 3> computeInverseJacobian3D(
        const std::array<std::array<double, 3>, 3>& jac);
    
    /**
     * @brief Transform derivatives from natural to physical coordinates
     * @param dNdxi Derivatives in natural coordinates
     * @param invJac Inverse Jacobian matrix
     * @param dim Dimension (2 or 3)
     * @return Derivatives in physical coordinates
     */
    static std::vector<double> transformDerivatives(
        const std::vector<double>& dNdxi, 
        const std::array<std::array<double, 3>, 3>& invJac, 
        int dim);

private:
    /**
     * @brief Precompute Gauss points and weights for 1D integration
     */
    static void initializeGaussPoints1D();
    
    /**
     * @brief Get precomputed Gauss points for 1D integration
     * @param npoints Number of points
     * @return Vector of points and weights
     */
    static std::vector<std::pair<double, double>> getPrecomputedGauss1D(int npoints);
    
    // Precomputed Gauss points and weights for 1D integration
    static std::vector<std::vector<std::pair<double, double>>> gaussPoints1D_;
    static bool initialized_;
};

} // namespace elmer