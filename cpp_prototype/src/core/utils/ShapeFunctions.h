#pragma once

#include "Types.h"
#include "Mesh.h"
#include <vector>
#include <array>
#include <functional>

namespace elmer {

/**
 * @brief Shape function system for finite element analysis
 * 
 * This module provides shape functions and their derivatives for various
 * element types, corresponding to Elmer FEM's ElementDescription module.
 */
class ShapeFunctions {
public:
    /**
     * @brief Shape function evaluation result
     */
    struct ShapeResult {
        std::vector<double> values;        ///< Shape function values
        std::vector<double> dNdxi;        ///< Derivatives in xi direction
        std::vector<double> dNdeta;       ///< Derivatives in eta direction
        std::vector<double> dNdzeta;      ///< Derivatives in zeta direction
        std::vector<double> dNdx;         ///< Derivatives in x direction
        std::vector<double> dNdy;         ///< Derivatives in y direction
        std::vector<double> dNdz;         ///< Derivatives in z direction
        double detJ;                      ///< Jacobian determinant
        
        ShapeResult(size_t nNodes) 
            : values(nNodes, 0.0), dNdxi(nNodes, 0.0), dNdeta(nNodes, 0.0),
              dNdzeta(nNodes, 0.0), dNdx(nNodes, 0.0), dNdy(nNodes, 0.0),
              dNdz(nNodes, 0.0), detJ(0.0) {}
    };

    /**
     * @brief Evaluate linear shape functions for 1D elements
     * @param xi Natural coordinate (-1 to 1)
     * @param nNodes Number of nodes (2 for linear, 3 for quadratic)
     * @return Shape function values
     */
    static std::vector<double> linear1D(double xi, int nNodes = 2);
    
    /**
     * @brief Evaluate derivatives of linear shape functions for 1D elements
     * @param xi Natural coordinate (-1 to 1)
     * @param nNodes Number of nodes (2 for linear, 3 for quadratic)
     * @return Shape function derivatives
     */
    static std::vector<double> linear1DDerivatives(double xi, int nNodes = 2);
    
    /**
     * @brief Evaluate linear shape functions for 2D quadrilateral elements
     * @param xi Natural coordinate in u direction (-1 to 1)
     * @param eta Natural coordinate in v direction (-1 to 1)
     * @param nNodes Number of nodes (4 for linear, 8 for quadratic)
     * @return Shape function values
     */
    static std::vector<double> linearQuadrilateral(double xi, double eta, int nNodes = 4);
    
    /**
     * @brief Evaluate derivatives of linear shape functions for 2D quadrilateral elements
     * @param xi Natural coordinate in u direction (-1 to 1)
     * @param eta Natural coordinate in v direction (-1 to 1)
     * @param nNodes Number of nodes (4 for linear, 8 for quadratic)
     * @return Shape function derivatives (dNdxi, dNdeta)
     */
    static std::pair<std::vector<double>, std::vector<double>> 
    linearQuadrilateralDerivatives(double xi, double eta, int nNodes = 4);
    
    /**
     * @brief Evaluate linear shape functions for 3D hexahedral elements
     * @param xi Natural coordinate in u direction (-1 to 1)
     * @param eta Natural coordinate in v direction (-1 to 1)
     * @param zeta Natural coordinate in w direction (-1 to 1)
     * @param nNodes Number of nodes (8 for linear, 20 for quadratic)
     * @return Shape function values
     */
    static std::vector<double> linearHexahedron(double xi, double eta, double zeta, int nNodes = 8);
    
    /**
     * @brief Evaluate derivatives of linear shape functions for 3D hexahedral elements
     * @param xi Natural coordinate in u direction (-1 to 1)
     * @param eta Natural coordinate in v direction (-1 to 1)
     * @param zeta Natural coordinate in w direction (-1 to 1)
     * @param nNodes Number of nodes (8 for linear, 20 for quadratic)
     * @return Shape function derivatives (dNdxi, dNdeta, dNdzeta)
     */
    static std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    linearHexahedronDerivatives(double xi, double eta, double zeta, int nNodes = 8);
    
    /**
     * @brief Evaluate linear shape functions for triangular elements
     * @param xi Natural coordinate in area coordinates
     * @param eta Natural coordinate in area coordinates
     * @param nNodes Number of nodes (3 for linear, 6 for quadratic)
     * @return Shape function values
     */
    static std::vector<double> linearTriangle(double xi, double eta, int nNodes = 3);
    
    /**
     * @brief Evaluate derivatives of linear shape functions for triangular elements
     * @param xi Natural coordinate in area coordinates
     * @param eta Natural coordinate in area coordinates
     * @param nNodes Number of nodes (3 for linear, 6 for quadratic)
     * @return Shape function derivatives (dNdxi, dNdeta)
     */
    static std::pair<std::vector<double>, std::vector<double>> 
    linearTriangleDerivatives(double xi, double eta, int nNodes = 3);
    
    /**
     * @brief Evaluate linear shape functions for tetrahedral elements
     * @param xi Natural coordinate in volume coordinates
     * @param eta Natural coordinate in volume coordinates
     * @param zeta Natural coordinate in volume coordinates
     * @param nNodes Number of nodes (4 for linear, 10 for quadratic)
     * @return Shape function values
     */
    static std::vector<double> linearTetrahedron(double xi, double eta, double zeta, int nNodes = 4);
    
    /**
     * @brief Evaluate derivatives of linear shape functions for tetrahedral elements
     * @param xi Natural coordinate in volume coordinates
     * @param eta Natural coordinate in volume coordinates
     * @param zeta Natural coordinate in volume coordinates
     * @param nNodes Number of nodes (4 for linear, 10 for quadratic)
     * @return Shape function derivatives (dNdxi, dNdeta, dNdzeta)
     */
    static std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    linearTetrahedronDerivatives(double xi, double eta, double zeta, int nNodes = 4);
    
    /**
     * @brief Compute shape functions and their global derivatives
     * @param elementType Type of element
     * @param nodes Node coordinates
     * @param xi Natural coordinate in u direction
     * @param eta Natural coordinate in v direction
     * @param zeta Natural coordinate in w direction
     * @return Complete shape function result with global derivatives
     */
    static ShapeResult computeShapeFunctions(ElementType elementType, 
                                            const std::vector<Node>& nodes,
                                            double xi, double eta, double zeta);
    
    /**
     * @brief Compute Jacobian matrix for coordinate transformation
     * @param nodes Node coordinates
     * @param dNdxi Derivatives in xi direction
     * @param dNdeta Derivatives in eta direction
     * @param dNdzeta Derivatives in zeta direction
     * @return 3x3 Jacobian matrix
     */
    static std::array<std::array<double, 3>, 3> 
    computeJacobianMatrix(const std::vector<Node>& nodes,
                         const std::vector<double>& dNdxi,
                         const std::vector<double>& dNdeta,
                         const std::vector<double>& dNdzeta);
    
    /**
     * @brief Compute inverse of Jacobian matrix
     * @param jac Jacobian matrix
     * @return Inverse Jacobian matrix
     */
    static std::array<std::array<double, 3>, 3> 
    computeInverseJacobian(const std::array<std::array<double, 3>, 3>& jac);
    
    /**
     * @brief Transform derivatives from natural to physical coordinates
     * @param invJac Inverse Jacobian matrix
     * @param dNdxi Derivatives in xi direction
     * @param dNdeta Derivatives in eta direction
     * @param dNdzeta Derivatives in zeta direction
     * @return Derivatives in physical coordinates (dNdx, dNdy, dNdz)
     */
    static std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    transformDerivatives(const std::array<std::array<double, 3>, 3>& invJac,
                        const std::vector<double>& dNdxi,
                        const std::vector<double>& dNdeta,
                        const std::vector<double>& dNdzeta);
    
    /**
     * @brief Compute determinant of Jacobian matrix
     * @param jac Jacobian matrix
     * @return Jacobian determinant
     */
    static double computeJacobianDeterminant(const std::array<std::array<double, 3>, 3>& jac);
    
    /**
     * @brief Get number of nodes for given element type and order
     * @param elementType Type of element
     * @param order Polynomial order (1 for linear, 2 for quadratic)
     * @return Number of nodes
     */
    static int getNumberOfNodes(ElementType elementType, int order = 1);
    
    /**
     * @brief Check if element type is supported
     * @param elementType Type of element
     * @return True if supported, false otherwise
     */
    static bool isSupported(ElementType elementType);
    
    /**
     * @brief Get dimension of element type
     * @param elementType Type of element
     * @return Dimension (1, 2, or 3)
     */
    static int getDimension(ElementType elementType);
    
private:
    /**
     * @brief Compute shape functions for 1D elements
     */
    static ShapeResult compute1DShapeFunctions(const std::vector<Node>& nodes,
                                              double xi);
    
    /**
     * @brief Compute shape functions for 2D quadrilateral elements
     */
    static ShapeResult compute2DQuadrilateralShapeFunctions(const std::vector<Node>& nodes,
                                                           double xi, double eta);
    
    /**
     * @brief Compute shape functions for 3D hexahedral elements
     */
    static ShapeResult compute3DHexahedronShapeFunctions(const std::vector<Node>& nodes,
                                                        double xi, double eta, double zeta);
    
    /**
     * @brief Compute shape functions for triangular elements
     */
    static ShapeResult computeTriangleShapeFunctions(const std::vector<Node>& nodes,
                                                    double xi, double eta);
    
    /**
     * @brief Compute shape functions for tetrahedral elements
     */
    static ShapeResult computeTetrahedronShapeFunctions(const std::vector<Node>& nodes,
                                                       double xi, double eta, double zeta);
    
    /**
     * @brief Compute shape functions for prism elements
     */
    static ShapeResult computePrismShapeFunctions(const std::vector<Node>& nodes,
                                                 double xi, double eta, double zeta);
    
    // Wedge/Prism shape functions (corresponding to Fortran WedgeNodalPBasis functions)
    
    /**
     * @brief Evaluate wedge nodal basis function at point (u,v,w)
     * @param node Node number (1-6)
     * @param u Natural coordinate u
     * @param v Natural coordinate v
     * @param w Natural coordinate w
     * @return Value of wedge nodal function
     */
    static double wedgeNodalPBasis(int node, double u, double v, double w);
    
    /**
     * @brief Evaluate all wedge nodal basis functions at point (u,v,w)
     * @param u Natural coordinate u
     * @param v Natural coordinate v
     * @param w Natural coordinate w
     * @return Vector of 6 shape function values
     */
    static std::vector<double> wedgeNodalPBasisAll(double u, double v, double w);
    
    /**
     * @brief Evaluate all wedge linear basis functions at point (u,v,w)
     * @param u Natural coordinate u
     * @param v Natural coordinate v
     * @param w Natural coordinate w
     * @return Vector of 6 linear shape function values
     */
    static std::vector<double> wedgeNodalLBasisAll(double u, double v, double w);
    
    /**
     * @brief Evaluate gradient of wedge nodal basis function
     * @param node Node number (1-6)
     * @param u Natural coordinate u
     * @param v Natural coordinate v
     * @param w Natural coordinate w
     * @return Gradient vector (3 components)
     */
    static std::vector<double> dWedgeNodalPBasis(int node, double u, double v, double w);
    
    /**
     * @brief Evaluate gradients of all wedge nodal basis functions
     * @param u Natural coordinate u
     * @param v Natural coordinate v
     * @param w Natural coordinate w
     * @return Matrix of gradients (6x3)
     */
    static std::vector<std::vector<double>> dWedgeNodalPBasisAll(double u, double v, double w);
    
    /**
     * @brief Evaluate gradients of all wedge linear basis functions
     * @param u Natural coordinate u
     * @param v Natural coordinate v
     * @param w Natural coordinate w
     * @return Matrix of gradients (6x3)
     */
    static std::vector<std::vector<double>> dWedgeNodalLBasisAll(double u, double v, double w);
    
    /**
     * @brief Evaluate second derivatives of wedge nodal basis function
     * @param node Node number (1-6)
     * @param u Natural coordinate u
     * @param v Natural coordinate v
     * @param w Natural coordinate w
     * @return Hessian matrix (3x3)
     */
    static std::vector<std::vector<double>> ddWedgeNodalPBasis(int node, double u, double v, double w);
    
    // Helper functions for wedge shape functions
    
    /**
     * @brief Evaluate wedge linear basis function
     * @param node Node number (1-3)
     * @param u Natural coordinate u
     * @param v Natural coordinate v
     * @return Value of wedge linear function
     */
    static double wedgeL(int node, double u, double v);
    
    /**
     * @brief Evaluate gradient of wedge linear basis function
     * @param node Node number (1-3)
     * @param u Natural coordinate u
     * @param v Natural coordinate v
     * @return Gradient vector (2 components)
     */
    static std::vector<double> dWedgeL(int node, double u, double v);
};

} // namespace elmer