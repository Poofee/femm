// ElementMatrix.h - Elmer FEM C++ Element Matrix Assembly Module
// Corresponds to Fortran module: ElementMatrix.F90

#pragma once

#include "core/base/Types.h"
#include "core/utils/ShapeFunctions.h"
#include "core/utils/GaussIntegration.h"
#include <vector>
#include <functional>

namespace elmer {

/**
 * @brief Element matrix assembly for finite element analysis
 * 
 * This module provides element matrix assembly for various physical fields
 * including stiffness matrices, mass matrices, and load vectors.
 */
class ElementMatrix {
public:
    /**
     * @brief Material properties for element matrix computation
     */
    struct MaterialProperties {
        double youngsModulus;      ///< Young's modulus (for elasticity)
        double poissonsRatio;      ///< Poisson's ratio (for elasticity)
        double density;            ///< Density (for mass matrix)
        double conductivity;       ///< Thermal conductivity (for heat transfer)
        double specificHeat;       ///< Specific heat capacity
        double viscosity;          ///< Viscosity (for fluid dynamics)
        
        MaterialProperties() 
            : youngsModulus(1.0), poissonsRatio(0.3), density(1.0), 
              conductivity(1.0), specificHeat(1.0), viscosity(1.0) {}
    };

    /**
     * @brief Compute element stiffness matrix for linear elasticity
     * @param elementType Type of element
     * @param nodes Element nodes
     * @param material Material properties
     * @param integrationOrder Integration order
     * @return Element stiffness matrix
     */
    static std::vector<std::vector<double>> computeStiffnessMatrix(
        ElementType elementType,
        const std::vector<Node>& nodes,
        const MaterialProperties& material,
        int integrationOrder = 2);

    /**
     * @brief Compute element mass matrix
     * @param elementType Type of element
     * @param nodes Element nodes
     * @param material Material properties
     * @param integrationOrder Integration order
     * @return Element mass matrix
     */
    static std::vector<std::vector<double>> computeMassMatrix(
        ElementType elementType,
        const std::vector<Node>& nodes,
        const MaterialProperties& material,
        int integrationOrder = 2);

    /**
     * @brief Compute element load vector for body forces
     * @param elementType Type of element
     * @param nodes Element nodes
     * @param bodyForce Body force vector (fx, fy, fz)
     * @param integrationOrder Integration order
     * @return Element load vector
     */
    static std::vector<double> computeLoadVector(
        ElementType elementType,
        const std::vector<Node>& nodes,
        const std::array<double, 3>& bodyForce,
        int integrationOrder = 2);

    /**
     * @brief Compute element conductivity matrix for heat transfer
     * @param elementType Type of element
     * @param nodes Element nodes
     * @param material Material properties
     * @param integrationOrder Integration order
     * @return Element conductivity matrix
     */
    static std::vector<std::vector<double>> computeConductivityMatrix(
        ElementType elementType,
        const std::vector<Node>& nodes,
        const MaterialProperties& material,
        int integrationOrder = 2);

    /**
     * @brief Compute element capacity matrix for transient problems
     * @param elementType Type of element
     * @param nodes Element nodes
     * @param material Material properties
     * @param integrationOrder Integration order
     * @return Element capacity matrix
     */
    static std::vector<std::vector<double>> computeCapacityMatrix(
        ElementType elementType,
        const std::vector<Node>& nodes,
        const MaterialProperties& material,
        int integrationOrder = 2);

    /**
     * @brief Compute element convection matrix for fluid dynamics
     * @param elementType Type of element
     * @param nodes Element nodes
     * @param velocity Velocity field at nodes
     * @param material Material properties
     * @param integrationOrder Integration order
     * @return Element convection matrix
     */
    static std::vector<std::vector<double>> computeConvectionMatrix(
        ElementType elementType,
        const std::vector<Node>& nodes,
        const std::vector<std::array<double, 3>>& velocity,
        const MaterialProperties& material,
        int integrationOrder = 2);

private:
    /**
     * @brief Compute B matrix for elasticity (strain-displacement matrix)
     * @param dNdx Shape function derivatives in global coordinates
     * @param nNodes Number of nodes
     * @param dim Spatial dimension
     * @return B matrix
     */
    static std::vector<std::vector<double>> computeBMatrix(
        const std::vector<double>& dNdx,
        const std::vector<double>& dNdy,
        const std::vector<double>& dNdz,
        int nNodes,
        int dim);

    /**
     * @brief Compute D matrix for elasticity (constitutive matrix)
     * @param youngsModulus Young's modulus
     * @param poissonsRatio Poisson's ratio
     * @param dim Spatial dimension
     * @return D matrix
     */
    static std::vector<std::vector<double>> computeDMatrix(
        double youngsModulus,
        double poissonsRatio,
        int dim);

    /**
     * @brief Compute element matrix using numerical integration
     * @param elementType Type of element
     * @param nodes Element nodes
     * @param integrationOrder Integration order
     * @param matrixFunction Function to compute matrix at integration point
     * @return Element matrix
     */
    static std::vector<std::vector<double>> integrateElementMatrix(
        ElementType elementType,
        const std::vector<Node>& nodes,
        int integrationOrder,
        std::function<std::vector<std::vector<double>>(double, double, double, const ShapeFunctions::ShapeResult&)> matrixFunction);

    /**
     * @brief Compute element vector using numerical integration
     * @param elementType Type of element
     * @param nodes Element nodes
     * @param integrationOrder Integration order
     * @param vectorFunction Function to compute vector at integration point
     * @return Element vector
     */
    static std::vector<double> integrateElementVector(
        ElementType elementType,
        const std::vector<Node>& nodes,
        int integrationOrder,
        std::function<std::vector<double>(double, double, double, const ShapeFunctions::ShapeResult&)> vectorFunction);
};

} // namespace elmer

