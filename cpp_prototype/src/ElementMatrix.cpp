// ElementMatrix.cpp - Elmer FEM C++ Element Matrix Assembly Module Implementation
// Corresponds to Fortran module: ElementMatrix.F90

#include "ElementMatrix.h"
#include <stdexcept>
#include <cmath>

using namespace elmer;

std::vector<std::vector<double>> ElementMatrix::computeStiffnessMatrix(
    ElementType elementType,
    const std::vector<Node>& nodes,
    const MaterialProperties& material,
    int integrationOrder) {
    
    int dim = 0;
    switch (elementType) {
        case ElementType::LINEAR:
            dim = 1;
            break;
        case ElementType::QUADRATIC:
            dim = 2;
            break;
        case ElementType::HEXAHEDRON:
            dim = 3;
            break;
        default:
            throw std::invalid_argument("Unsupported element type for stiffness matrix");
    }
    
    return integrateElementMatrix(elementType, nodes, integrationOrder,
        [&](double xi, double eta, double zeta, const ShapeFunctions::ShapeResult& shapeResult) {
            
            int nNodes = nodes.size();
            auto B = computeBMatrix(shapeResult.dNdx, shapeResult.dNdy, shapeResult.dNdz, nNodes, dim);
            auto D = computeDMatrix(material.youngsModulus, material.poissonsRatio, dim);
            
            // Compute B^T * D * B * detJ
            int strainSize = (dim == 1) ? 1 : (dim == 2) ? 3 : 6;
            int dofSize = nNodes * dim;
            
            std::vector<std::vector<double>> ke(dofSize, std::vector<double>(dofSize, 0.0));
            
            // B^T * D
            std::vector<std::vector<double>> BT_D(strainSize, std::vector<double>(dofSize, 0.0));
            for (int i = 0; i < strainSize; ++i) {
                for (int j = 0; j < dofSize; ++j) {
                    for (int k = 0; k < strainSize; ++k) {
                        BT_D[i][j] += B[k][j] * D[k][i];
                    }
                }
            }
            
            // (B^T * D) * B
            for (int i = 0; i < dofSize; ++i) {
                for (int j = 0; j < dofSize; ++j) {
                    for (int k = 0; k < strainSize; ++k) {
                        ke[i][j] += BT_D[k][i] * B[k][j];
                    }
                    ke[i][j] *= shapeResult.detJ;
                }
            }
            
            return ke;
        });
}

std::vector<std::vector<double>> ElementMatrix::computeMassMatrix(
    ElementType elementType,
    const std::vector<Node>& nodes,
    const MaterialProperties& material,
    int integrationOrder) {
    
    return integrateElementMatrix(elementType, nodes, integrationOrder,
        [&](double xi, double eta, double zeta, const ShapeFunctions::ShapeResult& shapeResult) {
            
            int nNodes = nodes.size();
            int dim = (elementType == ElementType::LINEAR) ? 1 : 
                     (elementType == ElementType::QUADRATIC) ? 2 : 3;
            int dofSize = nNodes * dim;
            
            std::vector<std::vector<double>> me(dofSize, std::vector<double>(dofSize, 0.0));
            
            // Consistent mass matrix: M = ρ ∫ N^T N dV
            for (int i = 0; i < nNodes; ++i) {
                for (int j = 0; j < nNodes; ++j) {
                    double mij = material.density * shapeResult.values[i] * shapeResult.values[j] * shapeResult.detJ;
                    
                    // For each degree of freedom
                    for (int d = 0; d < dim; ++d) {
                        int row = i * dim + d;
                        int col = j * dim + d;
                        me[row][col] = mij;
                    }
                }
            }
            
            return me;
        });
}

std::vector<double> ElementMatrix::computeLoadVector(
    ElementType elementType,
    const std::vector<Node>& nodes,
    const std::array<double, 3>& bodyForce,
    int integrationOrder) {
    
    return integrateElementVector(elementType, nodes, integrationOrder,
        [&](double xi, double eta, double zeta, const ShapeFunctions::ShapeResult& shapeResult) {
            
            int nNodes = nodes.size();
            int dim = (elementType == ElementType::LINEAR) ? 1 : 
                     (elementType == ElementType::QUADRATIC) ? 2 : 3;
            int dofSize = nNodes * dim;
            
            std::vector<double> fe(dofSize, 0.0);
            
            // Load vector: f = ∫ N^T b dV
            for (int i = 0; i < nNodes; ++i) {
                for (int d = 0; d < dim; ++d) {
                    int index = i * dim + d;
                    fe[index] = shapeResult.values[i] * bodyForce[d] * shapeResult.detJ;
                }
            }
            
            return fe;
        });
}

std::vector<std::vector<double>> ElementMatrix::computeConductivityMatrix(
    ElementType elementType,
    const std::vector<Node>& nodes,
    const MaterialProperties& material,
    int integrationOrder) {
    
    return integrateElementMatrix(elementType, nodes, integrationOrder,
        [&](double xi, double eta, double zeta, const ShapeFunctions::ShapeResult& shapeResult) {
            
            int nNodes = nodes.size();
            int dim = (elementType == ElementType::LINEAR) ? 1 : 
                     (elementType == ElementType::QUADRATIC) ? 2 : 3;
            
            std::vector<std::vector<double>> ke(nNodes, std::vector<double>(nNodes, 0.0));
            
            // Conductivity matrix: K = k ∫ (∇N)^T ∇N dV
            for (int i = 0; i < nNodes; ++i) {
                for (int j = 0; j < nNodes; ++j) {
                    double kij = 0.0;
                    
                    // Dot product of gradients
                    kij += shapeResult.dNdx[i] * shapeResult.dNdx[j];
                    if (dim >= 2) kij += shapeResult.dNdy[i] * shapeResult.dNdy[j];
                    if (dim >= 3) kij += shapeResult.dNdz[i] * shapeResult.dNdz[j];
                    
                    ke[i][j] = material.conductivity * kij * shapeResult.detJ;
                }
            }
            
            return ke;
        });
}

std::vector<std::vector<double>> ElementMatrix::computeCapacityMatrix(
    ElementType elementType,
    const std::vector<Node>& nodes,
    const MaterialProperties& material,
    int integrationOrder) {
    
    return integrateElementMatrix(elementType, nodes, integrationOrder,
        [&](double xi, double eta, double zeta, const ShapeFunctions::ShapeResult& shapeResult) {
            
            int nNodes = nodes.size();
            
            std::vector<std::vector<double>> ce(nNodes, std::vector<double>(nNodes, 0.0));
            
            // Capacity matrix: C = ρc ∫ N^T N dV
            for (int i = 0; i < nNodes; ++i) {
                for (int j = 0; j < nNodes; ++j) {
                    ce[i][j] = material.density * material.specificHeat * 
                              shapeResult.values[i] * shapeResult.values[j] * shapeResult.detJ;
                }
            }
            
            return ce;
        });
}

std::vector<std::vector<double>> ElementMatrix::computeConvectionMatrix(
    ElementType elementType,
    const std::vector<Node>& nodes,
    const std::vector<std::array<double, 3>>& velocity,
    const MaterialProperties& material,
    int integrationOrder) {
    
    return integrateElementMatrix(elementType, nodes, integrationOrder,
        [&](double xi, double eta, double zeta, const ShapeFunctions::ShapeResult& shapeResult) {
            
            int nNodes = nodes.size();
            int dim = (elementType == ElementType::LINEAR) ? 1 : 
                     (elementType == ElementType::QUADRATIC) ? 2 : 3;
            
            // Interpolate velocity at integration point
            std::array<double, 3> v = {0.0, 0.0, 0.0};
            for (int i = 0; i < nNodes; ++i) {
                v[0] += shapeResult.values[i] * velocity[i][0];
                v[1] += shapeResult.values[i] * velocity[i][1];
                v[2] += shapeResult.values[i] * velocity[i][2];
            }
            
            std::vector<std::vector<double>> ce(nNodes, std::vector<double>(nNodes, 0.0));
            
            // Convection matrix: C = ρ ∫ N^T (v·∇N) dV
            for (int i = 0; i < nNodes; ++i) {
                for (int j = 0; j < nNodes; ++j) {
                    double cij = shapeResult.values[i] * 
                                (v[0] * shapeResult.dNdx[j] + 
                                 v[1] * shapeResult.dNdy[j] + 
                                 v[2] * shapeResult.dNdz[j]);
                    
                    ce[i][j] = material.density * cij * shapeResult.detJ;
                }
            }
            
            return ce;
        });
}

// Private helper functions
std::vector<std::vector<double>> ElementMatrix::computeBMatrix(
    const std::vector<double>& dNdx,
    const std::vector<double>& dNdy,
    const std::vector<double>& dNdz,
    int nNodes,
    int dim) {
    
    int strainSize = (dim == 1) ? 1 : (dim == 2) ? 3 : 6;
    int dofSize = nNodes * dim;
    
    std::vector<std::vector<double>> B(strainSize, std::vector<double>(dofSize, 0.0));
    
    if (dim == 1) {
        // 1D: ε = du/dx
        for (int i = 0; i < nNodes; ++i) {
            B[0][i] = dNdx[i];
        }
    } else if (dim == 2) {
        // 2D: ε = [εx, εy, γxy]^T
        for (int i = 0; i < nNodes; ++i) {
            int uIndex = i * 2;
            int vIndex = i * 2 + 1;
            
            B[0][uIndex] = dNdx[i];  // εx
            B[1][vIndex] = dNdy[i];  // εy
            B[2][uIndex] = dNdy[i];  // γxy
            B[2][vIndex] = dNdx[i];  // γxy
        }
    } else if (dim == 3) {
        // 3D: ε = [εx, εy, εz, γxy, γyz, γzx]^T
        for (int i = 0; i < nNodes; ++i) {
            int uIndex = i * 3;
            int vIndex = i * 3 + 1;
            int wIndex = i * 3 + 2;
            
            B[0][uIndex] = dNdx[i];  // εx
            B[1][vIndex] = dNdy[i];  // εy
            B[2][wIndex] = dNdz[i];  // εz
            B[3][uIndex] = dNdy[i];  // γxy
            B[3][vIndex] = dNdx[i];  // γxy
            B[4][vIndex] = dNdz[i];  // γyz
            B[4][wIndex] = dNdy[i];  // γyz
            B[5][uIndex] = dNdz[i];  // γzx
            B[5][wIndex] = dNdx[i];  // γzx
        }
    }
    
    return B;
}

std::vector<std::vector<double>> ElementMatrix::computeDMatrix(
    double youngsModulus,
    double poissonsRatio,
    int dim) {
    
    int strainSize = (dim == 1) ? 1 : (dim == 2) ? 3 : 6;
    std::vector<std::vector<double>> D(strainSize, std::vector<double>(strainSize, 0.0));
    
    if (dim == 1) {
        // 1D: D = E
        D[0][0] = youngsModulus;
    } else if (dim == 2) {
        // 2D plane stress
        double factor = youngsModulus / (1.0 - poissonsRatio * poissonsRatio);
        D[0][0] = factor;
        D[0][1] = factor * poissonsRatio;
        D[1][0] = factor * poissonsRatio;
        D[1][1] = factor;
        D[2][2] = factor * (1.0 - poissonsRatio) / 2.0;
    } else if (dim == 3) {
        // 3D isotropic
        double lambda = youngsModulus * poissonsRatio / ((1.0 + poissonsRatio) * (1.0 - 2.0 * poissonsRatio));
        double mu = youngsModulus / (2.0 * (1.0 + poissonsRatio));
        
        D[0][0] = lambda + 2.0 * mu;
        D[0][1] = lambda;
        D[0][2] = lambda;
        
        D[1][0] = lambda;
        D[1][1] = lambda + 2.0 * mu;
        D[1][2] = lambda;
        
        D[2][0] = lambda;
        D[2][1] = lambda;
        D[2][2] = lambda + 2.0 * mu;
        
        D[3][3] = mu;
        D[4][4] = mu;
        D[5][5] = mu;
    }
    
    return D;
}

std::vector<std::vector<double>> ElementMatrix::integrateElementMatrix(
    ElementType elementType,
    const std::vector<Node>& nodes,
    int integrationOrder,
    std::function<std::vector<std::vector<double>>(double, double, double, const ShapeFunctions::ShapeResult&)> matrixFunction) {
    
    std::vector<GaussIntegration::IntegrationPoint> points;
    
    // Get appropriate integration points based on element type
    switch (elementType) {
        case ElementType::LINEAR:
            points = GaussIntegration::get1DPoints(integrationOrder);
            break;
        case ElementType::QUADRATIC:
            points = GaussIntegration::getQuadrilateralPoints(integrationOrder);
            break;
        case ElementType::HEXAHEDRON:
            points = GaussIntegration::getHexahedralPoints(integrationOrder);
            break;
        default:
            throw std::invalid_argument("Unsupported element type for integration");
    }
    
    // Initialize result matrix
    auto sampleResult = matrixFunction(0.0, 0.0, 0.0, ShapeFunctions::ShapeResult(nodes.size()));
    std::vector<std::vector<double>> result(sampleResult.size(), 
                                           std::vector<double>(sampleResult[0].size(), 0.0));
    
    // Perform numerical integration
    for (const auto& point : points) {
        auto shapeResult = ShapeFunctions::computeShapeFunctions(elementType, nodes, 
                                                               point.xi, point.eta, point.zeta);
        auto matrix = matrixFunction(point.xi, point.eta, point.zeta, shapeResult);
        
        // Accumulate weighted matrix
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[i].size(); ++j) {
                result[i][j] += point.weight * matrix[i][j];
            }
        }
    }
    
    return result;
}

std::vector<double> ElementMatrix::integrateElementVector(
    ElementType elementType,
    const std::vector<Node>& nodes,
    int integrationOrder,
    std::function<std::vector<double>(double, double, double, const ShapeFunctions::ShapeResult&)> vectorFunction) {
    
    std::vector<GaussIntegration::IntegrationPoint> points;
    
    // Get appropriate integration points based on element type
    switch (elementType) {
        case ElementType::LINEAR:
            points = GaussIntegration::get1DPoints(integrationOrder);
            break;
        case ElementType::QUADRATIC:
            points = GaussIntegration::getQuadrilateralPoints(integrationOrder);
            break;
        case ElementType::HEXAHEDRON:
            points = GaussIntegration::getHexahedralPoints(integrationOrder);
            break;
        default:
            throw std::invalid_argument("Unsupported element type for integration");
    }
    
    // Initialize result vector
    auto sampleResult = vectorFunction(0.0, 0.0, 0.0, ShapeFunctions::ShapeResult(nodes.size()));
    std::vector<double> result(sampleResult.size(), 0.0);
    
    // Perform numerical integration
    for (const auto& point : points) {
        auto shapeResult = ShapeFunctions::computeShapeFunctions(elementType, nodes, 
                                                               point.xi, point.eta, point.zeta);
        auto vector = vectorFunction(point.xi, point.eta, point.zeta, shapeResult);
        
        // Accumulate weighted vector
        for (size_t i = 0; i < vector.size(); ++i) {
            result[i] += point.weight * vector[i];
        }
    }
    
    return result;
}