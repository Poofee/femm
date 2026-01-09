// MatrixAssembly.h - Elmer FEM C++ Matrix Assembly Module
// Corresponds to Fortran module: MatrixAssembly.F90

#pragma once

#include "ElmerCpp.h"
#include <memory>
#include <vector>
#include <stdexcept>

namespace elmer {

// Matrix format enumeration
enum class MatrixFormat {
    CRS,        // Compressed Row Storage format
    LIST,       // List storage format
    BAND,       // Banded matrix format
    SBAND       // Symmetric banded matrix format
};

// Matrix assembly class
class MatrixAssembly {
public:
    // Constructor
    MatrixAssembly() = default;
    
    // Set matrix element
    static void SetMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                elmer::Integer i, elmer::Integer j, elmer::Real value);
    
    // Get matrix element
    static elmer::Real GetMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                elmer::Integer i, elmer::Integer j);
    
    // Modify matrix element and return old value
    static elmer::Real ChangeMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                   elmer::Integer i, elmer::Integer j, elmer::Real newValue);
    
    // Add to matrix element
    static void AddToMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                  elmer::Integer i, elmer::Integer j, elmer::Real value);
    
    // Complex matrix element operation
    static void AddToComplexMatrixElement(std::shared_ptr<Matrix> matrix,
                                         elmer::Integer rowId, elmer::Integer colId,
                                         elmer::Real re, elmer::Real im);
    
    // Move matrix element
    static void MoveMatrixElement(std::shared_ptr<Matrix> matrix,
                                 elmer::Integer i1, elmer::Integer j1, elmer::Integer i2, elmer::Integer j2);
    
    // Zero matrix
    static void ZeroMatrix(std::shared_ptr<Matrix> matrix);
    
    // Scale matrix
    static void ScaleMatrix(std::shared_ptr<Matrix> matrix, elmer::Real factor);
    
    // Matrix-vector multiplication
    static void MatrixVectorMultiply(std::shared_ptr<elmer::Matrix> matrix,
                                    std::shared_ptr<elmer::Vector> vector,
                                    std::shared_ptr<elmer::Vector> result);
    
    // Finite element matrix assembly related functions
    
    // Assemble element stiffness matrix to global matrix
    static void AssembleElementMatrix(std::shared_ptr<Matrix> globalMatrix,
                                     const std::vector<elmer::Integer>& dofIndices,
                                     const std::vector<std::vector<elmer::Real>>& elementMatrix);
    
    // Assemble element load vector to global vector
    static void AssembleElementVector(std::shared_ptr<elmer::Vector> globalVector,
                                     const std::vector<elmer::Integer>& dofIndices,
                                     const std::vector<elmer::Real>& elementVector);
    
    // Apply Dirichlet boundary conditions
    static void ApplyDirichletBC(std::shared_ptr<elmer::Matrix> matrix,
                                std::vector<elmer::Real>& rhs,
                                elmer::Integer dofIndex,
                                elmer::Real prescribedValue);
    
    // Apply Dirichlet boundary conditions (multiple DOFs)
    static void ApplyDirichletBC(std::shared_ptr<elmer::Matrix> matrix,
                                std::vector<elmer::Real>& rhs,
                                const std::vector<elmer::Integer>& dofIndices,
                                const std::vector<elmer::Real>& prescribedValues);
    
    // Matrix format conversion
    static std::shared_ptr<elmer::Matrix> ConvertMatrixFormat(std::shared_ptr<elmer::Matrix> source,
                                                      elmer::MatrixFormat targetFormat);
    
    // Matrix information query
    static elmer::Integer GetMatrixSize(std::shared_ptr<elmer::Matrix> matrix);
    static elmer::Integer GetMatrixNonzeros(std::shared_ptr<elmer::Matrix> matrix);
    static elmer::MatrixFormat GetMatrixFormat(std::shared_ptr<elmer::Matrix> matrix);
    
    // Get matrix row
    static std::vector<elmer::Real> GetMatrixRow(std::shared_ptr<elmer::Matrix> matrix, elmer::Integer row);
    
    // Set matrix row
    static void SetMatrixRow(std::shared_ptr<elmer::Matrix> matrix, elmer::Integer row,
                            const std::vector<elmer::Real>& values);
    
    // Matrix condition number estimation
    static elmer::Real EstimateConditionNumber(std::shared_ptr<elmer::Matrix> matrix);
    
    // Matrix symmetry check
    static bool IsSymmetric(std::shared_ptr<elmer::Matrix> matrix, elmer::Real tolerance = 1e-12);
    
    // Matrix positive definiteness check
    static bool IsPositiveDefinite(std::shared_ptr<elmer::Matrix> matrix);
    
private:
    // Internal helper functions
    
    // Check index validity
    static void CheckIndices(std::shared_ptr<elmer::Matrix> matrix, elmer::Integer i, elmer::Integer j);
    
    // CRS matrix specific operations
    static void CRS_SetMatrixElement(std::shared_ptr<Matrix> matrix,
                                    elmer::Integer i, elmer::Integer j, elmer::Real value);
    static elmer::Real CRS_GetMatrixElement(std::shared_ptr<Matrix> matrix,
                                    elmer::Integer i, elmer::Integer j);
    static void CRS_AddToMatrixElement(std::shared_ptr<Matrix> matrix,
                                      elmer::Integer i, elmer::Integer j, elmer::Real value);
    
    // Band matrix specific operations
    static void Band_SetMatrixElement(std::shared_ptr<Matrix> matrix,
                                     elmer::Integer i, elmer::Integer j, elmer::Real value);
    static elmer::Real Band_GetMatrixElement(std::shared_ptr<Matrix> matrix,
                                     elmer::Integer i, elmer::Integer j);
    static void Band_AddToMatrixElement(std::shared_ptr<Matrix> matrix,
                                       elmer::Integer i, elmer::Integer j, elmer::Real value);
    
    // List matrix specific operations
    static void List_SetMatrixElement(std::shared_ptr<Matrix> matrix,
                                     elmer::Integer i, elmer::Integer j, elmer::Real value);
    static elmer::Real List_GetMatrixElement(std::shared_ptr<Matrix> matrix,
                                     elmer::Integer i, elmer::Integer j);
    static void List_AddToMatrixElement(std::shared_ptr<Matrix> matrix,
                                       elmer::Integer i, elmer::Integer j, elmer::Real value);
};

} // namespace elmer