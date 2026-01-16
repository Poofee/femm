// CRSMatrix.h - Compressed Row Storage Matrix Implementation
// Corresponds to Elmer FEM's CRSMatrix.F90

#pragma once

#include "Types.h"
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace elmer {

/**
 * @brief High-performance CRS (Compressed Row Storage) matrix
 * 
 * This class implements the core sparse matrix functionality used in Elmer FEM.
 * It corresponds to the functionality in CRSMatrix.F90.
 */
class CRSMatrix : public Matrix {
private:
    Integer nrows_;                    ///< Number of rows
    Integer ncols_;                    ///< Number of columns
    std::vector<Real> values_;         ///< Non-zero values
    std::vector<Integer> col_indices_; ///< Column indices for non-zero values
    std::vector<Integer> row_pointers_; ///< Row pointers (start index for each row)
    
    // Internal storage for diagonal entries (optimization)
    std::vector<Real> diagonal_;
    bool diagonal_computed_;
    
public:
    /**
     * @brief Constructor
     * @param nrows Number of rows
     * @param ncols Number of columns
     */
    CRSMatrix(Integer nrows, Integer ncols);
    
    ~CRSMatrix() override = default;
    
    // Matrix interface implementation
    Integer GetNumRows() const override { return nrows_; }
    Integer GetNumCols() const override { return ncols_; }
    
    void Zero() override;
    
    Real GetElement(Integer i, Integer j) const override;
    void SetElement(Integer i, Integer j, Real value) override;
    void AddToElement(Integer i, Integer j, Real value) override;
    
    // Clone method for creating a copy of the matrix
    std::unique_ptr<Matrix> Clone() const;
    
    // CRS-specific methods
    /**
     * @brief Get the number of non-zero elements
     */
    Integer GetNonZeroCount() const { return static_cast<Integer>(values_.size()); }
    
    /**
     * @brief Get row pointer for row i
     */
    Integer GetRowPointer(Integer i) const;
    
    /**
     * @brief Get column indices for row i
     */
    std::vector<Integer> GetColumnIndices(Integer i) const;
    
    /**
     * @brief Get non-zero values for row i
     */
    std::vector<Real> GetRowValues(Integer i) const;
    
    /**
     * @brief Matrix-vector multiplication: y = A * x
     */
    void Multiply(const Vector& x, Vector& y) const;
    
    /**
     * @brief Matrix-vector multiplication with scaling: y = alpha * A * x + beta * y
     */
    void Multiply(Real alpha, const Vector& x, Real beta, Vector& y) const;
    
    /**
     * @brief Transpose matrix-vector multiplication: y = A^T * x
     */
    void MultiplyTranspose(const Vector& x, Vector& y) const;
    
    /**
     * @brief Get diagonal element (optimized access)
     */
    Real GetDiagonal(Integer i) const;
    
    /**
     * @brief Compute diagonal entries (for optimization)
     */
    void ComputeDiagonal();
    
    /**
     * @brief Check if matrix is symmetric (within tolerance)
     */
    bool IsSymmetric(Real tolerance = 1e-12) const;
    
    /**
     * @brief Get memory usage in bytes
     */
    std::size_t GetMemoryUsage() const;
    
    /**
     * @brief Print matrix structure for debugging
     */
    void PrintStructure() const;
    
    /**
     * @brief Get row pointers (for internal use and conversion)
     */
    const std::vector<Integer>& GetRowPointers() const { return row_pointers_; }
    
    /**
     * @brief Get column indices (for internal use and conversion)
     */
    const std::vector<Integer>& GetColumnIndices() const { return col_indices_; }
    
    /**
     * @brief Get non-zero values (for internal use and conversion)
     */
    const std::vector<Real>& GetValues() const { return values_; }
    
private:
    /**
     * @brief Find position of element (i,j) in CRS storage
     * @return Index in values_ array, or -1 if not found
     */
    Integer FindElementPosition(Integer i, Integer j) const;
    
    /**
     * @brief Insert element at specific position
     */
    void InsertElement(Integer i, Integer j, Real value, Integer pos);
    
    /**
     * @brief Validate indices
     */
    void CheckIndices(Integer i, Integer j) const;
};

/**
 * @brief Factory function for creating CRS matrices
 */
std::unique_ptr<Matrix> CreateCRSMatrix(Integer nrows, Integer ncols);

} // namespace elmer