// Matrix.cpp - Linear algebra matrix implementation
// Corresponds to Elmer FEM's CRSMatrix.F90 functionality

#include "ElmerCpp.h"
#include "../eigen-5.0.1/Eigen/Sparse"
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>

namespace elmer {

// Simple Vector implementation
class SimpleVector : public Vector {
public:
    SimpleVector(Integer size) : data_(size, 0.0) {}
    
    Integer Size() const override { return data_.size(); }
    
    Real& operator[](Integer i) override { return data_[i]; }
    
    const Real& operator[](Integer i) const override { return data_[i]; }
    
    void Zero() override { std::fill(data_.begin(), data_.end(), 0.0); }
    
private:
    std::vector<Real> data_;
};

// Vector factory method implementation
std::unique_ptr<Vector> Vector::Create(Integer size) {
    return std::make_unique<SimpleVector>(size);
}

// CRS matrix implementation (corresponds to CRSMatrix.F90)
class CRSMatrix : public Matrix {
public:
    CRSMatrix(Integer nrows, Integer ncols) 
        : nrows_(nrows), ncols_(ncols) {
        rows_.resize(nrows + 1, 0);
    }
    
    ~CRSMatrix() override = default;
    
    Integer GetNumRows() const override { return nrows_; }
    Integer GetNumCols() const override { return ncols_; }
    
    void Zero() override {
        std::fill(values_.begin(), values_.end(), 0.0);
    }
    
    void SetElement(Integer i, Integer j, Real value) override {
        // Find or insert element
        auto pos = FindElementPosition(i, j);
        if (pos != values_.end()) {
            *pos = value;
        } else {
            InsertElement(i, j, value);
        }
    }
    
    Real GetElement(Integer i, Integer j) const override {
        auto pos = FindElementPosition(i, j);
        return (pos != values_.end()) ? *pos : 0.0;
    }
    
    void AddToElement(Integer i, Integer j, Real value) override {
        auto pos = FindElementPosition(i, j);
        if (pos != values_.end()) {
            *pos += value;
        } else {
            InsertElement(i, j, value);
        }
    }
    
    // Get matrix-vector product (corresponds to Fortran matrix-vector operations)
    void Multiply(const Vector& x, Vector& y) const override {
        for (Integer i = 0; i < nrows_; ++i) {
            Real sum = 0.0;
            for (Integer k = rows_[i]; k < rows_[i + 1]; ++k) {
                sum += values_[k] * x[cols_[k]];
            }
            y[i] = sum;
        }
    }
    
    // Matrix-vector multiplication with scaling: y = alpha * A * x + beta * y
    void Multiply(Real alpha, const Vector& x, Real beta, Vector& y) const override {
        for (Integer i = 0; i < nrows_; ++i) {
            Real sum = 0.0;
            for (Integer k = rows_[i]; k < rows_[i + 1]; ++k) {
                sum += values_[k] * x[cols_[k]];
            }
            y[i] = alpha * sum + beta * y[i];
        }
    }
    
    // Accessors for CRS data (needed for Eigen conversion)
    const std::vector<Integer>& GetRowPointers() const { return rows_; }
    const std::vector<Integer>& GetColumnIndices() const { return cols_; }
    const std::vector<Real>& GetValues() const { return values_; }
    
    // Sort matrix (corresponds to CRS_SortMatrix)
    void SortRows() {
        for (Integer i = 0; i < nrows_; ++i) {
            Integer start = rows_[i];
            Integer end = rows_[i + 1];
            
            // Sort column indices for current row
            std::vector<Integer> indices(end - start);
            std::vector<Real> temp_values(end - start);
            
            for (Integer k = 0; k < end - start; ++k) {
                indices[k] = cols_[start + k];
                temp_values[k] = values_[start + k];
            }
            
            // Sort using column indices
            std::vector<Integer> sort_indices(end - start);
            for (Integer k = 0; k < end - start; ++k) sort_indices[k] = k;
            
            std::sort(sort_indices.begin(), sort_indices.end(), 
                     [&](Integer a, Integer b) { 
                         return indices[a] < indices[b]; 
                     });
            
            // Rearrange columns and values
            for (Integer k = 0; k < end - start; ++k) {
                cols_[start + k] = indices[sort_indices[k]];
                values_[start + k] = temp_values[sort_indices[k]];
            }
        }
    }
    
private:
    Integer nrows_, ncols_;
    std::vector<Integer> rows_;  // Row pointers
    std::vector<Integer> cols_;  // Column indices
    std::vector<Real> values_;   // Non-zero values
    
    // Find element position
    std::vector<Real>::iterator FindElementPosition(Integer i, Integer j) {
        // Check bounds
        if (i < 0 || i >= nrows_ || j < 0 || j >= ncols_) {
            return values_.end();
        }
        
        Integer start = rows_[i];
        Integer end = rows_[i + 1];
        
        // Check if start and end are valid
        if (start < 0 || end > cols_.size() || start > end) {
            return values_.end();
        }
        
        for (Integer k = start; k < end; ++k) {
            if (cols_[k] == j) {
                return values_.begin() + k;
            }
        }
        return values_.end();
    }
    
    std::vector<Real>::const_iterator FindElementPosition(Integer i, Integer j) const {
        // Check bounds
        if (i < 0 || i >= nrows_ || j < 0 || j >= ncols_) {
            return values_.end();
        }
        
        Integer start = rows_[i];
        Integer end = rows_[i + 1];
        
        // Check if start and end are valid
        if (start < 0 || end > cols_.size() || start > end) {
            return values_.end();
        }
        
        for (Integer k = start; k < end; ++k) {
            if (cols_[k] == j) {
                return values_.begin() + k;
            }
        }
        return values_.end();
    }
    
    // Insert new element
    void InsertElement(Integer i, Integer j, Real value) {
        // Check bounds
        if (i < 0 || i >= nrows_ || j < 0 || j >= ncols_) {
            return;
        }
        
        Integer start = rows_[i];
        Integer end = rows_[i + 1];
        
        // Find insertion position
        Integer insert_pos = start;
        while (insert_pos < end && cols_[insert_pos] < j) {
            ++insert_pos;
        }
        
        // Insert new element
        cols_.insert(cols_.begin() + insert_pos, j);
        values_.insert(values_.begin() + insert_pos, value);
        
        // Update row pointers
        for (Integer k = i + 1; k <= nrows_; ++k) {
            if (k < rows_.size()) {
                rows_[k]++;
            }
        }
    }
};

// Band matrix implementation (corresponds to BandMatrix.F90)
class BandMatrix : public Matrix {
public:
    BandMatrix(Integer nrows, Integer bandwidth) 
        : nrows_(nrows), bandwidth_(bandwidth) {
        data_.resize(nrows * (2 * bandwidth + 1), 0.0);
    }
    
    ~BandMatrix() override = default;
    
    Integer GetNumRows() const override { return nrows_; }
    Integer GetNumCols() const override { return nrows_; } // Band matrix is square
    
    void Zero() override {
        std::fill(data_.begin(), data_.end(), 0.0);
    }
    
    void SetElement(Integer i, Integer j, Real value) override {
        if (std::abs(i - j) <= bandwidth_) {
            data_[GetIndex(i, j)] = value;
        }
    }
    
    Real GetElement(Integer i, Integer j) const override {
        if (std::abs(i - j) <= bandwidth_) {
            return data_[GetIndex(i, j)];
        }
        return 0.0;
    }
    
    void AddToElement(Integer i, Integer j, Real value) override {
        if (std::abs(i - j) <= bandwidth_) {
            data_[GetIndex(i, j)] += value;
        }
    }
    
    // Matrix-vector multiplication: y = A * x
    void Multiply(const Vector& x, Vector& y) const override {
        if (x.Size() != nrows_ || y.Size() != nrows_) {
            throw std::runtime_error("Matrix-vector multiplication dimension mismatch");
        }
        
        for (Integer i = 0; i < nrows_; ++i) {
            Real sum = 0.0;
            for (Integer j = std::max(0, i - bandwidth_); j <= std::min(nrows_ - 1, i + bandwidth_); ++j) {
                if (std::abs(i - j) <= bandwidth_) {
                    sum += data_[GetIndex(i, j)] * x[j];
                }
            }
            y[i] = sum;
        }
    }
    
    // Matrix-vector multiplication with scaling: y = alpha * A * x + beta * y
    void Multiply(Real alpha, const Vector& x, Real beta, Vector& y) const override {
        if (x.Size() != nrows_ || y.Size() != nrows_) {
            throw std::runtime_error("Matrix-vector multiplication dimension mismatch");
        }
        
        for (Integer i = 0; i < nrows_; ++i) {
            Real sum = 0.0;
            for (Integer j = std::max(0, i - bandwidth_); j <= std::min(nrows_ - 1, i + bandwidth_); ++j) {
                if (std::abs(i - j) <= bandwidth_) {
                    sum += data_[GetIndex(i, j)] * x[j];
                }
            }
            y[i] = alpha * sum + beta * y[i];
        }
    }
    
private:
    Integer nrows_, bandwidth_;
    std::vector<Real> data_;
    
    Integer GetIndex(Integer i, Integer j) const {
        return i * (2 * bandwidth_ + 1) + (j - i + bandwidth_);
    }
};

// Matrix factory methods
std::unique_ptr<Matrix> Matrix::CreateCRS(Integer nrows, Integer ncols) {
    return std::make_unique<CRSMatrix>(nrows, ncols);
}

std::unique_ptr<Matrix> Matrix::CreateBand(Integer nrows, Integer bandwidth) {
    return std::make_unique<BandMatrix>(nrows, bandwidth);
}

// Simple dense matrix implementation (for testing and small problems)
class DenseMatrix : public Matrix {
public:
    DenseMatrix(Integer nrows, Integer ncols) 
        : nrows_(nrows), ncols_(ncols) {
        data_.resize(nrows * ncols, 0.0);
    }
    
    ~DenseMatrix() override = default;
    
    Integer GetNumRows() const override { return nrows_; }
    Integer GetNumCols() const override { return ncols_; }
    
    void Zero() override {
        std::fill(data_.begin(), data_.end(), 0.0);
    }
    
    void SetElement(Integer i, Integer j, Real value) override {
        data_[i * ncols_ + j] = value;
    }
    
    Real GetElement(Integer i, Integer j) const override {
        return data_[i * ncols_ + j];
    }
    
    void AddToElement(Integer i, Integer j, Real value) override {
        data_[i * ncols_ + j] += value;
    }
    
    // Matrix-vector multiplication: y = A * x
    void Multiply(const Vector& x, Vector& y) const override {
        if (x.Size() != ncols_ || y.Size() != nrows_) {
            throw std::runtime_error("Matrix-vector multiplication dimension mismatch");
        }
        
        for (Integer i = 0; i < nrows_; ++i) {
            Real sum = 0.0;
            for (Integer j = 0; j < ncols_; ++j) {
                sum += data_[i * ncols_ + j] * x[j];
            }
            y[i] = sum;
        }
    }
    
    // Matrix-vector multiplication with scaling: y = alpha * A * x + beta * y
    void Multiply(Real alpha, const Vector& x, Real beta, Vector& y) const override {
        if (x.Size() != ncols_ || y.Size() != nrows_) {
            throw std::runtime_error("Matrix-vector multiplication dimension mismatch");
        }
        
        for (Integer i = 0; i < nrows_; ++i) {
            Real sum = 0.0;
            for (Integer j = 0; j < ncols_; ++j) {
                sum += data_[i * ncols_ + j] * x[j];
            }
            y[i] = alpha * sum + beta * y[i];
        }
    }
    
private:
    Integer nrows_, ncols_;
    std::vector<Real> data_;
};

std::unique_ptr<Matrix> Matrix::CreateDense(Integer nrows, Integer ncols) {
    return std::make_unique<DenseMatrix>(nrows, ncols);
}

// Utility function to convert to Eigen sparse matrix
Eigen::SparseMatrix<Real> ConvertToEigen(const Matrix& matrix) {
    // For simplicity, assume it's a CRSMatrix
    const CRSMatrix* crs_matrix = dynamic_cast<const CRSMatrix*>(&matrix);
    if (!crs_matrix) {
        throw std::runtime_error("Only CRSMatrix can be converted to Eigen");
    }
    
    // Extract CRS data and build Eigen sparse matrix
    // Get the internal data from CRSMatrix using public accessors
    const std::vector<Integer>& rows = crs_matrix->GetRowPointers();
    const std::vector<Integer>& cols = crs_matrix->GetColumnIndices();
    const std::vector<Real>& values = crs_matrix->GetValues();
    
    Integer nrows = crs_matrix->GetNumRows();
    Integer ncols = crs_matrix->GetNumCols();
    
    // Count total non-zero elements
    Integer nnz = values.size();
    
    // Create Eigen sparse matrix
    Eigen::SparseMatrix<Real> eigen_matrix(nrows, ncols);
    eigen_matrix.reserve(nnz);
    
    // Fill the Eigen sparse matrix using triplets
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(nnz);
    
    for (Integer i = 0; i < nrows; ++i) {
        Integer start = rows[i];
        Integer end = rows[i + 1];
        
        for (Integer k = start; k < end; ++k) {
            triplets.emplace_back(i, cols[k], values[k]);
        }
    }
    
    eigen_matrix.setFromTriplets(triplets.begin(), triplets.end());
    
    return eigen_matrix;
}

} // namespace elmer