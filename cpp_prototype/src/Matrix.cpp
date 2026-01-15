// Matrix.cpp - Linear algebra matrix implementation
// Corresponds to Elmer FEM's CRSMatrix.F90 functionality

#include "Types.h"
#include "CRSMatrix.h"
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