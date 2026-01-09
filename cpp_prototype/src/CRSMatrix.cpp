// CRSMatrix.cpp - Compressed Row Storage Matrix Implementation
// Corresponds to Elmer FEM's CRSMatrix.F90

#include "CRSMatrix.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace elmer {

CRSMatrix::CRSMatrix(Integer nrows, Integer ncols) 
    : nrows_(nrows), ncols_(ncols), diagonal_computed_(false) {
    
    if (nrows <= 0 || ncols <= 0) {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }
    
    // Initialize row pointers (all rows start at position 0 initially)
    row_pointers_.resize(nrows + 1, 0);
    
    // Pre-allocate some space for efficiency
    values_.reserve(nrows * 10);  // Estimate 10 non-zeros per row
    col_indices_.reserve(nrows * 10);
    
    // Initialize diagonal storage
    diagonal_.resize(nrows, 0.0);
}

void CRSMatrix::Zero() {
    std::fill(values_.begin(), values_.end(), 0.0);
    diagonal_computed_ = false;
    std::fill(diagonal_.begin(), diagonal_.end(), 0.0);
}

Real CRSMatrix::GetElement(Integer i, Integer j) const {
    CheckIndices(i, j);
    
    Integer pos = FindElementPosition(i, j);
    if (pos >= 0) {
        return values_[pos];
    }
    return 0.0;
}

void CRSMatrix::SetElement(Integer i, Integer j, Real value) {
    CheckIndices(i, j);
    
    Integer pos = FindElementPosition(i, j);
    if (pos >= 0) {
        // Element exists, update value
        values_[pos] = value;
    } else {
        // Element doesn't exist, need to insert
        Integer row_start = row_pointers_[i];
        Integer row_end = row_pointers_[i + 1];
        
        // Find insertion position (keep columns sorted)
        Integer insert_pos = row_start;
        while (insert_pos < row_end && col_indices_[insert_pos] < j) {
            ++insert_pos;
        }
        
        // Insert new element
        values_.insert(values_.begin() + insert_pos, value);
        col_indices_.insert(col_indices_.begin() + insert_pos, j);
        
        // Update row pointers for all subsequent rows
        for (Integer k = i + 1; k <= nrows_; ++k) {
            row_pointers_[k]++;
        }
    }
    
    // Update diagonal if needed
    if (i == j) {
        diagonal_[i] = value;
        diagonal_computed_ = true;
    }
}

void CRSMatrix::AddToElement(Integer i, Integer j, Real value) {
    CheckIndices(i, j);
    
    Integer pos = FindElementPosition(i, j);
    if (pos >= 0) {
        values_[pos] += value;
    } else {
        SetElement(i, j, value);
        return;
    }
    
    // Update diagonal if needed
    if (i == j) {
        diagonal_[i] += value;
        diagonal_computed_ = true;
    }
}

Integer CRSMatrix::GetRowPointer(Integer i) const {
    if (i < 0 || i > nrows_) {
        throw std::out_of_range("Row index out of range");
    }
    return row_pointers_[i];
}

std::vector<Integer> CRSMatrix::GetColumnIndices(Integer i) const {
    CheckIndices(i, 0);  // Only check row index
    
    Integer start = row_pointers_[i];
    Integer end = row_pointers_[i + 1];
    
    std::vector<Integer> indices;
    indices.reserve(end - start);
    
    for (Integer j = start; j < end; ++j) {
        indices.push_back(col_indices_[j]);
    }
    
    return indices;
}

std::vector<Real> CRSMatrix::GetRowValues(Integer i) const {
    CheckIndices(i, 0);  // Only check row index
    
    Integer start = row_pointers_[i];
    Integer end = row_pointers_[i + 1];
    
    std::vector<Real> row_values;
    row_values.reserve(end - start);
    
    for (Integer j = start; j < end; ++j) {
        row_values.push_back(values_[j]);
    }
    
    return row_values;
}

void CRSMatrix::Multiply(const Vector& x, Vector& y) const {
    if (x.Size() != ncols_ || y.Size() != nrows_) {
        throw std::invalid_argument("Vector dimensions don't match matrix");
    }
    
    y.Zero();
    
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        Real sum = 0.0;
        for (Integer j = start; j < end; ++j) {
            sum += values_[j] * x[col_indices_[j]];
        }
        y[i] = sum;
    }
}

void CRSMatrix::Multiply(Real alpha, const Vector& x, Real beta, Vector& y) const {
    if (x.Size() != ncols_ || y.Size() != nrows_) {
        throw std::invalid_argument("Vector dimensions don't match matrix");
    }
    
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        Real sum = 0.0;
        for (Integer j = start; j < end; ++j) {
            sum += values_[j] * x[col_indices_[j]];
        }
        y[i] = alpha * sum + beta * y[i];
    }
}

void CRSMatrix::MultiplyTranspose(const Vector& x, Vector& y) const {
    if (x.Size() != nrows_ || y.Size() != ncols_) {
        throw std::invalid_argument("Vector dimensions don't match matrix transpose");
    }
    
    y.Zero();
    
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        Real xi = x[i];
        for (Integer j = start; j < end; ++j) {
            Integer col = col_indices_[j];
            y[col] += values_[j] * xi;
        }
    }
}

Real CRSMatrix::GetDiagonal(Integer i) const {
    CheckIndices(i, i);
    
    if (!diagonal_computed_) {
        // Compute diagonal on-the-fly
        const_cast<CRSMatrix*>(this)->ComputeDiagonal();
    }
    
    return diagonal_[i];
}

void CRSMatrix::ComputeDiagonal() {
    std::fill(diagonal_.begin(), diagonal_.end(), 0.0);
    
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        for (Integer j = start; j < end; ++j) {
            if (col_indices_[j] == i) {
                diagonal_[i] = values_[j];
                break;
            }
        }
    }
    
    diagonal_computed_ = true;
}

bool CRSMatrix::IsSymmetric(Real tolerance) const {
    if (nrows_ != ncols_) return false;
    
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        for (Integer j = start; j < end; ++j) {
            Integer col = col_indices_[j];
            Real a_ij = values_[j];
            
            // Find symmetric element
            Integer sym_pos = FindElementPosition(col, i);
            if (sym_pos < 0) {
                // Symmetric element doesn't exist
                if (std::abs(a_ij) > tolerance) return false;
            } else {
                Real a_ji = values_[sym_pos];
                if (std::abs(a_ij - a_ji) > tolerance) return false;
            }
        }
    }
    
    return true;
}

std::size_t CRSMatrix::GetMemoryUsage() const {
    return sizeof(*this) + 
           values_.size() * sizeof(Real) + 
           col_indices_.size() * sizeof(Integer) + 
           row_pointers_.size() * sizeof(Integer) + 
           diagonal_.size() * sizeof(Real);
}

void CRSMatrix::PrintStructure() const {
    std::cout << "CRS Matrix Structure:" << std::endl;
    std::cout << "Rows: " << nrows_ << ", Columns: " << ncols_ << std::endl;
    std::cout << "Non-zero elements: " << values_.size() << std::endl;
    
    for (Integer i = 0; i < nrows_; ++i) {
        std::cout << "Row " << i << ": ";
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        for (Integer j = start; j < end; ++j) {
            std::cout << "(" << col_indices_[j] << ", " << values_[j] << ") ";
        }
        std::cout << std::endl;
    }
}

Integer CRSMatrix::FindElementPosition(Integer i, Integer j) const {
    CheckIndices(i, j);
    
    Integer start = row_pointers_[i];
    Integer end = row_pointers_[i + 1];
    
    // Binary search for column j (columns are sorted)
    Integer left = start;
    Integer right = end - 1;
    
    while (left <= right) {
        Integer mid = left + (right - left) / 2;
        Integer col = col_indices_[mid];
        
        if (col == j) {
            return mid;
        } else if (col < j) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    
    return -1;  // Not found
}

void CRSMatrix::CheckIndices(Integer i, Integer j) const {
    if (i < 0 || i >= nrows_) {
        throw std::out_of_range("Row index out of range");
    }
    if (j < 0 || j >= ncols_) {
        throw std::out_of_range("Column index out of range");
    }
}

std::unique_ptr<Matrix> elmer::CreateCRSMatrix(Integer nrows, Integer ncols) {
    return std::make_unique<CRSMatrix>(nrows, ncols);
}

} // namespace elmer