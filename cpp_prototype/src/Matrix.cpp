// Matrix.cpp - 线性代数矩阵实现
// 对应Elmer FEM的CRSMatrix.F90功能

#include "ElmerCpp.h"
#include <Eigen/Sparse>
#include <vector>
#include <algorithm>
#include <memory>

namespace elmer {

// CRS矩阵实现（对应CRSMatrix.F90）
class CRSMatrix : public Matrix {
public:
    CRSMatrix(Integer nrows, Integer ncols) 
        : nrows_(nrows), ncols_(ncols) {
        rows_.resize(nrows + 1, 0);
    }
    
    ~CRSMatrix() override = default;
    
    void Zero() override {
        std::fill(values_.begin(), values_.end(), 0.0);
    }
    
    void SetElement(Integer i, Integer j, Real value) override {
        // 查找或插入元素
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
    
    // 获取矩阵-向量乘积（对应Fortran中的矩阵向量操作）
    void Multiply(const Vector& x, Vector& y) const {
        for (Integer i = 0; i < nrows_; ++i) {
            Real sum = 0.0;
            for (Integer k = rows_[i]; k < rows_[i + 1]; ++k) {
                sum += values_[k] * x[cols_[k]];
            }
            y[i] = sum;
        }
    }
    
    // 排序矩阵（对应CRS_SortMatrix）
    void SortRows() {
        for (Integer i = 0; i < nrows_; ++i) {
            Integer start = rows_[i];
            Integer end = rows_[i + 1];
            
            // 对当前行的列索引进行排序
            std::vector<Integer> indices(end - start);
            std::vector<Real> temp_values(end - start);
            
            for (Integer k = 0; k < end - start; ++k) {
                indices[k] = cols_[start + k];
                temp_values[k] = values_[start + k];
            }
            
            // 使用列索引排序
            std::vector<Integer> sort_indices(end - start);
            for (Integer k = 0; k < end - start; ++k) sort_indices[k] = k;
            
            std::sort(sort_indices.begin(), sort_indices.end(), 
                     [&](Integer a, Integer b) { 
                         return indices[a] < indices[b]; 
                     });
            
            // 重新排列列和值
            for (Integer k = 0; k < end - start; ++k) {
                cols_[start + k] = indices[sort_indices[k]];
                values_[start + k] = temp_values[sort_indices[k]];
            }
        }
    }
    
private:
    Integer nrows_, ncols_;
    std::vector<Integer> rows_;  // 行指针
    std::vector<Integer> cols_;  // 列索引
    std::vector<Real> values_;   // 非零值
    
    // 查找元素位置
    std::vector<Real>::iterator FindElementPosition(Integer i, Integer j) {
        Integer start = rows_[i];
        Integer end = rows_[i + 1];
        
        for (Integer k = start; k < end; ++k) {
            if (cols_[k] == j) {
                return values_.begin() + k;
            }
        }
        return values_.end();
    }
    
    std::vector<Real>::const_iterator FindElementPosition(Integer i, Integer j) const {
        Integer start = rows_[i];
        Integer end = rows_[i + 1];
        
        for (Integer k = start; k < end; ++k) {
            if (cols_[k] == j) {
                return values_.begin() + k;
            }
        }
        return values_.end();
    }
    
    // 插入新元素
    void InsertElement(Integer i, Integer j, Real value) {
        Integer start = rows_[i];
        Integer end = rows_[i + 1];
        
        // 找到插入位置
        Integer insert_pos = start;
        while (insert_pos < end && cols_[insert_pos] < j) {
            ++insert_pos;
        }
        
        // 插入新元素
        cols_.insert(cols_.begin() + insert_pos, j);
        values_.insert(values_.begin() + insert_pos, value);
        
        // 更新行指针
        for (Integer k = i + 1; k <= nrows_; ++k) {
            rows_[k]++;
        }
    }
};

// 带状矩阵实现（对应BandMatrix.F90）
class BandMatrix : public Matrix {
public:
    BandMatrix(Integer nrows, Integer bandwidth) 
        : nrows_(nrows), bandwidth_(bandwidth) {
        data_.resize(nrows * (2 * bandwidth + 1), 0.0);
    }
    
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
    
private:
    Integer nrows_, bandwidth_;
    std::vector<Real> data_;
    
    Integer GetIndex(Integer i, Integer j) const {
        return i * (2 * bandwidth_ + 1) + (j - i + bandwidth_);
    }
};

// 工厂方法实现
std::unique_ptr<Matrix> Matrix::CreateCRS(Integer nrows, Integer ncols) {
    return std::make_unique<CRSMatrix>(nrows, ncols);
}

std::unique_ptr<Matrix> Matrix::CreateBand(Integer nrows, Integer bandwidth) {
    return std::make_unique<BandMatrix>(nrows, bandwidth);
}

// Eigen矩阵包装器（可选，用于高性能计算）
class EigenMatrixWrapper : public Matrix {
public:
    EigenMatrixWrapper(Integer nrows, Integer ncols) 
        : matrix_(nrows, ncols) {}
    
    void Zero() override {
        matrix_.setZero();
    }
    
    void SetElement(Integer i, Integer j, Real value) override {
        matrix_.coeffRef(i, j) = value;
    }
    
    Real GetElement(Integer i, Integer j) const override {
        return matrix_.coeff(i, j);
    }
    
    void AddToElement(Integer i, Integer j, Real value) override {
        matrix_.coeffRef(i, j) += value;
    }
    
    // 转换为Eigen稀疏矩阵
    const Eigen::SparseMatrix<Real>& GetEigenMatrix() const { 
        return matrix_; 
    }
    
private:
    Eigen::SparseMatrix<Real> matrix_;
};

} // namespace elmer