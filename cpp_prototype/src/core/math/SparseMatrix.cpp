// SparseMatrix.cpp - 稀疏矩阵实现文件
// 实现CSR、CSC等稀疏矩阵格式的具体功能

#include "SparseMatrix.h"
#include "Types.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <iostream>

namespace elmer {

// ============================================================================
// CSRMatrix 实现
// ============================================================================

Real CSRMatrix::GetElement(Integer i, Integer j) const {
    CheckIndices(i, j);
    
    if (!is_assembled_) {
        // 在装配前，从临时存储中查找
        for (const auto& elem : temp_storage_[i]) {
            if (elem.first == j) {
                return elem.second;
            }
        }
        return 0.0;
    }
    
    // 在CSR存储中查找
    Integer start = row_pointers_[i];
    Integer end = row_pointers_[i + 1];
    
    for (Integer k = start; k < end; ++k) {
        if (col_indices_[k] == j) {
            return values_[k];
        }
    }
    
    return 0.0;
}

void CSRMatrix::SetElement(Integer i, Integer j, Real value) {
    CheckIndices(i, j);
    
    if (is_assembled_) {
        // 如果已经装配，需要重新构建
        throw std::runtime_error("矩阵已装配，不能直接设置元素。请使用AddToElement或在装配前设置。");
    }
    
    // 在临时存储中查找并设置
    bool found = false;
    for (auto& elem : temp_storage_[i]) {
        if (elem.first == j) {
            elem.second = value;
            found = true;
            break;
        }
    }
    
    if (!found) {
        temp_storage_[i].emplace_back(j, value);
    }
}

void CSRMatrix::AddToElement(Integer i, Integer j, Real value) {
    CheckIndices(i, j);
    
    if (is_assembled_) {
        // 在已装配的矩阵中累加
        Integer pos = FindElementPosition(i, j);
        if (pos >= 0) {
            values_[pos] += value;
        } else {
            throw std::runtime_error("无法在已装配的矩阵中添加新元素");
        }
    } else {
        // 在临时存储中累加
        bool found = false;
        for (auto& elem : temp_storage_[i]) {
            if (elem.first == j) {
                elem.second += value;
                found = true;
                break;
            }
        }
        
        if (!found) {
            temp_storage_[i].emplace_back(j, value);
        }
    }
}

void CSRMatrix::Zero() {
    if (is_assembled_) {
        std::fill(values_.begin(), values_.end(), 0.0);
        if (diagonal_computed_) {
            std::fill(diagonal_.begin(), diagonal_.end(), 0.0);
        }
    } else {
        for (auto& row : temp_storage_) {
            row.clear();
        }
    }
}

void CSRMatrix::Multiply(const Vector& x, Vector& y) const {
    if (!is_assembled_) {
        throw std::runtime_error("矩阵未装配，不能进行矩阵向量乘法");
    }
    
    if (x.Size() != ncols_ || y.Size() != nrows_) {
        throw std::invalid_argument("向量尺寸与矩阵不匹配");
    }
    
    // 清零输出向量
    y.Zero();
    
    // CSR格式的矩阵向量乘法
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        Real sum = 0.0;
        
        for (Integer k = start; k < end; ++k) {
            sum += values_[k] * x[col_indices_[k]];
        }
        
        y[i] = sum;
    }
}

void CSRMatrix::MultiplyTranspose(const Vector& x, Vector& y) const {
    if (!is_assembled_) {
        throw std::runtime_error("矩阵未装配，不能进行转置矩阵向量乘法");
    }
    
    if (x.Size() != nrows_ || y.Size() != ncols_) {
        throw std::invalid_argument("向量尺寸与矩阵不匹配");
    }
    
    // 清零输出向量
    y.Zero();
    
    // CSR格式的转置矩阵向量乘法
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        for (Integer k = start; k < end; ++k) {
            Integer j = col_indices_[k];
            y[j] += values_[k] * x[i];
        }
    }
}

void CSRMatrix::Assemble() {
    if (is_assembled_) {
        return; // 已经装配
    }
    
    // 计算每行的非零元素数量
    std::vector<Integer> row_nnz(nrows_, 0);
    for (Integer i = 0; i < nrows_; ++i) {
        row_nnz[i] = static_cast<Integer>(temp_storage_[i].size());
    }
    
    // 构建行指针
    row_pointers_[0] = 0;
    for (Integer i = 0; i < nrows_; ++i) {
        row_pointers_[i + 1] = row_pointers_[i] + row_nnz[i];
    }
    
    // 分配存储空间
    Integer total_nnz = row_pointers_[nrows_];
    values_.resize(total_nnz);
    col_indices_.resize(total_nnz);
    
    // 填充数据
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        
        // 对每行的元素按列索引排序（优化缓存性能）
        std::sort(temp_storage_[i].begin(), temp_storage_[i].end(),
                 [](const auto& a, const auto& b) { return a.first < b.first; });
        
        for (Integer k = 0; k < row_nnz[i]; ++k) {
            Integer pos = start + k;
            col_indices_[pos] = temp_storage_[i][k].first;
            values_[pos] = temp_storage_[i][k].second;
        }
    }
    
    // 清空临时存储
    for (auto& row : temp_storage_) {
        row.clear();
    }
    
    is_assembled_ = true;
    ComputeDiagonal();
}

void CSRMatrix::OptimizeStorage() {
    if (!is_assembled_) {
        throw std::runtime_error("矩阵未装配，不能优化存储");
    }
    
    // 压缩存储：移除零元素
    std::vector<Real> new_values;
    std::vector<Integer> new_col_indices;
    std::vector<Integer> new_row_pointers(nrows_ + 1, 0);
    
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        for (Integer k = start; k < end; ++k) {
            if (std::abs(values_[k]) > 1e-15) { // 忽略很小的值
                new_values.push_back(values_[k]);
                new_col_indices.push_back(col_indices_[k]);
                new_row_pointers[i + 1]++;
            }
        }
    }
    
    // 累积行指针
    for (Integer i = 0; i < nrows_; ++i) {
        new_row_pointers[i + 1] += new_row_pointers[i];
    }
    
    // 更新数据
    values_ = std::move(new_values);
    col_indices_ = std::move(new_col_indices);
    row_pointers_ = std::move(new_row_pointers);
    
    ComputeDiagonal();
}

std::size_t CSRMatrix::GetMemoryUsage() const {
    std::size_t usage = 0;
    
    usage += values_.capacity() * sizeof(Real);
    usage += col_indices_.capacity() * sizeof(Integer);
    usage += row_pointers_.capacity() * sizeof(Integer);
    usage += diagonal_.capacity() * sizeof(Real);
    
    if (!is_assembled_) {
        for (const auto& row : temp_storage_) {
            usage += row.capacity() * sizeof(std::pair<Integer, Real>);
        }
    }
    
    return usage;
}

std::unique_ptr<SparseMatrix> CSRMatrix::ConvertTo(SparseMatrixFormat new_format) const {
    if (!is_assembled_) {
        throw std::runtime_error("矩阵未装配，不能转换格式");
    }
    
    switch (new_format) {
        case SparseMatrixFormat::CSR:
            // 已经是CSR格式，返回副本
            return std::make_unique<CSRMatrix>(nrows_, ncols_, values_, col_indices_, row_pointers_);
            
        case SparseMatrixFormat::CSC:
            // 转换为CSC格式
            // TODO: 实现CSR到CSC的转换
            throw std::runtime_error("CSR到CSC转换尚未实现");
            
        case SparseMatrixFormat::COO:
            // 转换为COO格式
            // TODO: 实现CSR到COO的转换
            throw std::runtime_error("CSR到COO转换尚未实现");
            
        default:
            throw std::invalid_argument("不支持的矩阵格式");
    }
}

bool CSRMatrix::CheckSymmetry(Real tolerance) const {
    if (!is_assembled_) {
        throw std::runtime_error("矩阵未装配，不能检查对称性");
    }
    
    if (nrows_ != ncols_) {
        return false; // 非方阵不可能对称
    }
    
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        for (Integer k = start; k < end; ++k) {
            Integer j = col_indices_[k];
            Real a_ij = values_[k];
            Real a_ji = GetElement(j, i);
            
            if (std::abs(a_ij - a_ji) > tolerance) {
                return false;
            }
        }
    }
    
    return true;
}

Real CSRMatrix::GetDiagonal(Integer i) const {
    CheckIndices(i, i);
    
    if (diagonal_computed_) {
        return diagonal_[i];
    }
    
    // 如果没有预计算，直接查找
    return GetElement(i, i);
}

void CSRMatrix::SetDiagonal(Integer i, Real value) {
    CheckIndices(i, i);
    
    if (is_assembled_) {
        Integer pos = FindElementPosition(i, i);
        if (pos >= 0) {
            values_[pos] = value;
            if (diagonal_computed_) {
                diagonal_[i] = value;
            }
        } else {
            throw std::runtime_error("无法设置不存在的对角线元素");
        }
    } else {
        SetElement(i, i, value);
    }
}

bool CSRMatrix::SaveToFile(const std::string& filename) const {
    if (!is_assembled_) {
        throw std::runtime_error("矩阵未装配，不能保存到文件");
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    // 保存矩阵基本信息
    file << "%%MatrixMarket matrix coordinate real general" << std::endl;
    file << nrows_ << " " << ncols_ << " " << GetNonZeroCount() << std::endl;
    
    // 保存非零元素
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        for (Integer k = start; k < end; ++k) {
            file << (i + 1) << " " << (col_indices_[k] + 1) << " " 
                 << std::scientific << std::setprecision(15) << values_[k] << std::endl;
        }
    }
    
    file.close();
    return true;
}

bool CSRMatrix::LoadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    // 跳过注释行
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '%') {
            continue;
        }
        break;
    }
    
    // 读取矩阵尺寸
    Integer rows, cols, nnz;
    std::istringstream iss(line);
    iss >> rows >> cols >> nnz;
    
    if (rows <= 0 || cols <= 0 || nnz < 0) {
        return false;
    }
    
    // 重置矩阵
    nrows_ = rows;
    ncols_ = cols;
    row_pointers_.resize(nrows_ + 1, 0);
    temp_storage_.resize(nrows_);
    
    // 读取非零元素
    for (Integer k = 0; k < nnz; ++k) {
        Integer i, j;
        Real value;
        
        if (!std::getline(file, line)) {
            return false;
        }
        
        std::istringstream elem_iss(line);
        elem_iss >> i >> j >> value;
        
        // MatrixMarket使用1-based索引，转换为0-based
        SetElement(i - 1, j - 1, value);
    }
    
    file.close();
    
    // 装配矩阵
    Assemble();
    return true;
}

void CSRMatrix::PrintStructure() const {
    if (!is_assembled_) {
        std::cout << "矩阵未装配" << std::endl;
        return;
    }
    
    std::cout << "CSR矩阵结构 (" << nrows_ << "x" << ncols_ 
              << ", 非零元素: " << GetNonZeroCount() << ")" << std::endl;
    
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        std::cout << "行 " << i << ": ";
        for (Integer k = start; k < end; ++k) {
            std::cout << "(" << col_indices_[k] << ", " << values_[k] << ") ";
        }
        std::cout << std::endl;
    }
}

bool CSRMatrix::Validate() const {
    if (!is_assembled_) {
        return true; // 未装配的矩阵总是有效的
    }
    
    // 检查行指针
    if (row_pointers_.size() != static_cast<std::size_t>(nrows_ + 1)) {
        return false;
    }
    
    if (row_pointers_[0] != 0 || row_pointers_[nrows_] != GetNonZeroCount()) {
        return false;
    }
    
    // 检查列索引
    for (Integer j : col_indices_) {
        if (j < 0 || j >= ncols_) {
            return false;
        }
    }
    
    // 检查行指针单调性
    for (Integer i = 0; i < nrows_; ++i) {
        if (row_pointers_[i] > row_pointers_[i + 1]) {
            return false;
        }
    }
    
    return true;
}

std::vector<Integer> CSRMatrix::GetRowColumnIndices(Integer row) const {
    CheckIndices(row, 0);
    
    if (!is_assembled_) {
        throw std::runtime_error("矩阵未装配");
    }
    
    Integer start = row_pointers_[row];
    Integer end = row_pointers_[row + 1];
    
    return std::vector<Integer>(col_indices_.begin() + start, 
                                col_indices_.begin() + end);
}

std::vector<Real> CSRMatrix::GetRowValues(Integer row) const {
    CheckIndices(row, 0);
    
    if (!is_assembled_) {
        throw std::runtime_error("矩阵未装配");
    }
    
    Integer start = row_pointers_[row];
    Integer end = row_pointers_[row + 1];
    
    return std::vector<Real>(values_.begin() + start, 
                            values_.begin() + end);
}

Integer CSRMatrix::GetRowNonZeroCount(Integer row) const {
    CheckIndices(row, 0);
    
    if (!is_assembled_) {
        return static_cast<Integer>(temp_storage_[row].size());
    }
    
    return row_pointers_[row + 1] - row_pointers_[row];
}

void CSRMatrix::ComputeDiagonal() {
    diagonal_.resize(nrows_, 0.0);
    
    for (Integer i = 0; i < nrows_; ++i) {
        Integer start = row_pointers_[i];
        Integer end = row_pointers_[i + 1];
        
        for (Integer k = start; k < end; ++k) {
            if (col_indices_[k] == i) {
                diagonal_[i] = values_[k];
                break;
            }
        }
    }
    
    diagonal_computed_ = true;
}

Integer CSRMatrix::FindElementPosition(Integer i, Integer j) const {
    if (!is_assembled_) {
        return -1;
    }
    
    Integer start = row_pointers_[i];
    Integer end = row_pointers_[i + 1];
    
    // 二分查找（因为列索引已排序）
    auto it = std::lower_bound(col_indices_.begin() + start, 
                              col_indices_.begin() + end, j);
    
    if (it != col_indices_.begin() + end && *it == j) {
        return static_cast<Integer>(it - col_indices_.begin());
    }
    
    return -1;
}

void CSRMatrix::InsertElement(Integer i, Integer j, Real value, Integer pos) {
    // 在已装配的矩阵中插入元素（通常不推荐使用）
    values_.insert(values_.begin() + pos, value);
    col_indices_.insert(col_indices_.begin() + pos, j);
    
    // 更新行指针
    for (Integer k = i + 1; k <= nrows_; ++k) {
        row_pointers_[k]++;
    }
}

void CSRMatrix::CheckIndices(Integer i, Integer j) const {
    if (i < 0 || i >= nrows_ || j < 0 || j >= ncols_) {
        throw std::out_of_range("矩阵索引越界");
    }
}

// ============================================================================
// SparseMatrixUtils 实现
// ============================================================================

std::unique_ptr<CSRMatrix> SparseMatrixUtils::CreateDiagonallyDominantMatrix(Integer size, Real diagonal_value) {
    auto matrix = std::make_unique<CSRMatrix>(size, size);
    
    for (Integer i = 0; i < size; ++i) {
        // 对角线元素
        matrix->SetElement(i, i, diagonal_value);
        
        // 邻接元素
        if (i > 0) {
            matrix->SetElement(i, i - 1, -1.0);
        }
        if (i < size - 1) {
            matrix->SetElement(i, i + 1, -1.0);
        }
    }
    
    matrix->Assemble();
    return matrix;
}

std::unique_ptr<CSRMatrix> SparseMatrixUtils::CreateLaplacianMatrix(Integer size) {
    auto matrix = std::make_unique<CSRMatrix>(size, size);
    
    for (Integer i = 0; i < size; ++i) {
        // 对角线元素：度
        Integer degree = 0;
        
        if (i > 0) {
            matrix->SetElement(i, i - 1, -1.0);
            degree++;
        }
        if (i < size - 1) {
            matrix->SetElement(i, i + 1, -1.0);
            degree++;
        }
        
        matrix->SetElement(i, i, static_cast<Real>(degree));
    }
    
    matrix->Assemble();
    return matrix;
}

// 其他工具函数的实现...

} // namespace elmer