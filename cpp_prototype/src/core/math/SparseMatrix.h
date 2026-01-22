// SparseMatrix.h - 稀疏矩阵基础存储类
// 支持CSR和CSC格式，专为低频电磁有限元求解优化

#pragma once

#include "Types.h"
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <cstring>

namespace elmer {

/**
 * @brief 稀疏矩阵存储格式枚举
 */
enum class SparseMatrixFormat {
    CSR,    ///< 压缩行存储格式 (Compressed Row Storage)
    CSC,    ///< 压缩列存储格式 (Compressed Column Storage)
    COO     ///< 坐标格式 (Coordinate Format)
};

/**
 * @brief 稀疏矩阵对称性枚举
 */
enum class MatrixSymmetry {
    GENERAL,        ///< 一般矩阵
    SYMMETRIC,      ///< 对称矩阵
    HERMITIAN,      ///< 埃尔米特矩阵（复数对称）
    SKEW_SYMMETRIC  ///< 反对称矩阵
};

/**
 * @brief 稀疏矩阵基类
 * 
 * 提供稀疏矩阵的基本接口和通用功能，
 * 支持CSR、CSC等不同存储格式
 */
class SparseMatrix {
protected:
    Integer nrows_;                 ///< 矩阵行数
    Integer ncols_;                 ///< 矩阵列数
    SparseMatrixFormat format_;     ///< 存储格式
    MatrixSymmetry symmetry_;       ///< 对称性
    bool is_assembled_;             ///< 是否已完成装配
    
public:
    /**
     * @brief 构造函数
     * @param nrows 行数
     * @param ncols 列数
     * @param format 存储格式
     * @param symmetry 对称性
     */
    SparseMatrix(Integer nrows, Integer ncols, 
                 SparseMatrixFormat format = SparseMatrixFormat::CSR,
                 MatrixSymmetry symmetry = MatrixSymmetry::GENERAL)
        : nrows_(nrows), ncols_(ncols), 
          format_(format), symmetry_(symmetry), 
          is_assembled_(false) {}
    
    virtual ~SparseMatrix() = default;
    
    // 基本信息获取
    Integer GetNumRows() const { return nrows_; }
    Integer GetNumCols() const { return ncols_; }
    SparseMatrixFormat GetFormat() const { return format_; }
    MatrixSymmetry GetSymmetry() const { return symmetry_; }
    bool IsAssembled() const { return is_assembled_; }
    
    // 纯虚函数接口
    virtual Integer GetNonZeroCount() const = 0;
    virtual Real GetElement(Integer i, Integer j) const = 0;
    virtual void SetElement(Integer i, Integer j, Real value) = 0;
    virtual void AddToElement(Integer i, Integer j, Real value) = 0;
    virtual void Zero() = 0;
    
    // 矩阵-向量乘法
    virtual void Multiply(const Vector& x, Vector& y) const = 0;
    virtual void MultiplyTranspose(const Vector& x, Vector& y) const = 0;
    
    // 矩阵装配和优化
    virtual void Assemble() = 0;
    virtual void OptimizeStorage() = 0;
    
    // 内存使用统计
    virtual std::size_t GetMemoryUsage() const = 0;
    
    // 格式转换
    virtual std::unique_ptr<SparseMatrix> ConvertTo(SparseMatrixFormat new_format) const = 0;
    
    // 对称性检查
    virtual bool CheckSymmetry(Real tolerance = 1e-12) const = 0;
    
    // 对角线操作
    virtual Real GetDiagonal(Integer i) const = 0;
    virtual void SetDiagonal(Integer i, Real value) = 0;
    
    // 文件I/O
    virtual bool SaveToFile(const std::string& filename) const = 0;
    virtual bool LoadFromFile(const std::string& filename) = 0;
    
    // 调试和验证
    virtual void PrintStructure() const = 0;
    virtual bool Validate() const = 0;
};

/**
 * @brief CSR格式稀疏矩阵实现
 * 
 * 压缩行存储格式，特别适合有限元矩阵的存储和操作
 */
class CSRMatrix : public SparseMatrix {
private:
    std::vector<Real> values_;         ///< 非零元素值
    std::vector<Integer> col_indices_; ///< 列索引
    std::vector<Integer> row_pointers_; ///< 行指针
    
    // 对角线优化
    std::vector<Real> diagonal_;
    bool diagonal_computed_;
    
    // 临时存储用于装配
    std::vector<std::vector<std::pair<Integer, Real>>> temp_storage_;
    
public:
    /**
     * @brief 构造函数
     * @param nrows 行数
     * @param ncols 列数
     * @param symmetry 对称性
     */
    CSRMatrix(Integer nrows, Integer ncols, 
              MatrixSymmetry symmetry = MatrixSymmetry::GENERAL)
        : SparseMatrix(nrows, ncols, SparseMatrixFormat::CSR, symmetry),
          diagonal_computed_(false) {
        
        row_pointers_.resize(nrows_ + 1, 0);
        temp_storage_.resize(nrows_);
    }
    
    /**
     * @brief 从已有数据构造CSR矩阵
     * @param nrows 行数
     * @param ncols 列数
     * @param values 非零值数组
     * @param col_indices 列索引数组
     * @param row_pointers 行指针数组
     */
    CSRMatrix(Integer nrows, Integer ncols,
              const std::vector<Real>& values,
              const std::vector<Integer>& col_indices,
              const std::vector<Integer>& row_pointers)
        : SparseMatrix(nrows, ncols, SparseMatrixFormat::CSR, MatrixSymmetry::GENERAL),
          values_(values), col_indices_(col_indices), row_pointers_(row_pointers),
          diagonal_computed_(false) {
        
        if (row_pointers_.size() != static_cast<std::size_t>(nrows_ + 1)) {
            throw std::invalid_argument("行指针数组大小必须为nrows+1");
        }
        
        is_assembled_ = true;
        ComputeDiagonal();
    }
    
    // SparseMatrix接口实现
    Integer GetNonZeroCount() const override {
        return static_cast<Integer>(values_.size());
    }
    
    Real GetElement(Integer i, Integer j) const override;
    void SetElement(Integer i, Integer j, Real value) override;
    void AddToElement(Integer i, Integer j, Real value) override;
    void Zero() override;
    
    void Multiply(const Vector& x, Vector& y) const override;
    void MultiplyTranspose(const Vector& x, Vector& y) const override;
    
    void Assemble() override;
    void OptimizeStorage() override;
    
    std::size_t GetMemoryUsage() const override;
    
    std::unique_ptr<SparseMatrix> ConvertTo(SparseMatrixFormat new_format) const override;
    
    bool CheckSymmetry(Real tolerance = 1e-12) const override;
    
    Real GetDiagonal(Integer i) const override;
    void SetDiagonal(Integer i, Real value) override;
    
    bool SaveToFile(const std::string& filename) const override;
    bool LoadFromFile(const std::string& filename) override;
    
    void PrintStructure() const override;
    bool Validate() const override;
    
    // CSR特定方法
    /**
     * @brief 获取行指针
     */
    const std::vector<Integer>& GetRowPointers() const { return row_pointers_; }
    
    /**
     * @brief 获取列索引
     */
    const std::vector<Integer>& GetColumnIndices() const { return col_indices_; }
    
    /**
     * @brief 获取非零值
     */
    const std::vector<Real>& GetValues() const { return values_; }
    
    /**
     * @brief 获取指定行的列索引
     */
    std::vector<Integer> GetRowColumnIndices(Integer row) const;
    
    /**
     * @brief 获取指定行的非零值
     */
    std::vector<Real> GetRowValues(Integer row) const;
    
    /**
     * @brief 获取指定行的非零元素数量
     */
    Integer GetRowNonZeroCount(Integer row) const;
    
private:
    /**
     * @brief 计算对角线元素（优化访问）
     */
    void ComputeDiagonal();
    
    /**
     * @brief 查找元素位置
     */
    Integer FindElementPosition(Integer i, Integer j) const;
    
    /**
     * @brief 插入元素到指定位置
     */
    void InsertElement(Integer i, Integer j, Real value, Integer pos);
    
    /**
     * @brief 索引验证
     */
    void CheckIndices(Integer i, Integer j) const;
};

/**
 * @brief CSC格式稀疏矩阵实现
 * 
 * 压缩列存储格式，适合某些求解器的需求
 */
class CSCMatrix : public SparseMatrix {
private:
    std::vector<Real> values_;         ///< 非零元素值
    std::vector<Integer> row_indices_; ///< 行索引
    std::vector<Integer> col_pointers_; ///< 列指针
    
public:
    CSCMatrix(Integer nrows, Integer ncols, 
              MatrixSymmetry symmetry = MatrixSymmetry::GENERAL)
        : SparseMatrix(nrows, ncols, SparseMatrixFormat::CSC, symmetry) {
        
        col_pointers_.resize(ncols_ + 1, 0);
    }
    
    // SparseMatrix接口实现（类似CSRMatrix，但基于列存储）
    // 实现细节省略，结构与CSRMatrix类似
    
    // CSC特定方法
    const std::vector<Integer>& GetColPointers() const { return col_pointers_; }
    const std::vector<Integer>& GetRowIndices() const { return row_indices_; }
    const std::vector<Real>& GetValues() const { return values_; }
};

/**
 * @brief 稀疏矩阵工具类
 */
class SparseMatrixUtils {
public:
    /**
     * @brief 创建对角占优矩阵（用于测试）
     */
    static std::unique_ptr<CSRMatrix> CreateDiagonallyDominantMatrix(Integer size, Real diagonal_value = 2.0);
    
    /**
     * @brief 创建拉普拉斯矩阵（有限元离散常用）
     */
    static std::unique_ptr<CSRMatrix> CreateLaplacianMatrix(Integer size);
    
    /**
     * @brief 矩阵转置
     */
    static std::unique_ptr<SparseMatrix> Transpose(const SparseMatrix& matrix);
    
    /**
     * @brief 矩阵相加
     */
    static std::unique_ptr<SparseMatrix> Add(const SparseMatrix& A, const SparseMatrix& B, Real alpha = 1.0, Real beta = 1.0);
    
    /**
     * @brief 矩阵相乘
     */
    static std::unique_ptr<SparseMatrix> Multiply(const SparseMatrix& A, const SparseMatrix& B);
    
    /**
     * @brief 矩阵范数计算
     */
    static Real FrobeniusNorm(const SparseMatrix& matrix);
    static Real OneNorm(const SparseMatrix& matrix);
    static Real InfinityNorm(const SparseMatrix& matrix);
    
    /**
     * @brief 条件数估计
     */
    static Real EstimateConditionNumber(const SparseMatrix& matrix);
    
    /**
     * @brief 矩阵重排序（优化求解器性能）
     */
    static void ReorderMatrix(SparseMatrix& matrix, const std::vector<Integer>& permutation);
};

} // namespace elmer