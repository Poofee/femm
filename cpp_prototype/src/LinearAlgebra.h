// LinearAlgebra.h - Elmer FEM C++线性代数接口定义
// 对应Fortran模块: CRSMatrix.F90, IterSolve.F90, MatrixAssembly.F90

#pragma once

#include <memory>
#include <vector>
#include <string>
#include "ElmerCpp.h"

namespace elmer {

// ===== 求解器方法枚举 =====
enum class SolverMethod {
    CG,          // 共轭梯度法
    BiCGStab,    // 稳定双共轭梯度法
    GMRES,       // 广义最小残差法
    TFQMR,       // 无转置拟最小残差法
    MINRES       // 最小残差法
};

// ===== 求解器参数结构 =====
struct SolverParameters {
    SolverMethod method = SolverMethod::CG;
    double tolerance = 1e-8;
    int maxIterations = 1000;
    std::string preconditioner = "none"; // "none", "jacobi", "ilu"
    bool verbose = false;
};

// ===== 求解器结果结构 =====
struct SolverResult {
    int iterations;
    double residual;
    bool converged;
};



// ===== CRS矩阵实现 =====
class CRSMatrixImpl : public elmer::Matrix {
private:
    std::vector<double> values;
    std::vector<int> col_indices;
    std::vector<int> row_pointers;
    int rows, cols;
    
public:
    CRSMatrixImpl(int num_rows, int num_cols) : rows(num_rows), cols(num_cols) {
        row_pointers.resize(rows + 1, 0);
    }
    
    void SetElement(Integer i, Integer j, Real value) override {
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw std::out_of_range("Matrix index out of range");
        }
        
        // 查找插入位置
        int start = row_pointers[i];
        int end = row_pointers[i + 1];
        
        // 查找列索引
        auto it = std::lower_bound(col_indices.begin() + start, 
                                  col_indices.begin() + end, j);
        int pos = std::distance(col_indices.begin(), it);
        
        if (pos < end && col_indices[pos] == j) {
            // 更新现有元素
            values[pos] = value;
        } else {
            // 插入新元素
            values.insert(values.begin() + pos, value);
            col_indices.insert(col_indices.begin() + pos, j);
            
            // 更新行指针
            for (int k = i + 1; k <= rows; ++k) {
                row_pointers[k]++;
            }
        }
    }
    
    void AddToElement(Integer i, Integer j, Real value) override {
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw std::out_of_range("Matrix index out of range");
        }
        
        // 查找插入位置
        int start = row_pointers[i];
        int end = row_pointers[i + 1];
        
        // 查找列索引
        auto it = std::lower_bound(col_indices.begin() + start, 
                                  col_indices.begin() + end, j);
        int pos = std::distance(col_indices.begin(), it);
        
        if (pos < end && col_indices[pos] == j) {
            // 更新现有元素
            values[pos] += value;
        } else {
            // 插入新元素
            values.insert(values.begin() + pos, value);
            col_indices.insert(col_indices.begin() + pos, j);
            
            // 更新行指针
            for (int k = i + 1; k <= rows; ++k) {
                row_pointers[k]++;
            }
        }
    }
    
    Real GetElement(Integer i, Integer j) const override {
        if (i < 0 || i >= rows || j < 0 || j >= cols) {
            throw std::out_of_range("Matrix index out of range");
        }
        
        int start = row_pointers[i];
        int end = row_pointers[i + 1];
        
        // 二分查找列索引
        auto it = std::lower_bound(col_indices.begin() + start, 
                                  col_indices.begin() + end, j);
        
        if (it != col_indices.begin() + end && *it == j) {
            int pos = std::distance(col_indices.begin(), it);
            return values[pos];
        }
        
        return 0.0;
    }
    
    void Zero() override {
        values.clear();
        col_indices.clear();
        row_pointers.assign(rows + 1, 0);
    }
    
    Integer GetNumRows() const override {
        return rows;
    }
    
    Integer GetNumCols() const override {
        return cols;
    }
    
    void Multiply(const Vector& x, Vector& y) const override {
        // 实现矩阵-向量乘法
        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int j = row_pointers[i]; j < row_pointers[i + 1]; ++j) {
                sum += values[j] * x[col_indices[j]];
            }
            y[i] = sum;
        }
    }
    
    void Multiply(Real alpha, const Vector& x, Real beta, Vector& y) const override {
        // 实现带缩放的矩阵-向量乘法
        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int j = row_pointers[i]; j < row_pointers[i + 1]; ++j) {
                sum += values[j] * x[col_indices[j]];
            }
            y[i] = alpha * sum + beta * y[i];
        }
    }
};

// ===== 迭代求解器接口 =====
class IIterativeSolver {
public:
    virtual ~IIterativeSolver() = default;
    
    // 设置求解器参数
    virtual void setParameters(const SolverParameters& parameters) = 0;
    
    // 求解线性系统 Ax = b
    virtual SolverResult solve(const std::vector<double>& rhs) = 0;
    
    // 获取解向量
    virtual std::vector<double> getSolution() const = 0;
    
    // 获取残差历史
    virtual std::vector<double> getResidualHistory() const = 0;
};

// ===== 矩阵组装器接口 =====

class MatrixAssembler {
public:
    virtual ~MatrixAssembler() = default;
    
    // 添加单元刚度矩阵到全局矩阵
    virtual void addElementMatrix(const std::vector<int>& indices,
                                 const std::vector<std::vector<double>>& element_matrix) = 0;
    
    // 添加单元载荷向量到全局向量
    virtual void addElementVector(const std::vector<int>& indices,
                                 const std::vector<double>& element_vector,
                                 std::vector<double>& global_vector) = 0;
    
    // 应用边界条件
    virtual void applyDirichletBC(int dof_index, double value,
                                 std::shared_ptr<Matrix> matrix,
                                 std::vector<double>& rhs) = 0;
    
    // 应用多个边界条件
    virtual void applyDirichletBCs(const std::vector<int>& dof_indices,
                                  const std::vector<double>& values,
                                  std::shared_ptr<Matrix> matrix,
                                  std::vector<double>& rhs) = 0;
};

// ===== 工厂函数声明 =====

// 创建CRS矩阵
std::shared_ptr<elmer::Matrix> createCRSMatrix(int rows, int cols);

// 创建迭代求解器
std::shared_ptr<IIterativeSolver> createIterativeSolver(std::shared_ptr<Matrix> matrix);

// 创建矩阵组装器
std::shared_ptr<MatrixAssembler> createMatrixAssembler(std::shared_ptr<Matrix> matrix);

// ===== 辅助函数 =====

// 向量点积
double dotProduct(const std::vector<double>& a, const std::vector<double>& b);

// 向量范数
double norm(const std::vector<double>& v);

// 向量加法
std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b);

// 向量标量乘法
std::vector<double> vectorScalarMultiply(const std::vector<double>& v, double scalar);

// 打印向量（调试用）
void printVector(const std::vector<double>& v, const std::string& name = "");

// 打印矩阵（调试用）
void printMatrix(std::shared_ptr<Matrix> matrix, const std::string& name = "");

} // namespace elmer