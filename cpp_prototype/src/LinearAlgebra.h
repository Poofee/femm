// LinearAlgebra.h - Elmer FEM C++线性代数接口定义
// 对应Fortran模块: CRSMatrix.F90, IterSolve.F90, MatrixAssembly.F90

#pragma once

#include <memory>
#include <vector>
#include <string>

namespace ElmerCpp {

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

// ===== 矩阵接口 =====
class IMatrix {
public:
    virtual ~IMatrix() = default;
    
    // 设置矩阵元素
    virtual void setElement(int i, int j, double value) = 0;
    
    // 设置矩阵元素（别名，用于边界条件）
    virtual void set(int i, int j, double value) = 0;
    
    // 获取矩阵元素
    virtual double getElement(int i, int j) const = 0;
    
    // 获取矩阵元素（别名，用于边界条件）
    virtual double get(int i, int j) const = 0;
    
    // 矩阵清零
    virtual void zero() = 0;
    
    // 清零矩阵行（用于边界条件）
    virtual void zeroRow(int row) = 0;
    
    // 获取矩阵维度
    virtual int getRows() const = 0;
    virtual int getCols() const = 0;
    
    // 矩阵-向量乘法
    virtual std::vector<double> multiply(const std::vector<double>& x) const = 0;
    
    // 获取矩阵非零元素数量
    virtual int getNonzeros() const = 0;
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
class IMatrixAssembler {
public:
    virtual ~IMatrixAssembler() = default;
    
    // 添加单元刚度矩阵到全局矩阵
    virtual void addElementMatrix(const std::vector<int>& indices,
                                 const std::vector<std::vector<double>>& element_matrix) = 0;
    
    // 添加单元载荷向量到全局向量
    virtual void addElementVector(const std::vector<int>& indices,
                                 const std::vector<double>& element_vector,
                                 std::vector<double>& global_vector) = 0;
    
    // 应用边界条件
    virtual void applyDirichletBC(int dof_index, double value,
                                 std::shared_ptr<IMatrix> matrix,
                                 std::vector<double>& rhs) = 0;
    
    // 应用多个边界条件
    virtual void applyDirichletBCs(const std::vector<int>& dof_indices,
                                  const std::vector<double>& values,
                                  std::shared_ptr<IMatrix> matrix,
                                  std::vector<double>& rhs) = 0;
};

// ===== 工厂函数声明 =====

// 创建CRS矩阵
std::shared_ptr<IMatrix> createCRSMatrix(int rows, int cols);

// 创建迭代求解器
std::shared_ptr<IIterativeSolver> createIterativeSolver(std::shared_ptr<IMatrix> matrix);

// 创建矩阵组装器
std::shared_ptr<IMatrixAssembler> createMatrixAssembler(std::shared_ptr<IMatrix> matrix);

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
void printMatrix(std::shared_ptr<IMatrix> matrix, const std::string& name = "");

} // namespace ElmerCpp