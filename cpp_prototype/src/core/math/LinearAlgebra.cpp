// LinearAlgebra.cpp - Elmer FEM C++线性代数核心模块
// 对应Fortran模块: CRSMatrix.F90, IterSolve.F90, MatrixAssembly.F90

#include "LinearAlgebra.h"
#include "Types.h"
#include <vector>
#include <memory>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cmath>

using elmer::Integer;
using elmer::Real;

namespace elmer {



// ===== 迭代求解器实现 =====

class IterativeSolverImpl : public IIterativeSolver {
private:
    std::shared_ptr<Matrix> matrix;
    std::vector<double> solution;
    SolverParameters params;
    
public:
    IterativeSolverImpl(std::shared_ptr<Matrix> mat) : matrix(mat) {
        solution.resize(matrix->GetNumRows(), 0.0);
    }
    
    void setParameters(const SolverParameters& parameters) override {
        params = parameters;
    }
    
    SolverResult solve(const std::vector<double>& rhs) override {
        if (rhs.size() != static_cast<size_t>(matrix->GetNumRows())) {
            throw std::invalid_argument("Right-hand side size mismatch");
        }
        
        SolverResult result;
        
        switch (params.method) {
            case SolverMethod::CG:
                result = solveCG(rhs);
                break;
            case SolverMethod::BiCGStab:
                result = solveBiCGStab(rhs);
                break;
            case SolverMethod::GMRES:
                result = solveGMRES(rhs);
                break;
            default:
                throw std::invalid_argument("Unsupported solver method");
        }
        
        return result;
    }
    
    std::vector<double> getSolution() const override {
        return solution;
    }
    
    std::vector<double> getResidualHistory() const override {
        // 简化实现 - 实际求解器中应该记录残差历史
        return {};
    }
    
private:
    // 共轭梯度法 (对应Fortran的CG求解器)
    SolverResult solveCG(const std::vector<double>& rhs) {
        int n = matrix->GetNumRows();
        std::vector<double> x(n, 0.0);  // 初始解
        std::vector<double> r = rhs;    // 残差
        std::vector<double> p = r;      // 搜索方向
        
        double r_norm_sq = dotProduct(r, r);
        double initial_residual = std::sqrt(r_norm_sq);
        
        for (int iter = 0; iter < params.maxIterations; ++iter) {
            auto Ap = Vector::Create(matrix->GetNumRows());
            auto p_vec = Vector::Create(n);
            
            // 将p复制到Vector对象中
            for (int i = 0; i < n; ++i) {
                (*p_vec)[i] = p[i];
            }
            
            matrix->Multiply(*p_vec, *Ap);
            
            // 计算混合类型的点积
            double p_dot_Ap = 0.0;
            for (int i = 0; i < n; ++i) {
                p_dot_Ap += p[i] * (*Ap)[i];
            }
            double alpha = r_norm_sq / p_dot_Ap;
            
            // 更新解和残差
            for (int i = 0; i < n; ++i) {
                x[i] += alpha * p[i];
                r[i] -= alpha * (*Ap)[i];
            }
            
            double new_r_norm_sq = dotProduct(r, r);
            double residual_norm = std::sqrt(new_r_norm_sq);
            
            // 检查收敛
            if (residual_norm < params.tolerance * initial_residual) {
                solution = x;
                return {iter + 1, residual_norm, true};
            }
            
            // 更新搜索方向
            double beta = new_r_norm_sq / r_norm_sq;
            for (int i = 0; i < n; ++i) {
                p[i] = r[i] + beta * p[i];
            }
            
            r_norm_sq = new_r_norm_sq;
        }
        
        solution = x;
        return {params.maxIterations, std::sqrt(r_norm_sq), false};
    }
    
    // BiCGStab求解器 (对应Fortran的BiCGStab求解器)
    SolverResult solveBiCGStab(const std::vector<double>& rhs) {
        int n = matrix->GetNumRows();
        std::vector<double> x(n, 0.0);
        std::vector<double> r = rhs;
        std::vector<double> r0 = r;
        std::vector<double> p = r;
        
        double rho_prev = 1.0, alpha = 1.0, omega = 1.0;
        double initial_residual = norm(r);
        
        for (int iter = 0; iter < params.maxIterations; ++iter) {
            double rho = dotProduct(r0, r);
            
            if (std::abs(rho) < 1e-15) {
                solution = x;
                return {iter, norm(r), false};
            }
            
            if (iter == 0) {
                p = r;
            } else {
                double beta = (rho / rho_prev) * (alpha / omega);
                std::vector<double> temp = matrixMultiply(matrix, p);
                for (int i = 0; i < n; ++i) {
                    p[i] = r[i] + beta * (p[i] - omega * temp[i]);
                }
            }
            
            std::vector<double> v = matrixMultiply(matrix, p);
            alpha = rho / dotProduct(r0, v);
            
            std::vector<double> s(n);
            for (int i = 0; i < n; ++i) {
                s[i] = r[i] - alpha * v[i];
            }
            
            std::vector<double> t = matrixMultiply(matrix, s);
            omega = dotProduct(t, s) / dotProduct(t, t);
            
            // 更新解和残差
            for (int i = 0; i < n; ++i) {
                x[i] += alpha * p[i] + omega * s[i];
                r[i] = s[i] - omega * t[i];
            }
            
            double residual_norm = norm(r);
            if (residual_norm < params.tolerance * initial_residual) {
                solution = x;
                return {iter + 1, residual_norm, true};
            }
            
            rho_prev = rho;
        }
        
        solution = x;
        return {params.maxIterations, norm(r), false};
    }
    
    // GMRES求解器 (简化实现)
    SolverResult solveGMRES(const std::vector<double>& rhs) {
        // 简化实现 - 实际GMRES需要更复杂的Krylov子空间构造
        int n = matrix->GetNumRows();
        std::vector<double> x(n, 0.0);
        std::vector<double> r = rhs;
        
        double initial_residual = norm(r);
        
        for (int iter = 0; iter < params.maxIterations; ++iter) {
            // 简化版 - 实际应使用Arnoldi过程
            std::vector<double> direction = r;
            std::vector<double> Ad = matrixMultiply(matrix, direction);
            
            double alpha = dotProduct(r, Ad) / dotProduct(Ad, Ad);
            
            for (int i = 0; i < n; ++i) {
                x[i] += alpha * direction[i];
                r[i] -= alpha * Ad[i];
            }
            
            double residual_norm = norm(r);
            if (residual_norm < params.tolerance * initial_residual) {
                solution = x;
                return {iter + 1, residual_norm, true};
            }
        }
        
        solution = x;
        return {params.maxIterations, norm(r), false};
    }
    
    // 辅助函数
    double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vector sizes must match for dot product");
        }
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            result += a[i] * b[i];
        }
        return result;
    }
    
    double dotProduct(const Vector& a, const Vector& b) {
        if (a.Size() != b.Size()) {
            throw std::invalid_argument("Vector sizes must match for dot product");
        }
        double result = 0.0;
        for (int i = 0; i < a.Size(); ++i) {
            result += a[i] * b[i];
        }
        return result;
    }
    
    // 矩阵-向量乘法辅助函数
    std::vector<double> matrixMultiply(const std::shared_ptr<elmer::Matrix>& mat, const std::vector<double>& vec) {
        auto result_vec = elmer::Vector::Create(vec.size());
        auto input_vec = elmer::Vector::Create(vec.size());
        
        // 将输入向量复制到Vector对象中
        for (size_t i = 0; i < vec.size(); ++i) {
            (*input_vec)[i] = vec[i];
        }
        
        // 执行矩阵-向量乘法
        mat->Multiply(*input_vec, *result_vec);
        
        // 将结果转换回std::vector<double>
        std::vector<double> result(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            result[i] = (*result_vec)[i];
        }
        
        return result;
    }
    
    double norm(const std::vector<double>& v) {
        return std::sqrt(dotProduct(v, v));
    }
};

// ===== 矩阵组装工具 =====

class MatrixAssemblerImpl : public elmer::MatrixAssembler {
private:
    std::shared_ptr<elmer::Matrix> matrix;
    
public:
    MatrixAssemblerImpl(std::shared_ptr<elmer::Matrix> mat) : matrix(mat) {}
    
    void addElementMatrix(const std::vector<int>& indices,
                         const std::vector<std::vector<double>>& element_matrix) override {
        int n = indices.size();
        
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                int global_i = indices[i];
                int global_j = indices[j];
                
                if (global_i >= 0 && global_j >= 0) {
                    double current_value = matrix->GetElement(global_i, global_j);
                    matrix->SetElement(global_i, global_j, current_value + element_matrix[i][j]);
                }
            }
        }
    }
    
    void addElementVector(const std::vector<int>& indices,
                         const std::vector<double>& element_vector,
                         std::vector<double>& global_vector) override {
        int n = indices.size();
        
        for (int i = 0; i < n; ++i) {
            int global_i = indices[i];
            if (global_i >= 0 && global_i < static_cast<int>(global_vector.size())) {
                global_vector[global_i] += element_vector[i];
            }
        }
    }
    
    void applyDirichletBC(int dof_index, double value,
                          std::shared_ptr<elmer::Matrix> inputMatrix,
                          std::vector<double>& rhs) override {
        if (dof_index < 0 || dof_index >= inputMatrix->GetNumRows()) return;
        
        // 修改矩阵：将对应行清零，对角线设为1
        for (int j = 0; j < inputMatrix->GetNumCols(); ++j) {
            if (j == dof_index) {
                inputMatrix->SetElement(dof_index, j, 1.0);
            } else {
                inputMatrix->SetElement(dof_index, j, 0.0);
            }
        }
        
        // 修改右端项
        rhs[dof_index] = value;
    }
    
    // 应用多个边界条件
    void applyDirichletBCs(const std::vector<int>& dof_indices,
                          const std::vector<double>& values,
                          std::shared_ptr<elmer::Matrix> inputMatrix,
                          std::vector<double>& rhs) override {
        if (dof_indices.size() != values.size()) {
            throw std::invalid_argument("Number of DOF indices and values must match");
        }
        
        for (size_t i = 0; i < dof_indices.size(); ++i) {
            applyDirichletBC(dof_indices[i], values[i], matrix, rhs);
        }
    }
};

// ===== 工厂函数 =====

// ===== Vector实现类 =====

class VectorImpl : public elmer::Vector {
private:
    std::vector<double> data;
    
public:
    VectorImpl(Integer size) : data(size, 0.0) {}
    
    Integer Size() const override { return static_cast<Integer>(data.size()); }
    
    double& operator[](Integer i) override {
        if (i < 0 || i >= static_cast<Integer>(data.size())) {
            throw std::out_of_range("Vector index out of range");
        }
        return data[i];
    }
    
    const double& operator[](Integer i) const override {
        if (i < 0 || i >= static_cast<Integer>(data.size())) {
            throw std::out_of_range("Vector index out of range");
        }
        return data[i];
    }
    
    void Zero() override {
        std::fill(data.begin(), data.end(), 0.0);
    }
};

// Vector工厂方法实现
std::unique_ptr<elmer::Vector> elmer::Vector::Create(Integer size) {
    return std::make_unique<VectorImpl>(size);
}

std::shared_ptr<elmer::Matrix> createCRSMatrix(int rows, int cols) {
    return std::make_shared<elmer::CRSMatrixImpl>(rows, cols);
}

std::shared_ptr<IIterativeSolver> createIterativeSolver(std::shared_ptr<Matrix> matrix) {
    return std::make_shared<IterativeSolverImpl>(matrix);
}

std::shared_ptr<MatrixAssembler> createMatrixAssembler(std::shared_ptr<Matrix> matrix) {
    return std::make_shared<MatrixAssemblerImpl>(matrix);
}

// ===== 辅助函数实现 =====

double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes must match for dot product");
    }
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

double norm(const std::vector<double>& v) {
    return std::sqrt(dotProduct(v, v));
}

std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes must match for addition");
    }
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> vectorScalarMultiply(const std::vector<double>& v, double scalar) {
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i] * scalar;
    }
    return result;
}

void printVector(const std::vector<double>& v, const std::string& name) {
    if (!name.empty()) {
        std::cout << name << " = ";
    }
    std::cout << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << v[i];
        if (i < v.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}

void printMatrix(std::shared_ptr<elmer::Matrix> matrix, const std::string& name) {
    if (!name.empty()) {
        std::cout << name << " = " << std::endl;
    }
    for (int i = 0; i < matrix->GetNumRows(); ++i) {
        std::cout << "[";
        for (int j = 0; j < matrix->GetNumCols(); ++j) {
            std::cout << matrix->GetElement(i, j);
            if (j < matrix->GetNumCols() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
}

} // namespace elmer