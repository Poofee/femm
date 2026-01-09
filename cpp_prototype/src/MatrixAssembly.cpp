// MatrixAssembly.cpp - Elmer FEM C++矩阵组装模块实现
// 对应Fortran模块: MatrixAssembly.F90

#include "MatrixAssembly.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace elmer {

// 设置矩阵元素
void MatrixAssembly::SetMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                     elmer::Integer i, elmer::Integer j, elmer::Real value) {
    CheckIndices(matrix, i, j);
    
    // 根据矩阵格式调用相应的实现
    // 这里需要根据实际的矩阵类型进行分派
    // 目前假设所有矩阵都支持SetElement方法
    matrix->SetElement(i, j, value);
}

// 获取矩阵元素
elmer::Real MatrixAssembly::GetMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                     elmer::Integer i, elmer::Integer j) {
    CheckIndices(matrix, i, j);
    
    // 根据矩阵格式调用相应的实现
    return matrix->GetElement(i, j);
}

// 修改矩阵元素并返回旧值
elmer::Real MatrixAssembly::ChangeMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                        elmer::Integer i, elmer::Integer j, elmer::Real newValue) {
    CheckIndices(matrix, i, j);
    
    elmer::Real oldValue = GetMatrixElement(matrix, i, j);
    SetMatrixElement(matrix, i, j, newValue);
    return oldValue;
}

// 添加到矩阵元素
void MatrixAssembly::AddToMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                       elmer::Integer i, elmer::Integer j, elmer::Real value) {
    CheckIndices(matrix, i, j);
    
    // 根据矩阵格式调用相应的实现
    // 这里需要根据实际的矩阵类型进行分派
    // 目前假设所有矩阵都支持AddToElement方法
    matrix->AddToElement(i, j, value);
}

// 复数矩阵元素操作
void MatrixAssembly::AddToComplexMatrixElement(std::shared_ptr<elmer::Matrix> matrix,
                                              elmer::Integer rowId, elmer::Integer colId,
                                              elmer::Real re, elmer::Real im) {
    // 复数矩阵存储为2x2块
    AddToMatrixElement(matrix, rowId, colId, re);
    AddToMatrixElement(matrix, rowId, colId + 1, -im);
    AddToMatrixElement(matrix, rowId + 1, colId, im);
    AddToMatrixElement(matrix, rowId + 1, colId + 1, re);
}

// 移动矩阵元素
void MatrixAssembly::MoveMatrixElement(std::shared_ptr<elmer::Matrix> matrix,
                                      elmer::Integer i1, elmer::Integer j1, elmer::Integer i2, elmer::Integer j2) {
    CheckIndices(matrix, i1, j1);
    CheckIndices(matrix, i2, j2);
    
    elmer::Real value = ChangeMatrixElement(matrix, i1, j1, 0.0);
    AddToMatrixElement(matrix, i2, j2, value);
}

// 矩阵清零
void MatrixAssembly::ZeroMatrix(std::shared_ptr<elmer::Matrix> matrix) {
    matrix->Zero();
}

// 矩阵缩放
void MatrixAssembly::ScaleMatrix(std::shared_ptr<elmer::Matrix> matrix, elmer::Real factor) {
    // 这里需要根据矩阵格式实现缩放操作
    // 目前假设矩阵支持遍历所有元素
    elmer::Integer size = GetMatrixSize(matrix);
    
    for (elmer::Integer i = 0; i < size; ++i) {
        for (elmer::Integer j = 0; j < size; ++j) {
            elmer::Real value = GetMatrixElement(matrix, i, j);
            if (value != 0.0) {
                SetMatrixElement(matrix, i, j, value * factor);
            }
        }
    }
}

// 矩阵向量乘法
void MatrixAssembly::MatrixVectorMultiply(std::shared_ptr<elmer::Matrix> matrix,
                                         const std::vector<elmer::Real>& x,
                                         std::vector<elmer::Real>& y) {
    elmer::Integer size = GetMatrixSize(matrix);
    
    if (x.size() != static_cast<size_t>(size) || y.size() != static_cast<size_t>(size)) {
        throw std::invalid_argument("Matrix and vector sizes do not match");
    }
    
    // 清零结果向量
    std::fill(y.begin(), y.end(), 0.0);
    
    // 执行矩阵向量乘法
    for (elmer::Integer i = 0; i < size; ++i) {
        for (elmer::Integer j = 0; j < size; ++j) {
            elmer::Real value = GetMatrixElement(matrix, i, j);
            if (value != 0.0) {
                y[i] += value * x[j];
            }
        }
    }
}

// 组装单元刚度矩阵到全局矩阵
void MatrixAssembly::AssembleElementMatrix(std::shared_ptr<elmer::Matrix> globalMatrix,
                                          const std::vector<elmer::Integer>& dofIndices,
                                          const std::vector<std::vector<elmer::Real>>& elementMatrix) {
    elmer::Integer numDofs = dofIndices.size();
    
    if (elementMatrix.size() != static_cast<size_t>(numDofs) || 
        (numDofs > 0 && elementMatrix[0].size() != static_cast<size_t>(numDofs))) {
        throw std::invalid_argument("Element matrix dimensions do not match DOF indices");
    }
    
    // 组装单元矩阵到全局矩阵
    for (elmer::Integer i = 0; i < numDofs; ++i) {
        elmer::Integer globalI = dofIndices[i];
        for (elmer::Integer j = 0; j < numDofs; ++j) {
            elmer::Integer globalJ = dofIndices[j];
            elmer::Real value = elementMatrix[i][j];
            
            if (value != 0.0) {
                AddToMatrixElement(globalMatrix, globalI, globalJ, value);
            }
        }
    }
}

// 组装单元载荷向量到全局向量
void MatrixAssembly::AssembleElementVector(std::vector<elmer::Real>& globalVector,
                                          const std::vector<elmer::Integer>& dofIndices,
                                          const std::vector<elmer::Real>& elementVector) {
    elmer::Integer numDofs = dofIndices.size();
    
    if (elementVector.size() != static_cast<size_t>(numDofs)) {
        throw std::invalid_argument("Element vector size does not match DOF indices");
    }
    
    // 组装单元向量到全局向量
    for (elmer::Integer i = 0; i < numDofs; ++i) {
        elmer::Integer globalI = dofIndices[i];
        globalVector[globalI] += elementVector[i];
    }
}

// 应用Dirichlet边界条件
void MatrixAssembly::ApplyDirichletBC(std::shared_ptr<elmer::Matrix> matrix,
                                     std::vector<elmer::Real>& rhs,
                                     elmer::Integer dofIndex, elmer::Real prescribedValue) {
    CheckIndices(matrix, dofIndex, dofIndex);
    
    elmer::Integer size = GetMatrixSize(matrix);
    
    // 修改矩阵：将对应行和列清零，对角线设为1
    for (elmer::Integer j = 0; j < size; ++j) {
        if (j != dofIndex) {
            SetMatrixElement(matrix, dofIndex, j, 0.0);
            SetMatrixElement(matrix, j, dofIndex, 0.0);
        }
    }
    SetMatrixElement(matrix, dofIndex, dofIndex, 1.0);
    
    // 修改右端项
    rhs[dofIndex] = prescribedValue;
}

// 应用Dirichlet边界条件（多个自由度）
void MatrixAssembly::ApplyDirichletBC(std::shared_ptr<elmer::Matrix> matrix,
                                     std::vector<elmer::Real>& rhs,
                                     const std::vector<elmer::Integer>& dofIndices,
                                     const std::vector<elmer::Real>& prescribedValues) {
    if (dofIndices.size() != prescribedValues.size()) {
        throw std::invalid_argument("Number of DOF indices and prescribed values must match");
    }
    
    for (size_t idx = 0; idx < dofIndices.size(); ++idx) {
        ApplyDirichletBC(matrix, rhs, dofIndices[idx], prescribedValues[idx]);
    }
}

// 矩阵格式转换（简化实现）
std::shared_ptr<elmer::Matrix> MatrixAssembly::ConvertMatrixFormat(std::shared_ptr<elmer::Matrix> source,
                                                       elmer::MatrixFormat targetFormat) {
    // 这里需要根据源矩阵和目标格式实现转换
    // 目前返回源矩阵（简化实现）
    return source;
}

// 矩阵信息查询
elmer::Integer MatrixAssembly::GetMatrixSize(std::shared_ptr<elmer::Matrix> matrix) {
    // 这里需要根据矩阵类型返回大小
    // 目前假设矩阵是方阵，通过查询对角线元素数量来确定大小
    elmer::Integer size = 0;
    while (true) {
        try {
            matrix->GetElement(size, size);
            ++size;
        } catch (...) {
            break;
        }
    }
    return size;
}

elmer::Integer MatrixAssembly::GetMatrixNonzeros(std::shared_ptr<elmer::Matrix> matrix) {
    elmer::Integer size = GetMatrixSize(matrix);
    elmer::Integer nonzeros = 0;
    
    for (elmer::Integer i = 0; i < size; ++i) {
        for (elmer::Integer j = 0; j < size; ++j) {
            if (GetMatrixElement(matrix, i, j) != 0.0) {
                ++nonzeros;
            }
        }
    }
    
    return nonzeros;
}
// 矩阵格式获取
elmer::MatrixFormat MatrixAssembly::GetMatrixFormat(std::shared_ptr<elmer::Matrix> matrix) {
    // 这里需要根据矩阵类型确定格式
    // 目前返回CRS格式（简化实现）
    return elmer::MatrixFormat::CRS;
}

// 矩阵条件数估计（简化实现）
elmer::Real MatrixAssembly::EstimateConditionNumber(std::shared_ptr<elmer::Matrix> matrix) {
    // 这里需要实现更精确的条件数估计
    // 目前返回一个估计值
    elmer::Integer size = GetMatrixSize(matrix);
    
    // 计算Frobenius范数作为条件数的粗略估计
    elmer::Real norm = 0.0;
    for (elmer::Integer i = 0; i < size; ++i) {
        for (elmer::Integer j = 0; j < size; ++j) {
            elmer::Real value = GetMatrixElement(matrix, i, j);
            norm += value * value;
        }
    }
    norm = std::sqrt(norm);
    
    // 简化条件数估计
    return norm * norm; // 近似条件数
}

// 矩阵对称性检查
bool MatrixAssembly::IsSymmetric(std::shared_ptr<elmer::Matrix> matrix, elmer::Real tolerance) {
    elmer::Integer size = GetMatrixSize(matrix);
    
    // 检查矩阵是否对称
    for (elmer::Integer i = 0; i < size; ++i) {
        for (elmer::Integer j = i + 1; j < size; ++j) {
            elmer::Real a_ij = GetMatrixElement(matrix, i, j);
            elmer::Real a_ji = GetMatrixElement(matrix, j, i);
            
            if (std::abs(a_ij - a_ji) > tolerance) {
                return false;
            }
        }
    }
    
    return true;
}

// 矩阵正定性检查（简化实现）
bool MatrixAssembly::IsPositiveDefinite(std::shared_ptr<elmer::Matrix> matrix) {
    // 这里需要实现更精确的正定性检查
    // 目前检查对角线元素是否为正
    elmer::Integer size = GetMatrixSize(matrix);
    
    for (elmer::Integer i = 0; i < size; ++i) {
        elmer::Real diag = GetMatrixElement(matrix, i, i);
        if (diag <= 0.0) {
            return false;
        }
    }
    
    return true;
}

// 索引检查
void MatrixAssembly::CheckIndices(std::shared_ptr<elmer::Matrix> matrix, elmer::Integer i, elmer::Integer j) {
    elmer::Integer size = GetMatrixSize(matrix);
    
    if (i < 0 || i >= size || j < 0 || j >= size) {
        throw std::out_of_range("Matrix indices out of range");
    }
}

// 矩阵操作列表函数（简化实现）
elmer::Real MatrixAssembly::List_GetMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                             elmer::Integer i, elmer::Integer j) {
    return GetMatrixElement(matrix, i, j);
}

void MatrixAssembly::List_SetMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                            elmer::Integer i, elmer::Integer j, elmer::Real value) {
    SetMatrixElement(matrix, i, j, value);
}

void MatrixAssembly::List_AddToMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                             elmer::Integer i, elmer::Integer j, elmer::Real value) {
    AddToMatrixElement(matrix, i, j, value);
}

elmer::Real MatrixAssembly::List_ChangeMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                                  elmer::Integer i, elmer::Integer j, elmer::Real newValue) {
    return ChangeMatrixElement(matrix, i, j, newValue);
}

void MatrixAssembly::List_MoveMatrixElement(std::shared_ptr<Matrix> matrix, 
                                           elmer::Integer i1, elmer::Integer j1, 
                                           elmer::Integer i2, elmer::Integer j2) {
    MoveMatrixElement(matrix, i1, j1, i2, j2);
}

void MatrixAssembly::List_ZeroMatrix(std::shared_ptr<Matrix> matrix) {
    ZeroMatrix(matrix);
}

void MatrixAssembly::List_ScaleMatrix(std::shared_ptr<Matrix> matrix, elmer::Real factor) {
    ScaleMatrix(matrix, factor);
}

void MatrixAssembly::List_MatrixVectorMultiply(std::shared_ptr<Matrix> matrix, 
                                              const std::vector<elmer::Real>& x, 
                                              std::vector<elmer::Real>& y) {
    MatrixVectorMultiply(matrix, x, y);
}

void MatrixAssembly::List_AssembleElementMatrix(std::shared_ptr<Matrix> globalMatrix, 
                                               const std::vector<elmer::Integer>& dofIndices, 
                                               const std::vector<std::vector<elmer::Real>>& elementMatrix) {
    AssembleElementMatrix(globalMatrix, dofIndices, elementMatrix);
}

void MatrixAssembly::List_AssembleElementVector(std::vector<elmer::Real>& globalVector, 
                                               const std::vector<elmer::Integer>& dofIndices, 
                                               const std::vector<elmer::Real>& elementVector) {
    AssembleElementVector(globalVector, dofIndices, elementVector);
}

void MatrixAssembly::List_ApplyDirichletBC(std::shared_ptr<Matrix> matrix, 
                                          std::vector<elmer::Real>& rhs, 
                                          elmer::Integer dofIndex, elmer::Real prescribedValue) {
    ApplyDirichletBC(matrix, rhs, dofIndex, prescribedValue);
}

void MatrixAssembly::List_ApplyDirichletBC_Multi(std::shared_ptr<Matrix> matrix, 
                                                std::vector<elmer::Real>& rhs, 
                                                const std::vector<elmer::Integer>& dofIndices, 
                                                const std::vector<elmer::Real>& prescribedValues) {
    ApplyDirichletBC(matrix, rhs, dofIndices, prescribedValues);
}

elmer::Integer MatrixAssembly::List_GetMatrixSize(std::shared_ptr<Matrix> matrix) {
    return GetMatrixSize(matrix);
}

elmer::Integer MatrixAssembly::List_GetMatrixNonzeros(std::shared_ptr<elmer::Matrix> matrix) {
    return GetMatrixNonzeros(matrix);
}

elmer::MatrixFormat MatrixAssembly::List_GetMatrixFormat(std::shared_ptr<elmer::Matrix> matrix) {
    return GetMatrixFormat(matrix);
}

elmer::Real MatrixAssembly::List_EstimateConditionNumber(std::shared_ptr<elmer::Matrix> matrix) {
    return EstimateConditionNumber(matrix);
}

bool MatrixAssembly::List_IsSymmetric(std::shared_ptr<Matrix> matrix, elmer::Real tolerance) {
    return IsSymmetric(matrix, tolerance);
}

bool MatrixAssembly::List_IsPositiveDefinite(std::shared_ptr<Matrix> matrix) {
    return IsPositiveDefinite(matrix);
}

} // namespace elmer

// 检查索引有效性
void MatrixAssembly::CheckIndices(std::shared_ptr<Matrix> matrix, Integer i, Integer j) {
    Integer size = GetMatrixSize(matrix);
    
    if (i < 0 || i >= size || j < 0 || j >= size) {
        throw std::out_of_range("Matrix index out of range");
    }
}

// CRS矩阵特定操作（需要CRS矩阵的具体实现）
void MatrixAssembly::CRS_SetMatrixElement(std::shared_ptr<Matrix> matrix,
                                         Integer i, Integer j, Real value) {
    // 这里需要CRS矩阵的具体实现
    matrix->SetElement(i, j, value);
}

Real MatrixAssembly::CRS_GetMatrixElement(std::shared_ptr<Matrix> matrix,
                                         Integer i, Integer j) {
    // 这里需要CRS矩阵的具体实现
    return matrix->GetElement(i, j);
}

void MatrixAssembly::CRS_AddToMatrixElement(std::shared_ptr<Matrix> matrix,
                                           Integer i, Integer j, Real value) {
    // 这里需要CRS矩阵的具体实现
    matrix->AddToElement(i, j, value);
}

// 带状矩阵特定操作（需要带状矩阵的具体实现）
void MatrixAssembly::Band_SetMatrixElement(std::shared_ptr<Matrix> matrix,
                                          Integer i, Integer j, Real value) {
    // 这里需要带状矩阵的具体实现
    matrix->SetElement(i, j, value);
}

Real MatrixAssembly::Band_GetMatrixElement(std::shared_ptr<Matrix> matrix,
                                          Integer i, Integer j) {
    // 这里需要带状矩阵的具体实现
    return matrix->GetElement(i, j);
}

void MatrixAssembly::Band_AddToMatrixElement(std::shared_ptr<Matrix> matrix,
                                            Integer i, Integer j, Real value) {
    // 这里需要带状矩阵的具体实现
    matrix->AddToElement(i, j, value);
}

// 列表矩阵特定操作（需要列表矩阵的具体实现）
void MatrixAssembly::List_SetMatrixElement(std::shared_ptr<Matrix> matrix,
                                          Integer i, Integer j, Real value) {
    // 这里需要列表矩阵的具体实现
    matrix->SetElement(i, j, value);
}

Real MatrixAssembly::List_GetMatrixElement(std::shared_ptr<Matrix> matrix,
                                          Integer i, Integer j) {
    // 这里需要列表矩阵的具体实现
    return matrix->GetElement(i, j);
}

void MatrixAssembly::List_AddToMatrixElement(std::shared_ptr<Matrix> matrix,
                                            Integer i, Integer j, Real value) {
    // 这里需要列表矩阵的具体实现
    matrix->AddToElement(i, j, value);
}

} // namespace elmer