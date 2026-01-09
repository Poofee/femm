// MatrixAssembly.h - Elmer FEM C++矩阵组装模块
// 对应Fortran模块: MatrixAssembly.F90

#pragma once

#include "ElmerCpp.h"
#include <memory>
#include <vector>
#include <stdexcept>

namespace elmer {

// 矩阵格式枚举
enum class MatrixFormat {
    CRS,        // 压缩行存储格式
    LIST,       // 列表存储格式
    BAND,       // 带状矩阵格式
    SBAND       // 对称带状矩阵格式
};

// 矩阵组装类
class MatrixAssembly {
public:
    // 构造函数
    MatrixAssembly() = default;
    
    // 设置矩阵元素
    static void SetMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                elmer::Integer i, elmer::Integer j, elmer::Real value);
    
    // 获取矩阵元素
    static elmer::Real GetMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                elmer::Integer i, elmer::Integer j);
    
    // 修改矩阵元素并返回旧值
    static elmer::Real ChangeMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                   elmer::Integer i, elmer::Integer j, elmer::Real newValue);
    
    // 添加到矩阵元素
    static void AddToMatrixElement(std::shared_ptr<elmer::Matrix> matrix, 
                                  elmer::Integer i, elmer::Integer j, elmer::Real value);
    
    // 复数矩阵元素操作
    static void AddToComplexMatrixElement(std::shared_ptr<Matrix> matrix,
                                         elmer::Integer rowId, elmer::Integer colId,
                                         elmer::Real re, elmer::Real im);
    
    // 移动矩阵元素
    static void MoveMatrixElement(std::shared_ptr<Matrix> matrix,
                                 elmer::Integer i1, elmer::Integer j1, elmer::Integer i2, elmer::Integer j2);
    
    // 矩阵清零
    static void ZeroMatrix(std::shared_ptr<Matrix> matrix);
    
    // 矩阵缩放
    static void ScaleMatrix(std::shared_ptr<Matrix> matrix, elmer::Real factor);
    
    // 矩阵-向量乘法
    static void MatrixVectorMultiply(std::shared_ptr<elmer::Matrix> matrix,
                                    std::shared_ptr<elmer::Vector> vector,
                                    std::shared_ptr<elmer::Vector> result);
    
    // 有限元矩阵组装相关函数
    
    // 组装单元刚度矩阵到全局矩阵
    static void AssembleElementMatrix(std::shared_ptr<Matrix> globalMatrix,
                                     const std::vector<elmer::Integer>& dofIndices,
                                     const std::vector<std::vector<elmer::Real>>& elementMatrix);
    
    // 组装单元载荷向量到全局向量
    static void AssembleElementVector(std::shared_ptr<elmer::Vector> globalVector,
                                     const std::vector<elmer::Integer>& dofIndices,
                                     const std::vector<elmer::Real>& elementVector);
    
    // 应用Dirichlet边界条件
    static void ApplyDirichletBC(std::shared_ptr<elmer::Matrix> matrix,
                                std::vector<elmer::Real>& rhs,
                                elmer::Integer dofIndex,
                                elmer::Real prescribedValue);
    
    // 应用Dirichlet边界条件（多个自由度）
    static void ApplyDirichletBC(std::shared_ptr<elmer::Matrix> matrix,
                                std::vector<elmer::Real>& rhs,
                                const std::vector<elmer::Integer>& dofIndices,
                                const std::vector<elmer::Real>& prescribedValues);
    
    // 矩阵格式转换
    static std::shared_ptr<elmer::Matrix> ConvertMatrixFormat(std::shared_ptr<elmer::Matrix> source,
                                                      elmer::MatrixFormat targetFormat);
    
    // 矩阵信息查询
    static elmer::Integer GetMatrixSize(std::shared_ptr<elmer::Matrix> matrix);
    static elmer::Integer GetMatrixNonzeros(std::shared_ptr<elmer::Matrix> matrix);
    static elmer::MatrixFormat GetMatrixFormat(std::shared_ptr<elmer::Matrix> matrix);
    
    // 获取矩阵行
    static std::vector<elmer::Real> GetMatrixRow(std::shared_ptr<elmer::Matrix> matrix, elmer::Integer row);
    
    // 设置矩阵行
    static void SetMatrixRow(std::shared_ptr<elmer::Matrix> matrix, elmer::Integer row, 
                            const std::vector<elmer::Real>& values);
    
    // 矩阵条件数估计
    static elmer::Real EstimateConditionNumber(std::shared_ptr<elmer::Matrix> matrix);
    
    // 矩阵对称性检查
    static bool IsSymmetric(std::shared_ptr<elmer::Matrix> matrix, elmer::Real tolerance = 1e-12);
    
    // 矩阵正定性检查
    static bool IsPositiveDefinite(std::shared_ptr<elmer::Matrix> matrix);
    
private:
    // 内部辅助函数
    
    // 检查索引有效性
    static void CheckIndices(std::shared_ptr<elmer::Matrix> matrix, elmer::Integer i, elmer::Integer j);
    
    // CRS矩阵特定操作
    static void CRS_SetMatrixElement(std::shared_ptr<Matrix> matrix,
                                    elmer::Integer i, elmer::Integer j, elmer::Real value);
    static elmer::Real CRS_GetMatrixElement(std::shared_ptr<Matrix> matrix,
                                    elmer::Integer i, elmer::Integer j);
    static void CRS_AddToMatrixElement(std::shared_ptr<Matrix> matrix,
                                      elmer::Integer i, elmer::Integer j, elmer::Real value);
    
    // 带状矩阵特定操作
    static void Band_SetMatrixElement(std::shared_ptr<Matrix> matrix,
                                     elmer::Integer i, elmer::Integer j, elmer::Real value);
    static elmer::Real Band_GetMatrixElement(std::shared_ptr<Matrix> matrix,
                                     elmer::Integer i, elmer::Integer j);
    static void Band_AddToMatrixElement(std::shared_ptr<Matrix> matrix,
                                       elmer::Integer i, elmer::Integer j, elmer::Real value);
    
    // 列表矩阵特定操作
    static void List_SetMatrixElement(std::shared_ptr<Matrix> matrix,
                                     elmer::Integer i, elmer::Integer j, elmer::Real value);
    static elmer::Real List_GetMatrixElement(std::shared_ptr<Matrix> matrix,
                                     elmer::Integer i, elmer::Integer j);
    static void List_AddToMatrixElement(std::shared_ptr<Matrix> matrix,
                                       elmer::Integer i, elmer::Integer j, elmer::Real value);
};

} // namespace elmer