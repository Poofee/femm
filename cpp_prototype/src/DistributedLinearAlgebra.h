// DistributedLinearAlgebra.h - 分布式线性代数数据结构定义
// 对应Fortran模块: ParallelMatrix.F90, ParallelVector.F90

#pragma once

#include "MPIConfig.h"
#include "LinearAlgebra.h"
#include <memory>
#include <vector>
#include <map>

namespace elmer {

// ===== 分布式向量类 =====
class DistributedVector : public Vector {
private:
    std::shared_ptr<MPICommunicator> comm_;  // MPI通信器
    int localSize_;                          // 本地向量大小
    int globalSize_;                         // 全局向量大小
    int offset_;                             // 本地向量在全局向量中的偏移
    std::vector<double> localData_;          // 本地数据
    
    // 重叠区域管理
    std::map<int, std::vector<double>> ghostData_;  // 幽灵数据（来自相邻进程）
    std::vector<int> ghostIndices_;                 // 幽灵索引
    
public:
    /**
     * @brief 构造函数
     * @param globalSize 全局向量大小
     * @param comm MPI通信器
     */
    DistributedVector(int globalSize, std::shared_ptr<MPICommunicator> comm = nullptr);
    
    /**
     * @brief 构造函数（指定本地大小和偏移）
     */
    DistributedVector(int localSize, int globalSize, int offset, 
                     std::shared_ptr<MPICommunicator> comm = nullptr);
    
    // 向量接口实现
    void SetElement(Integer i, Real value) override;
    Real GetElement(Integer i) const override;
    void AddToElement(Integer i, Real value) override;
    Integer Size() const override;
    
    // 分布式特定方法
    int getLocalSize() const { return localSize_; }
    int getGlobalSize() const { return globalSize_; }
    int getOffset() const { return offset_; }
    
    // 重叠区域管理
    void addGhostIndex(int globalIndex);
    void updateGhostData();
    double getGhostValue(int globalIndex) const;
    
    // 通信操作
    void gatherToMaster(std::vector<double>& globalVector) const;
    void scatterFromMaster(const std::vector<double>& globalVector);
    void allGather(std::vector<double>& globalVector) const;
    
    // 向量操作
    void setZero();
    void scale(double factor);
    double norm() const;
    double dot(const DistributedVector& other) const;
    
    // 获取本地数据指针
    const double* getLocalData() const { return localData_.data(); }
    double* getLocalData() { return localData_.data(); }
    
    // 获取通信器
    std::shared_ptr<MPICommunicator> getCommunicator() const { return comm_; }
};

// ===== 分布式矩阵类 =====
class DistributedMatrix : public Matrix {
private:
    std::shared_ptr<MPICommunicator> comm_;  // MPI通信器
    int localRows_;                          // 本地行数
    int globalRows_;                         // 全局行数
    int localCols_;                          // 本地列数
    int globalCols_;                         // 全局列数
    int rowOffset_;                          // 行偏移
    int colOffset_;                          // 列偏移
    
    // 本地矩阵存储
    std::shared_ptr<Matrix> localMatrix_;    // 本地矩阵
    
    // 重叠区域管理
    std::vector<int> ghostRows_;             // 幽灵行索引
    std::vector<int> ghostCols_;             // 幽灵列索引
    
public:
    /**
     * @brief 构造函数
     * @param globalRows 全局行数
     * @param globalCols 全局列数
     * @param comm MPI通信器
     */
    DistributedMatrix(int globalRows, int globalCols, 
                     std::shared_ptr<MPICommunicator> comm = nullptr);
    
    /**
     * @brief 构造函数（指定本地大小和偏移）
     */
    DistributedMatrix(int localRows, int localCols, int globalRows, int globalCols,
                     int rowOffset, int colOffset, 
                     std::shared_ptr<MPICommunicator> comm = nullptr);
    
    // 矩阵接口实现
    void SetElement(Integer i, Integer j, Real value) override;
    Real GetElement(Integer i, Integer j) const override;
    void AddToElement(Integer i, Integer j, Real value) override;
    Integer GetRows() const override;
    Integer GetCols() const override;
    
    // 分布式特定方法
    int getLocalRows() const { return localRows_; }
    int getGlobalRows() const { return globalRows_; }
    int getLocalCols() const { return localCols_; }
    int getGlobalCols() const { return globalCols_; }
    int getRowOffset() const { return rowOffset_; }
    int getColOffset() const { return colOffset_; }
    
    // 矩阵-向量乘法
    void multiply(const DistributedVector& x, DistributedVector& y) const;
    
    // 通信操作
    void assemble();  // 组装重叠区域数据
    
    // 获取本地矩阵
    std::shared_ptr<Matrix> getLocalMatrix() const { return localMatrix_; }
    
    // 获取通信器
    std::shared_ptr<MPICommunicator> getCommunicator() const { return comm_; }
};

// ===== 分布式线性系统类 =====
class DistributedLinearSystem {
private:
    std::shared_ptr<DistributedMatrix> A_;   // 分布式矩阵
    std::shared_ptr<DistributedVector> b_;   // 分布式右端向量
    std::shared_ptr<DistributedVector> x_;   // 分布式解向量
    std::shared_ptr<MPICommunicator> comm_;  // MPI通信器
    
public:
    /**
     * @brief 构造函数
     * @param globalSize 全局系统大小
     * @param comm MPI通信器
     */
    DistributedLinearSystem(int globalSize, std::shared_ptr<MPICommunicator> comm = nullptr);
    
    /**
     * @brief 构造函数（指定本地大小和偏移）
     */
    DistributedLinearSystem(int localSize, int globalSize, int offset,
                           std::shared_ptr<MPICommunicator> comm = nullptr);
    
    // 获取系统组件
    std::shared_ptr<DistributedMatrix> getMatrix() const { return A_; }
    std::shared_ptr<DistributedVector> getRHS() const { return b_; }
    std::shared_ptr<DistributedVector> getSolution() const { return x_; }
    
    // 系统操作
    void setZero();
    void assemble();
    
    // 获取通信器
    std::shared_ptr<MPICommunicator> getCommunicator() const { return comm_; }
};

// ===== 分布式线性代数工具函数 =====
namespace DistributedLinearAlgebraUtils {
    
    // 创建均匀分布的向量
    std::shared_ptr<DistributedVector> createUniformDistributedVector(int globalSize, 
                                                                     std::shared_ptr<MPICommunicator> comm = nullptr);
    
    // 创建均匀分布的矩阵
    std::shared_ptr<DistributedMatrix> createUniformDistributedMatrix(int globalRows, int globalCols,
                                                                     std::shared_ptr<MPICommunicator> comm = nullptr);
    
    // 创建均匀分布的线性系统
    std::shared_ptr<DistributedLinearSystem> createUniformDistributedLinearSystem(int globalSize,
                                                                                 std::shared_ptr<MPICommunicator> comm = nullptr);
    
    // 计算全局向量范数
    double globalNorm(const DistributedVector& vec);
    
    // 计算全局向量点积
    double globalDot(const DistributedVector& vec1, const DistributedVector& vec2);
    
    // 全局向量复制
    void globalCopy(const DistributedVector& src, DistributedVector& dest);
    
    // 全局向量缩放
    void globalScale(DistributedVector& vec, double factor);
    
    // 全局向量加
    void globalAdd(DistributedVector& vec1, const DistributedVector& vec2);
    
    // 全局向量加（带缩放）
    void globalAddScaled(DistributedVector& vec1, double alpha, const DistributedVector& vec2);
    
    // 检查全局收敛性
    bool checkGlobalConvergence(const DistributedVector& residual, double tolerance);
    
} // namespace DistributedLinearAlgebraUtils

} // namespace elmer