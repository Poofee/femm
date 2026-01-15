// DistributedLinearAlgebra.cpp - 分布式线性代数数据结构实现
// 对应Fortran模块: ParallelMatrix.F90, ParallelVector.F90

#include "DistributedLinearAlgebra.h"
#include "CRSMatrix.h"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace elmer {

// ===== DistributedVector实现 =====
DistributedVector::DistributedVector(int globalSize, std::shared_ptr<MPICommunicator> comm)
    : comm_(comm ? comm : MPIUtils::getDefaultComm()), globalSize_(globalSize) {
    
    if (!comm_) {
        // 如果没有MPI通信器，使用串行模式
        localSize_ = globalSize;
        offset_ = 0;
    } else {
        // 计算本地大小和偏移
        int rank = comm_->getRank();
        int size = comm_->getSize();
        
        localSize_ = globalSize / size;
        int remainder = globalSize % size;
        
        if (rank < remainder) {
            localSize_++;
            offset_ = rank * localSize_;
        } else {
            offset_ = rank * localSize_ + remainder;
        }
    }
    
    localData_.resize(localSize_, 0.0);
}

DistributedVector::DistributedVector(int localSize, int globalSize, int offset, 
                                   std::shared_ptr<MPICommunicator> comm)
    : comm_(comm ? comm : MPIUtils::getDefaultComm()), 
      localSize_(localSize), globalSize_(globalSize), offset_(offset) {
    
    localData_.resize(localSize_, 0.0);
}

void DistributedVector::SetElement(Integer i, Real value) {
    if (i < offset_ || i >= offset_ + localSize_) {
        throw std::out_of_range("DistributedVector::SetElement: index out of local range");
    }
    localData_[i - offset_] = value;
}

Real DistributedVector::GetElement(Integer i) const {
    if (i >= offset_ && i < offset_ + localSize_) {
        return localData_[i - offset_];
    } else {
        // 检查是否是幽灵数据
        auto it = ghostData_.find(i);
        if (it != ghostData_.end()) {
            return it->second[0]; // 返回幽灵数据
        }
        throw std::out_of_range("DistributedVector::GetElement: index not found");
    }
}

void DistributedVector::AddToElement(Integer i, Real value) {
    if (i < offset_ || i >= offset_ + localSize_) {
        throw std::out_of_range("DistributedVector::AddToElement: index out of local range");
    }
    localData_[i - offset_] += value;
}

Integer DistributedVector::Size() const {
    return globalSize_;
}

void DistributedVector::addGhostIndex(int globalIndex) {
    ghostIndices_.push_back(globalIndex);
}

void DistributedVector::updateGhostData() {
    if (!comm_ || ghostIndices_.empty()) return;
    
    // 简化实现：通过主进程收集所有数据，然后广播
    // 在实际应用中，应该使用更高效的点对点通信
    
    if (comm_->isMaster()) {
        // 主进程收集所有数据
        std::vector<double> globalData(globalSize_);
        gatherToMaster(globalData);
        
        // 广播幽灵数据
        for (int index : ghostIndices_) {
            if (index >= 0 && index < globalSize_) {
                double value = globalData[index];
                comm_->broadcast(&value, 1, MPI_DOUBLE);
            }
        }
    } else {
        // 从进程收集数据到主进程
        gatherToMaster(std::vector<double>());
        
        // 接收幽灵数据
        for (int index : ghostIndices_) {
            if (index >= 0 && index < globalSize_) {
                double value;
                comm_->broadcast(&value, 1, MPI_DOUBLE);
                ghostData_[index] = {value};
            }
        }
    }
}

double DistributedVector::getGhostValue(int globalIndex) const {
    auto it = ghostData_.find(globalIndex);
    if (it != ghostData_.end()) {
        return it->second[0];
    }
    throw std::runtime_error("Ghost value not found for index: " + std::to_string(globalIndex));
}

void DistributedVector::gatherToMaster(std::vector<double>& globalVector) const {
    if (!comm_) {
        globalVector = localData_;
        return;
    }
    
    if (comm_->isMaster()) {
        globalVector.resize(globalSize_);
        
        // 收集所有进程的数据
        std::vector<int> recvCounts(comm_->getSize());
        std::vector<int> displs(comm_->getSize());
        
        // 首先收集每个进程的本地大小
        std::vector<int> localSizes(comm_->getSize());
        comm_->gather(&localSize_, 1, MPI_INT, localSizes.data(), 0);
        
        // 计算位移
        displs[0] = 0;
        for (int i = 1; i < comm_->getSize(); ++i) {
            displs[i] = displs[i-1] + localSizes[i-1];
        }
        
        // 收集数据
        comm_->gather(localData_.data(), localSize_, MPI_DOUBLE, 
                     globalVector.data(), localSizes.data(), displs.data(), 0);
    } else {
        // 从进程发送数据
        comm_->gather(localData_.data(), localSize_, MPI_DOUBLE, nullptr, 0);
    }
}

void DistributedVector::scatterFromMaster(const std::vector<double>& globalVector) {
    if (!comm_) {
        if (globalVector.size() == localData_.size()) {
            localData_ = globalVector;
        }
        return;
    }
    
    if (comm_->isMaster()) {
        if (globalVector.size() != globalSize_) {
            throw std::runtime_error("Global vector size mismatch");
        }
        
        // 分散数据到各进程
        std::vector<int> sendCounts(comm_->getSize());
        std::vector<int> displs(comm_->getSize());
        
        // 计算发送计数和位移
        int remainder = globalSize_ % comm_->getSize();
        displs[0] = 0;
        for (int i = 0; i < comm_->getSize(); ++i) {
            sendCounts[i] = globalSize_ / comm_->getSize();
            if (i < remainder) sendCounts[i]++;
            if (i > 0) displs[i] = displs[i-1] + sendCounts[i-1];
        }
        
        comm_->scatter(globalVector.data(), sendCounts.data(), displs.data(), 
                      localData_.data(), localSize_, MPI_DOUBLE, 0);
    } else {
        comm_->scatter(nullptr, nullptr, nullptr, 
                      localData_.data(), localSize_, MPI_DOUBLE, 0);
    }
}

void DistributedVector::allGather(std::vector<double>& globalVector) const {
    if (!comm_) {
        globalVector = localData_;
        return;
    }
    
    globalVector.resize(globalSize_);
    
    // 收集所有进程的本地大小
    std::vector<int> localSizes(comm_->getSize());
    comm_->allgather(&localSize_, 1, MPI_INT, localSizes.data());
    
    // 计算位移
    std::vector<int> displs(comm_->getSize());
    displs[0] = 0;
    for (int i = 1; i < comm_->getSize(); ++i) {
        displs[i] = displs[i-1] + localSizes[i-1];
    }
    
    // 收集数据
    comm_->allgather(localData_.data(), localSize_, MPI_DOUBLE, 
                    globalVector.data(), localSizes.data(), displs.data());
}

void DistributedVector::setZero() {
    std::fill(localData_.begin(), localData_.end(), 0.0);
}

void DistributedVector::scale(double factor) {
    for (auto& value : localData_) {
        value *= factor;
    }
}

double DistributedVector::norm() const {
    double localNorm = 0.0;
    for (double value : localData_) {
        localNorm += value * value;
    }
    
    if (comm_) {
        double globalNorm;
        comm_->allreduce(&localNorm, &globalNorm, 1, MPI_DOUBLE, MPI_SUM);
        return std::sqrt(globalNorm);
    }
    
    return std::sqrt(localNorm);
}

double DistributedVector::dot(const DistributedVector& other) const {
    if (localSize_ != other.localSize_ || offset_ != other.offset_) {
        throw std::runtime_error("Vector size mismatch in dot product");
    }
    
    double localDot = 0.0;
    for (int i = 0; i < localSize_; ++i) {
        localDot += localData_[i] * other.localData_[i];
    }
    
    if (comm_) {
        double globalDot;
        comm_->allreduce(&localDot, &globalDot, 1, MPI_DOUBLE, MPI_SUM);
        return globalDot;
    }
    
    return localDot;
}

// ===== DistributedMatrix实现 =====
DistributedMatrix::DistributedMatrix(int globalRows, int globalCols, 
                                   std::shared_ptr<MPICommunicator> comm)
    : comm_(comm ? comm : MPIUtils::getDefaultComm()), 
      globalRows_(globalRows), globalCols_(globalCols) {
    
    if (!comm_) {
        // 串行模式
        localRows_ = globalRows;
        localCols_ = globalCols;
        rowOffset_ = 0;
        colOffset_ = 0;
    } else {
        // 计算本地大小和偏移
        int rank = comm_->getRank();
        int size = comm_->getSize();
        
        localRows_ = globalRows / size;
        int remainder = globalRows % size;
        
        if (rank < remainder) {
            localRows_++;
            rowOffset_ = rank * localRows_;
        } else {
            rowOffset_ = rank * localRows_ + remainder;
        }
        
        // 列分布（假设列也是均匀分布）
        localCols_ = globalCols / size;
        remainder = globalCols % size;
        
        if (rank < remainder) {
            localCols_++;
            colOffset_ = rank * localCols_;
        } else {
            colOffset_ = rank * localCols_ + remainder;
        }
    }
    
    // 创建本地矩阵
    localMatrix_ = std::make_shared<CRSMatrix>(localRows_, localCols_);
}

DistributedMatrix::DistributedMatrix(int localRows, int localCols, int globalRows, int globalCols,
                                   int rowOffset, int colOffset, 
                                   std::shared_ptr<MPICommunicator> comm)
    : comm_(comm ? comm : MPIUtils::getDefaultComm()),
      localRows_(localRows), globalRows_(globalRows),
      localCols_(localCols), globalCols_(globalCols),
      rowOffset_(rowOffset), colOffset_(colOffset) {
    
    localMatrix_ = std::make_shared<CRSMatrix>(localRows_, localCols_);
}

void DistributedMatrix::SetElement(Integer i, Integer j, Real value) {
    if (i < rowOffset_ || i >= rowOffset_ + localRows_ ||
        j < colOffset_ || j >= colOffset_ + localCols_) {
        throw std::out_of_range("DistributedMatrix::SetElement: indices out of local range");
    }
    
    localMatrix_->SetElement(i - rowOffset_, j - colOffset_, value);
}

Real DistributedMatrix::GetElement(Integer i, Integer j) const {
    if (i >= rowOffset_ && i < rowOffset_ + localRows_ &&
        j >= colOffset_ && j < colOffset_ + localCols_) {
        return localMatrix_->GetElement(i - rowOffset_, j - colOffset_);
    }
    throw std::out_of_range("DistributedMatrix::GetElement: indices not found");
}

void DistributedMatrix::AddToElement(Integer i, Integer j, Real value) {
    if (i < rowOffset_ || i >= rowOffset_ + localRows_ ||
        j < colOffset_ || j >= colOffset_ + localCols_) {
        throw std::out_of_range("DistributedMatrix::AddToElement: indices out of local range");
    }
    
    localMatrix_->AddToElement(i - rowOffset_, j - colOffset_, value);
}

Integer DistributedMatrix::GetRows() const {
    return globalRows_;
}

Integer DistributedMatrix::GetCols() const {
    return globalCols_;
}

void DistributedMatrix::multiply(const DistributedVector& x, DistributedVector& y) const {
    if (globalCols_ != x.Size() || globalRows_ != y.Size()) {
        throw std::runtime_error("Matrix-vector multiplication size mismatch");
    }
    
    // 本地矩阵-向量乘法
    for (int i = 0; i < localRows_; ++i) {
        double sum = 0.0;
        for (int j = 0; j < localCols_; ++j) {
            sum += localMatrix_->GetElement(i, j) * x.GetElement(colOffset_ + j);
        }
        y.SetElement(rowOffset_ + i, sum);
    }
    
    // 组装重叠区域（简化实现）
    y.updateGhostData();
}

void DistributedMatrix::assemble() {
    // 简化实现：在实际应用中需要实现重叠区域的通信
    // 这里只是占位符实现
    if (comm_) {
        comm_->barrier();
    }
}

// ===== DistributedLinearSystem实现 =====
DistributedLinearSystem::DistributedLinearSystem(int globalSize, std::shared_ptr<MPICommunicator> comm)
    : comm_(comm ? comm : MPIUtils::getDefaultComm()) {
    
    A_ = std::make_shared<DistributedMatrix>(globalSize, globalSize, comm_);
    b_ = std::make_shared<DistributedVector>(globalSize, comm_);
    x_ = std::make_shared<DistributedVector>(globalSize, comm_);
}

DistributedLinearSystem::DistributedLinearSystem(int localSize, int globalSize, int offset,
                                               std::shared_ptr<MPICommunicator> comm)
    : comm_(comm ? comm : MPIUtils::getDefaultComm()) {
    
    A_ = std::make_shared<DistributedMatrix>(localSize, localSize, globalSize, globalSize,
                                            offset, offset, comm_);
    b_ = std::make_shared<DistributedVector>(localSize, globalSize, offset, comm_);
    x_ = std::make_shared<DistributedVector>(localSize, globalSize, offset, comm_);
}

void DistributedLinearSystem::setZero() {
    A_->getLocalMatrix()->SetZero();
    b_->setZero();
    x_->setZero();
}

void DistributedLinearSystem::assemble() {
    A_->assemble();
    // 在实际应用中，还需要组装右端向量
}

// ===== DistributedLinearAlgebraUtils实现 =====
namespace DistributedLinearAlgebraUtils {
    
std::shared_ptr<DistributedVector> createUniformDistributedVector(int globalSize, 
                                                                 std::shared_ptr<MPICommunicator> comm) {
    return std::make_shared<DistributedVector>(globalSize, comm);
}

std::shared_ptr<DistributedMatrix> createUniformDistributedMatrix(int globalRows, int globalCols,
                                                                 std::shared_ptr<MPICommunicator> comm) {
    return std::make_shared<DistributedMatrix>(globalRows, globalCols, comm);
}

std::shared_ptr<DistributedLinearSystem> createUniformDistributedLinearSystem(int globalSize,
                                                                             std::shared_ptr<MPICommunicator> comm) {
    return std::make_shared<DistributedLinearSystem>(globalSize, comm);
}

double globalNorm(const DistributedVector& vec) {
    return vec.norm();
}

double globalDot(const DistributedVector& vec1, const DistributedVector& vec2) {
    return vec1.dot(vec2);
}

void globalCopy(const DistributedVector& src, DistributedVector& dest) {
    if (src.getGlobalSize() != dest.getGlobalSize()) {
        throw std::runtime_error("Vector size mismatch in global copy");
    }
    
    // 简化实现：通过主进程进行复制
    if (src.getCommunicator() && src.getCommunicator()->isMaster()) {
        std::vector<double> globalData;
        src.allGather(globalData);
        dest.scatterFromMaster(globalData);
    }
}

void globalScale(DistributedVector& vec, double factor) {
    vec.scale(factor);
}

void globalAdd(DistributedVector& vec1, const DistributedVector& vec2) {
    if (vec1.getLocalSize() != vec2.getLocalSize()) {
        throw std::runtime_error("Vector size mismatch in global add");
    }
    
    for (int i = 0; i < vec1.getLocalSize(); ++i) {
        double value = vec1.getLocalData()[i] + vec2.getLocalData()[i];
        vec1.SetElement(vec1.getOffset() + i, value);
    }
}

void globalAddScaled(DistributedVector& vec1, double alpha, const DistributedVector& vec2) {
    if (vec1.getLocalSize() != vec2.getLocalSize()) {
        throw std::runtime_error("Vector size mismatch in global add scaled");
    }
    
    for (int i = 0; i < vec1.getLocalSize(); ++i) {
        double value = vec1.getLocalData()[i] + alpha * vec2.getLocalData()[i];
        vec1.SetElement(vec1.getOffset() + i, value);
    }
}

bool checkGlobalConvergence(const DistributedVector& residual, double tolerance) {
    double norm = residual.norm();
    
    if (residual.getCommunicator()) {
        // 所有进程检查收敛性
        bool localConverged = (norm < tolerance);
        bool globalConverged;
        
        residual.getCommunicator()->allreduce(&localConverged, &globalConverged, 
                                             1, MPI_C_BOOL, MPI_LAND);
        return globalConverged;
    }
    
    return norm < tolerance;
}

} // namespace DistributedLinearAlgebraUtils

} // namespace elmer