/**
 * @file MPICommunicator.h
 * @brief MPI通信器封装类
 * 
 * 提供统一的MPI通信接口，支持串行和并行模式
 */

#pragma once

#include <iostream>
#include <vector>
#include <memory>

#ifdef USE_MPI
#include <mpi.h>
#else
// MPI模拟类型定义（当MPI不可用时）
typedef int MPI_Datatype;
typedef int MPI_Op;
#define MPI_INT 1
#define MPI_FLOAT 2
#define MPI_DOUBLE 3
#define MPI_LONG 4
#define MPI_UNSIGNED 5
#define MPI_BYTE 6
#define MPI_SUM 0
#define MPI_MAX 1
#define MPI_MIN 2
#define MPI_PROD 3
#define MPI_OP_NULL -1
#endif

namespace elmer {

/**
 * @brief MPI通信器类
 * 
 * 封装MPI通信功能，支持串行和并行模式
 */
class MPICommunicator {
private:
#ifdef USE_MPI
    MPI_Comm comm_;                    ///< MPI通信器
#endif
    int rank_ = 0;                     ///< 当前进程编号
    int size_ = 1;                     ///< 进程总数
    bool initialized_ = false;         ///< 是否已初始化
    bool useMPI_ = false;              ///< 是否使用MPI
    
public:
    /**
     * @brief 构造函数
     */
    MPICommunicator() = default;
    
    /**
     * @brief 析构函数
     */
    ~MPICommunicator() {
        cleanup();
    }
    
    /**
     * @brief 初始化MPI
     * @return 是否初始化成功
     */
    bool initialize() {
#ifdef USE_MPI
        int mpiInitialized = 0;
        MPI_Initialized(&mpiInitialized);
        
        if (!mpiInitialized) {
            if (MPI_Init(nullptr, nullptr) != MPI_SUCCESS) {
                std::cerr << "错误：MPI初始化失败" << std::endl;
                return false;
            }
        }
        
        comm_ = MPI_COMM_WORLD;
        MPI_Comm_rank(comm_, &rank_);
        MPI_Comm_size(comm_, &size_);
        
        useMPI_ = true;
        initialized_ = true;
        
        return true;
#else
        std::cerr << "警告：编译时未启用MPI支持" << std::endl;
        return false;
#endif
    }
    
    /**
     * @brief 初始化串行模式
     */
    void initializeSerial() {
        rank_ = 0;
        size_ = 1;
        useMPI_ = false;
        initialized_ = true;
    }
    
    /**
     * @brief 清理资源
     */
    void cleanup() {
#ifdef USE_MPI
        if (useMPI_ && initialized_) {
            int mpiFinalized = 0;
            MPI_Finalized(&mpiFinalized);
            
            if (!mpiFinalized) {
                MPI_Finalize();
            }
        }
#endif
        initialized_ = false;
    }
    
    /**
     * @brief 获取当前进程编号
     * @return 进程编号
     */
    int getRank() const {
        return rank_;
    }
    
    /**
     * @brief 获取进程总数
     * @return 进程总数
     */
    int getSize() const {
        return size_;
    }
    
    /**
     * @brief 检查是否为主进程
     * @return 是否为主进程
     */
    bool isMaster() const {
        return rank_ == 0;
    }
    
    /**
     * @brief 检查是否使用MPI
     * @return 是否使用MPI
     */
    bool isUsingMPI() const {
        return useMPI_;
    }
    
    /**
     * @brief 检查是否已初始化
     * @return 是否已初始化
     */
    bool isInitialized() const {
        return initialized_;
    }
    
    /**
     * @brief 同步所有进程
     */
    void barrier() const {
#ifdef USE_MPI
        if (useMPI_ && initialized_) {
            MPI_Barrier(comm_);
        }
#endif
    }
    
    /**
     * @brief 广播数据
     * @tparam T 数据类型
     * @param data 数据指针
     * @param count 数据数量
     * @param root 根进程编号
     */
    template<typename T>
    void broadcast(T* data, int count, int root = 0) const {
#ifdef USE_MPI
        if (useMPI_ && initialized_) {
            MPI_Datatype mpiType = getMPIType<T>();
            MPI_Bcast(data, count, mpiType, root, comm_);
        }
#endif
    }
    
    /**
     * @brief 收集数据
     * @tparam T 数据类型
     * @param sendData 发送数据
     * @param recvData 接收数据
     * @param count 每个进程发送的数据数量
     * @param root 根进程编号
     */
    template<typename T>
    void gather(const T* sendData, T* recvData, int count, int root = 0) const {
#ifdef USE_MPI
        if (useMPI_ && initialized_) {
            MPI_Datatype mpiType = getMPIType<T>();
            MPI_Gather(sendData, count, mpiType, recvData, count, mpiType, root, comm_);
        } else {
            // 串行模式，直接复制数据
            if (rank_ == root && sendData != recvData) {
                std::copy(sendData, sendData + count, recvData);
            }
        }
#endif
    }
    
    /**
     * @brief 分发数据
     * @tparam T 数据类型
     * @param sendData 发送数据
     * @param recvData 接收数据
     * @param count 每个进程接收的数据数量
     * @param root 根进程编号
     */
    template<typename T>
    void scatter(const T* sendData, T* recvData, int count, int root = 0) const {
#ifdef USE_MPI
        if (useMPI_ && initialized_) {
            MPI_Datatype mpiType = getMPIType<T>();
            MPI_Scatter(sendData, count, mpiType, recvData, count, mpiType, root, comm_);
        } else {
            // 串行模式，直接复制数据
            if (rank_ == root && sendData != recvData) {
                std::copy(sendData, sendData + count, recvData);
            }
        }
#endif
    }
    
    /**
     * @brief 归约操作
     * @tparam T 数据类型
     * @param sendData 发送数据
     * @param recvData 接收数据
     * @param count 数据数量
     * @param op 归约操作类型
     * @param root 根进程编号
     */
    template<typename T>
    void reduce(const T* sendData, T* recvData, int count, const std::string& op = "sum", int root = 0) const {
#ifdef USE_MPI
        if (useMPI_ && initialized_) {
            MPI_Datatype mpiType = getMPIType<T>();
            MPI_Op mpiOp = getMPIOp(op);
            MPI_Reduce(sendData, recvData, count, mpiType, mpiOp, root, comm_);
        } else {
            // 串行模式，直接复制数据
            if (rank_ == root && sendData != recvData) {
                std::copy(sendData, sendData + count, recvData);
            }
        }
#endif
    }
    
    /**
     * @brief 全局归约操作
     * @tparam T 数据类型
     * @param sendData 发送数据
     * @param recvData 接收数据
     * @param count 数据数量
     * @param op 归约操作类型
     */
    template<typename T>
    void allReduce(const T* sendData, T* recvData, int count, const std::string& op = "sum") const {
#ifdef USE_MPI
        if (useMPI_ && initialized_) {
            MPI_Datatype mpiType = getMPIType<T>();
            MPI_Op mpiOp = getMPIOp(op);
            MPI_Allreduce(sendData, recvData, count, mpiType, mpiOp, comm_);
        } else {
            // 串行模式，直接复制数据
            if (sendData != recvData) {
                std::copy(sendData, sendData + count, recvData);
            }
        }
#endif
    }
    
private:
    /**
     * @brief 获取MPI数据类型
     * @tparam T C++数据类型
     * @return MPI数据类型
     */
    template<typename T>
    static MPI_Datatype getMPIType() {
#ifdef USE_MPI
        if constexpr (std::is_same_v<T, int>) {
            return MPI_INT;
        } else if constexpr (std::is_same_v<T, float>) {
            return MPI_FLOAT;
        } else if constexpr (std::is_same_v<T, double>) {
            return MPI_DOUBLE;
        } else if constexpr (std::is_same_v<T, long>) {
            return MPI_LONG;
        } else if constexpr (std::is_same_v<T, unsigned>) {
            return MPI_UNSIGNED;
        } else {
            return MPI_BYTE;
        }
#else
        return MPI_DATATYPE_NULL;
#endif
    }
    
    /**
     * @brief 获取MPI归约操作
     * @param op 操作名称
     * @return MPI操作
     */
    static MPI_Op getMPIOp(const std::string& op) {
#ifdef USE_MPI
        if (op == "sum") {
            return MPI_SUM;
        } else if (op == "max") {
            return MPI_MAX;
        } else if (op == "min") {
            return MPI_MIN;
        } else if (op == "prod") {
            return MPI_PROD;
        } else {
            return MPI_SUM;
        }
#else
        return MPI_OP_NULL;
#endif
    }
};

} // namespace elmer