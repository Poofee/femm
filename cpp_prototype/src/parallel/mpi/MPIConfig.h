// MPIConfig.h - MPI配置和通信包装器定�?// 对应Fortran模块: MPI.F90, Parallel.F90

#pragma once

#ifdef MPI_CXX_FOUND
#include <mpi.h>
#else
// MPI模拟类型定义（当MPI不可用时�?typedef int MPI_Comm;
typedef int MPI_Datatype;
#define MPI_COMM_WORLD 0
#define MPI_COMM_NULL -1
#define MPI_INT 1
#define MPI_DOUBLE 2
#define MPI_FLOAT 3
#define MPI_CHAR 4
#endif

#include <memory>
#include <vector>
#include <string>
#include <stdexcept>

namespace elmer {

// ===== MPI配置�?=====
class MPIConfig {
private:
    static MPIConfig* instance_;  // 单例实例
    bool initialized_;            // MPI初始化标�?    int rank_;                   // 当前进程排名
    int size_;                   // 进程总数
    
    // 私有构造函�?    MPIConfig();
    
public:
    // 禁止拷贝和赋�?    MPIConfig(const MPIConfig&) = delete;
    MPIConfig& operator=(const MPIConfig&) = delete;
    
    // 获取单例实例
    static MPIConfig& getInstance();
    
    // 初始化MPI环境
    void initialize(int* argc = nullptr, char*** argv = nullptr);
    
    // 清理MPI环境
    void finalize();
    
    // 获取进程信息
    int getRank() const { return rank_; }
    int getSize() const { return size_; }
    bool isMaster() const { return rank_ == 0; }
    
    // 检查MPI是否可用
    bool isMPIAvailable() const { return initialized_; }
    
    // 析构函数
    ~MPIConfig();
};

// ===== MPI通信包装器类 =====
class MPICommunicator {
private:
    int rank_;                   // 当前进程排名
    int size_;                   // 进程总数
    MPI_Comm comm_;              // MPI通信�?    bool ownsComm_;              // 是否拥有通信�?    
public:
    // 构造函�?    MPICommunicator(MPI_Comm comm = MPI_COMM_WORLD, bool ownsComm = false);
    
    // 析构函数
    ~MPICommunicator();
    
    // 获取进程信息
    int getRank() const { return rank_; }
    int getSize() const { return size_; }
    bool isMaster() const { return rank_ == 0; }
    MPI_Comm getComm() const { return comm_; }
    
    // 同步操作
    void barrier() const;
    
    // 点对点通信
    void send(const void* data, int count, MPI_Datatype type, int dest, int tag = 0) const;
    void recv(void* data, int count, MPI_Datatype type, int source, int tag = 0, MPI_Status* status = MPI_STATUS_IGNORE) const;
    
    // 广播操作
    void broadcast(void* data, int count, MPI_Datatype type, int root = 0) const;
    
    // 收集操作
    void gather(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root = 0) const;
    void allgather(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) const;
    
    // 分散操作
    void scatter(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root = 0) const;
    
    // 归约操作
    void reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, MPI_Op op, int root = 0) const;
    void allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, MPI_Op op) const;
    
    // 自定义归约操�?    void sum(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root = 0) const;
    void allsum(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) const;
    
    void max(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root = 0) const;
    void allmax(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) const;
    
    void min(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root = 0) const;
    void allmin(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) const;
    
    // 创建子通信�?    std::shared_ptr<MPICommunicator> split(int color, int key) const;
    
    // 错误处理
    static void checkMPIError(int errorCode, const std::string& operation = "");
};

// ===== MPI数据类型映射 =====
template<typename T>
struct MPI_TypeMap;

template<>
struct MPI_TypeMap<int> {
    static MPI_Datatype getType() { return MPI_INT; }
};

template<>
struct MPI_TypeMap<double> {
    static MPI_Datatype getType() { return MPI_DOUBLE; }
};

template<>
struct MPI_TypeMap<float> {
    static MPI_Datatype getType() { return MPI_FLOAT; }
};

template<>
struct MPI_TypeMap<char> {
    static MPI_Datatype getType() { return MPI_CHAR; }
};

// ===== MPI工具函数 =====
namespace MPIUtils {
    
    // 获取全局MPI配置
    MPIConfig& getConfig();
    
    // 获取默认通信�?    std::shared_ptr<MPICommunicator> getDefaultComm();
    
    // 检查是否在MPI环境中运�?    bool isMPIEnabled();
    
    // 获取进程总数
    int getWorldSize();
    
    // 获取当前进程排名
    int getWorldRank();
    
    // 检查是否是主进�?    bool isMasterProcess();
    
    // 同步所有进�?    void barrier();
    
    // 打印MPI信息
    void printMPIInfo();
}

} // namespace elmer
