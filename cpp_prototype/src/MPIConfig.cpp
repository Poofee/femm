// MPIConfig.cpp - MPI配置和通信包装器实现
// 对应Fortran模块: MPI.F90, Parallel.F90

#include "MPIConfig.h"
#include <iostream>

namespace elmer {

// ===== MPIConfig单例实现 =====
MPIConfig* MPIConfig::instance_ = nullptr;

MPIConfig::MPIConfig() : initialized_(false), rank_(0), size_(1) {
    // 默认情况下，假设没有MPI环境
}

MPIConfig& MPIConfig::getInstance() {
    if (!instance_) {
        instance_ = new MPIConfig();
    }
    return *instance_;
}

void MPIConfig::initialize(int* argc, char*** argv) {
#ifdef MPI_CXX_FOUND
    if (!initialized_) {
        int provided;
        MPI_Init_thread(argc, argv, MPI_THREAD_SINGLE, &provided);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
        initialized_ = true;
        
        std::cout << "MPI initialized: rank " << rank_ << "/" << size_ << std::endl;
    }
#else
    // 如果没有MPI支持，使用单进程模式
    rank_ = 0;
    size_ = 1;
    initialized_ = false;
#endif
}

void MPIConfig::finalize() {
#ifdef MPI_CXX_FOUND
    if (initialized_) {
        MPI_Finalize();
        initialized_ = false;
    }
#endif
}

MPIConfig::~MPIConfig() {
    finalize();
    if (instance_) {
        delete instance_;
        instance_ = nullptr;
    }
}

// ===== MPICommunicator实现 =====
MPICommunicator::MPICommunicator(MPI_Comm comm, bool ownsComm) 
    : comm_(comm), ownsComm_(ownsComm) {
#ifdef MPI_CXX_FOUND
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
#else
    rank_ = 0;
    size_ = 1;
#endif
}

MPICommunicator::~MPICommunicator() {
#ifdef MPI_CXX_FOUND
    if (ownsComm_ && comm_ != MPI_COMM_NULL) {
        MPI_Comm_free(&comm_);
    }
#endif
}

void MPICommunicator::barrier() const {
#ifdef MPI_CXX_FOUND
    MPI_Barrier(comm_);
#endif
}

void MPICommunicator::send(const void* data, int count, MPI_Datatype type, int dest, int tag) const {
#ifdef MPI_CXX_FOUND
    checkMPIError(MPI_Send(const_cast<void*>(data), count, type, dest, tag, comm_), "MPI_Send");
#endif
}

void MPICommunicator::recv(void* data, int count, MPI_Datatype type, int source, int tag, MPI_Status* status) const {
#ifdef MPI_CXX_FOUND
    checkMPIError(MPI_Recv(data, count, type, source, tag, comm_, status), "MPI_Recv");
#endif
}

void MPICommunicator::broadcast(void* data, int count, MPI_Datatype type, int root) const {
#ifdef MPI_CXX_FOUND
    checkMPIError(MPI_Bcast(data, count, type, root, comm_), "MPI_Bcast");
#endif
}

void MPICommunicator::gather(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root) const {
#ifdef MPI_CXX_FOUND
    checkMPIError(MPI_Gather(const_cast<void*>(sendbuf), count, type, recvbuf, count, type, root, comm_), "MPI_Gather");
#endif
}

void MPICommunicator::allgather(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) const {
#ifdef MPI_CXX_FOUND
    checkMPIError(MPI_Allgather(const_cast<void*>(sendbuf), count, type, recvbuf, count, type, comm_), "MPI_Allgather");
#endif
}

void MPICommunicator::scatter(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root) const {
#ifdef MPI_CXX_FOUND
    checkMPIError(MPI_Scatter(const_cast<void*>(sendbuf), count, type, recvbuf, count, type, root, comm_), "MPI_Scatter");
#endif
}

void MPICommunicator::reduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, MPI_Op op, int root) const {
#ifdef MPI_CXX_FOUND
    checkMPIError(MPI_Reduce(const_cast<void*>(sendbuf), recvbuf, count, type, op, root, comm_), "MPI_Reduce");
#endif
}

void MPICommunicator::allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, MPI_Op op) const {
#ifdef MPI_CXX_FOUND
    checkMPIError(MPI_Allreduce(const_cast<void*>(sendbuf), recvbuf, count, type, op, comm_), "MPI_Allreduce");
#endif
}

void MPICommunicator::sum(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root) const {
    reduce(sendbuf, recvbuf, count, type, MPI_SUM, root);
}

void MPICommunicator::allsum(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) const {
    allreduce(sendbuf, recvbuf, count, type, MPI_SUM);
}

void MPICommunicator::max(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root) const {
    reduce(sendbuf, recvbuf, count, type, MPI_MAX, root);
}

void MPICommunicator::allmax(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) const {
    allreduce(sendbuf, recvbuf, count, type, MPI_MAX);
}

void MPICommunicator::min(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type, int root) const {
    reduce(sendbuf, recvbuf, count, type, MPI_MIN, root);
}

void MPICommunicator::allmin(const void* sendbuf, void* recvbuf, int count, MPI_Datatype type) const {
    allreduce(sendbuf, recvbuf, count, type, MPI_MIN);
}

std::shared_ptr<MPICommunicator> MPICommunicator::split(int color, int key) const {
#ifdef MPI_CXX_FOUND
    MPI_Comm newComm;
    checkMPIError(MPI_Comm_split(comm_, color, key, &newComm), "MPI_Comm_split");
    return std::make_shared<MPICommunicator>(newComm, true);
#else
    return std::make_shared<MPICommunicator>(MPI_COMM_WORLD, false);
#endif
}

void MPICommunicator::checkMPIError(int errorCode, const std::string& operation) {
#ifdef MPI_CXX_FOUND
    if (errorCode != MPI_SUCCESS) {
        char errorString[MPI_MAX_ERROR_STRING];
        int resultLen;
        MPI_Error_string(errorCode, errorString, &resultLen);
        throw std::runtime_error("MPI error in " + operation + ": " + std::string(errorString));
    }
#endif
}

// ===== MPIUtils实现 =====
namespace MPIUtils {
    
MPIConfig& getConfig() {
    return MPIConfig::getInstance();
}

std::shared_ptr<MPICommunicator> getDefaultComm() {
    return std::make_shared<MPICommunicator>(MPI_COMM_WORLD, false);
}

bool isMPIEnabled() {
#ifdef MPI_CXX_FOUND
    return getConfig().isMPIAvailable();
#else
    return false;
#endif
}

int getWorldSize() {
    return getConfig().getSize();
}

int getWorldRank() {
    return getConfig().getRank();
}

bool isMasterProcess() {
    return getConfig().isMaster();
}

void barrier() {
    if (isMPIEnabled()) {
        getDefaultComm()->barrier();
    }
}

void printMPIInfo() {
    if (isMPIEnabled()) {
        auto& config = getConfig();
        std::cout << "MPI Info: " << config.getSize() << " processes, current rank: " 
                  << config.getRank() << std::endl;
    } else {
        std::cout << "MPI not available - running in serial mode" << std::endl;
    }
}

} // namespace MPIUtils

} // namespace elmer