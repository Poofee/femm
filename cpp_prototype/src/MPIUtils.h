/**
 * @file MPIUtils.h
 * @brief MPI工具函数
 * 
 * 提供MPI相关的工具函数和宏定义
 */

#pragma once

#include <iostream>
#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace elmer {

/**
 * @brief MPI工具类
 */
class MPIUtils {
public:
    /**
     * @brief 检查是否启用MPI支持
     * @return 是否启用MPI
     */
    static bool isMPIEnabled() {
#ifdef USE_MPI
        return true;
#else
        return false;
#endif
    }
    
    /**
     * @brief 获取MPI版本信息
     * @return MPI版本字符串
     */
    static std::string getMPIVersion() {
#ifdef USE_MPI
        int version, subversion;
        MPI_Get_version(&version, &subversion);
        return std::to_string(version) + "." + std::to_string(subversion);
#else
        return "MPI not enabled";
#endif
    }
    
    /**
     * @brief 打印MPI环境信息
     */
    static void printMPIInfo() {
#ifdef USE_MPI
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        
        if (rank == 0) {
            std::cout << "MPI环境信息:" << std::endl;
            std::cout << "  MPI版本: " << getMPIVersion() << std::endl;
            std::cout << "  进程总数: " << size << std::endl;
            std::cout << "  当前进程: " << rank << std::endl;
        }
#else
        std::cout << "MPI未启用，使用串行模式" << std::endl;
#endif
    }
    
    /**
     * @brief 检查MPI错误
     * @param errorCode MPI错误代码
     * @param functionName 函数名称
     * @return 是否成功
     */
    static bool checkMPIError(int errorCode, const std::string& functionName) {
#ifdef USE_MPI
        if (errorCode != MPI_SUCCESS) {
            char errorString[MPI_MAX_ERROR_STRING];
            int resultLen;
            MPI_Error_string(errorCode, errorString, &resultLen);
            
            std::cerr << "MPI错误在 " << functionName << ": " << errorString << std::endl;
            return false;
        }
#endif
        return true;
    }
    
    /**
     * @brief 获取当前时间（用于性能测量）
     * @return 当前时间（秒）
     */
    static double getTime() {
#ifdef USE_MPI
        return MPI_Wtime();
#else
        // 使用系统时钟作为替代
        static auto startTime = std::chrono::high_resolution_clock::now();
        auto currentTime = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(currentTime - startTime).count();
#endif
    }
    
    /**
     * @brief 同步所有进程并测量时间
     * @param comm MPI通信器
     * @return 同步时间（秒）
     */
    static double barrierAndTime() {
#ifdef USE_MPI
        double startTime = getTime();
        MPI_Barrier(MPI_COMM_WORLD);
        double endTime = getTime();
        return endTime - startTime;
#else
        return 0.0;
#endif
    }
    
    /**
     * @brief 获取进程名称
     * @return 进程名称字符串
     */
    static std::string getProcessorName() {
#ifdef USE_MPI
        char name[MPI_MAX_PROCESSOR_NAME];
        int nameLen;
        MPI_Get_processor_name(name, &nameLen);
        return std::string(name, nameLen);
#else
        return "localhost";
#endif
    }
    
    /**
     * @brief 打印进程分布信息
     */
    static void printProcessorDistribution() {
#ifdef USE_MPI
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        
        // 收集所有进程的处理器名称
        std::vector<std::string> processorNames(size);
        std::string myName = getProcessorName();
        
        // 主进程收集所有信息
        if (rank == 0) {
            processorNames[0] = myName;
            for (int i = 1; i < size; i++) {
                char buffer[MPI_MAX_PROCESSOR_NAME];
                MPI_Recv(buffer, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                processorNames[i] = std::string(buffer);
            }
            
            // 打印分布信息
            std::cout << "进程分布:" << std::endl;
            for (int i = 0; i < size; i++) {
                std::cout << "  进程 " << i << ": " << processorNames[i] << std::endl;
            }
        } else {
            // 从进程发送信息
            MPI_Send(myName.c_str(), myName.length() + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
#else
        std::cout << "进程分布: 单进程运行于 " << getProcessorName() << std::endl;
#endif
    }
};

/**
 * @brief MPI范围保护类
 * 
 * 用于自动管理MPI环境的初始化和清理
 */
class MPIScope {
private:
    bool initializedHere_ = false;
    
public:
    /**
     * @brief 构造函数
     * @param initialize 是否初始化MPI
     */
    MPIScope(bool initialize = true) {
#ifdef USE_MPI
        if (initialize) {
            int mpiInitialized;
            MPI_Initialized(&mpiInitialized);
            
            if (!mpiInitialized) {
                MPI_Init(nullptr, nullptr);
                initializedHere_ = true;
            }
        }
#endif
    }
    
    /**
     * @brief 析构函数
     */
    ~MPIScope() {
#ifdef USE_MPI
        if (initializedHere_) {
            int mpiFinalized;
            MPI_Finalized(&mpiFinalized);
            
            if (!mpiFinalized) {
                MPI_Finalize();
            }
        }
#endif
    }
    
    /**
     * @brief 检查是否在此处初始化了MPI
     * @return 是否在此处初始化
     */
    bool initializedHere() const {
        return initializedHere_;
    }
};

/**
 * @brief MPI性能测量类
 */
class MPITimer {
private:
    double startTime_;
    std::string operationName_;
    
public:
    /**
     * @brief 构造函数
     * @param operationName 操作名称
     */
    explicit MPITimer(const std::string& operationName = "")
        : operationName_(operationName) {
        startTime_ = MPIUtils::getTime();
    }
    
    /**
     * @brief 析构函数（自动打印耗时）
     */
    ~MPITimer() {
        double endTime = MPIUtils::getTime();
        double duration = endTime - startTime_;
        
#ifdef USE_MPI
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        if (rank == 0 && !operationName_.empty()) {
            std::cout << operationName_ << " 耗时: " << duration << " 秒" << std::endl;
        }
#else
        if (!operationName_.empty()) {
            std::cout << operationName_ << " 耗时: " << duration << " 秒" << std::endl;
        }
#endif
    }
    
    /**
     * @brief 获取经过的时间
     * @return 经过的时间（秒）
     */
    double elapsed() const {
        return MPIUtils::getTime() - startTime_;
    }
    
    /**
     * @brief 重置计时器
     */
    void reset() {
        startTime_ = MPIUtils::getTime();
    }
};

// MPI相关的宏定义
#ifdef USE_MPI
#define MPI_SAFE_CALL(call) \
    do { \
        int mpiError = (call); \
        if (!MPIUtils::checkMPIError(mpiError, #call)) { \
            MPI_Abort(MPI_COMM_WORLD, mpiError); \
        } \
    } while(0)

#define MPI_MASTER_ONLY(code) \
    do { \
        int rank; \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); \
        if (rank == 0) { \
            code \
        } \
    } while(0)

#define MPI_TIMER(name) MPITimer timer##__LINE__(name)

#else
#define MPI_SAFE_CALL(call) call
#define MPI_MASTER_ONLY(code) code
#define MPI_TIMER(name) MPITimer timer##__LINE__(name)
#endif

} // namespace elmer