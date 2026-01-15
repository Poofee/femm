/**
 * @file LoggerFactory.cpp
 * @brief 日志工厂类实现
 * 
 * 负责创建和管理日志器实例，实现日志工厂的具体功能。
 */

#include "LoggerFactory.h"
#include "SpdlogAdapter.h"
#include <memory>
#include <unordered_map>
#include <mutex>

namespace elmer {

namespace {
    // 默认日志级别
    LogLevel g_defaultLevel = LogLevel::INFO;
    
    // 默认日志文件名
    std::string g_logFile;
    
    // 默认日志器实例
    std::shared_ptr<ILogger> g_defaultLogger;
    
    // 日志器注册表
    std::unordered_map<std::string, std::shared_ptr<ILogger>> g_loggers;
    
    // 互斥锁用于线程安全
    std::mutex g_loggerMutex;
}

std::shared_ptr<ILogger> LoggerFactory::getDefaultLogger() {
    std::lock_guard<std::mutex> lock(g_loggerMutex);
    
    if (!g_defaultLogger) {
        // 创建默认日志器
        g_defaultLogger = std::make_shared<SpdlogAdapter>("default");
        g_defaultLogger->setLevel(g_defaultLevel);
        
        // 如果设置了日志文件，添加文件输出
        if (!g_logFile.empty()) {
            auto spdlogAdapter = std::dynamic_pointer_cast<SpdlogAdapter>(g_defaultLogger);
            if (spdlogAdapter) {
                spdlogAdapter->addFileSink(g_logFile);
            }
        }
        
        // 注册到日志器表
        g_loggers["default"] = g_defaultLogger;
    }
    
    return g_defaultLogger;
}

std::shared_ptr<ILogger> LoggerFactory::createLogger(const std::string& name) {
    std::lock_guard<std::mutex> lock(g_loggerMutex);
    
    // 检查是否已存在同名日志器
    auto it = g_loggers.find(name);
    if (it != g_loggers.end()) {
        return it->second;
    }
    
    // 创建新日志器
    auto logger = std::make_shared<SpdlogAdapter>(name);
    logger->setLevel(g_defaultLevel);
    
    // 如果设置了日志文件，添加文件输出
    if (!g_logFile.empty()) {
        auto spdlogAdapter = std::dynamic_pointer_cast<SpdlogAdapter>(logger);
        if (spdlogAdapter) {
            spdlogAdapter->addFileSink(g_logFile);
        }
    }
    
    // 注册到日志器表
    g_loggers[name] = logger;
    
    return logger;
}

void LoggerFactory::setDefaultLevel(LogLevel level) {
    std::lock_guard<std::mutex> lock(g_loggerMutex);
    
    g_defaultLevel = level;
    
    // 更新所有已创建的日志器
    for (auto& pair : g_loggers) {
        pair.second->setLevel(level);
    }
}

void LoggerFactory::setLogFile(const std::string& filename) {
    std::lock_guard<std::mutex> lock(g_loggerMutex);
    
    g_logFile = filename;
    
    // 为所有已创建的日志器添加文件输出
    for (auto& pair : g_loggers) {
        auto spdlogAdapter = std::dynamic_pointer_cast<SpdlogAdapter>(pair.second);
        if (spdlogAdapter) {
            spdlogAdapter->addFileSink(filename);
        }
    }
}

void LoggerFactory::flushAll() {
    std::lock_guard<std::mutex> lock(g_loggerMutex);
    
    for (auto& pair : g_loggers) {
        pair.second->flush();
    }
}

} // namespace elmer