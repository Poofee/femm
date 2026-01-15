/**
 * @file LoggerFactory.h
 * @brief 日志工厂类的声明
 * 
 * 负责创建和管理日志器实例。
 */

#pragma once

#include "LoggerInterface.h"
#include <memory>
#include <string>

namespace elmer {

/**
 * @brief 日志工厂类
 * 
 * 负责创建和管理日志器实例。
 */
class LoggerFactory {
public:
    /**
     * @brief 获取默认日志器
     * @return 默认日志器实例
     */
    static std::shared_ptr<ILogger> getDefaultLogger();
    
    /**
     * @brief 创建指定名称的日志器
     * @param name 日志器名称
     * @return 日志器实例
     */
    static std::shared_ptr<ILogger> createLogger(const std::string& name);
    
    /**
     * @brief 设置默认日志级别
     * @param level 默认日志级别
     */
    static void setDefaultLevel(LogLevel level);
    
    /**
     * @brief 设置日志输出文件
     * @param filename 日志文件名
     */
    static void setLogFile(const std::string& filename);
    
    /**
     * @brief 刷新所有日志器
     */
    static void flushAll();
};

/**
 * @brief 全局日志宏
 * 
 * 提供便捷的全局日志访问方式。
 */
#define ELMER_TRACE(...)   elmer::LoggerFactory::getDefaultLogger()->trace(__VA_ARGS__)
#define ELMER_DEBUG(...)   elmer::LoggerFactory::getDefaultLogger()->debug(__VA_ARGS__)
#define ELMER_INFO(...)    elmer::LoggerFactory::getDefaultLogger()->info(__VA_ARGS__)
#define ELMER_WARN(...)    elmer::LoggerFactory::getDefaultLogger()->warn(__VA_ARGS__)
#define ELMER_ERROR(...)   elmer::LoggerFactory::getDefaultLogger()->error(__VA_ARGS__)
#define ELMER_CRITICAL(...) elmer::LoggerFactory::getDefaultLogger()->critical(__VA_ARGS__)

} // namespace elmer