/**
 * @file LoggerInterface.h
 * @brief 日志接口类，隐藏具体日志库的实现
 * 
 * 提供统一的日志接口，支持多种日志级别和格式化输出。
 * 当前实现基于spdlog，但可以轻松切换到其他日志库。
 */

#pragma once

// Windows头文件保护，避免ERROR宏冲突
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef ERROR
#endif

#include <memory>
#include <string>

namespace elmer {

/**
 * @brief 日志级别枚举
 */
enum class LogLevel {
    TRACE,    ///< 跟踪级别，最详细的日志
    DEBUG,    ///< 调试级别，用于调试信息
    INFO,     ///< 信息级别，常规信息
    WARN,     ///< 警告级别，潜在问题
    ERROR,    ///< 错误级别，错误信息
    CRITICAL  ///< 严重级别，严重错误
};

/**
 * @brief 日志接口类
 * 
 * 提供统一的日志接口，隐藏具体日志库的实现细节。
 * 支持多种日志级别和格式化输出。
 */
class ILogger {
public:
    virtual ~ILogger() = default;
    
    /**
     * @brief 设置日志级别
     * @param level 日志级别
     */
    virtual void setLevel(LogLevel level) = 0;
    
    /**
     * @brief 获取当前日志级别
     * @return 当前日志级别
     */
    virtual LogLevel getLevel() const = 0;
    
    /**
     * @brief 跟踪级别日志
     * @param message 日志消息
     */
    virtual void trace(const std::string& message) = 0;
    
    /**
     * @brief 调试级别日志
     * @param message 日志消息
     */
    virtual void debug(const std::string& message) = 0;
    
    /**
     * @brief 信息级别日志
     * @param message 日志消息
     */
    virtual void info(const std::string& message) = 0;
    
    /**
     * @brief 警告级别日志
     * @param message 日志消息
     */
    virtual void warn(const std::string& message) = 0;
    
    /**
     * @brief 错误级别日志
     * @param message 日志消息
     */
    virtual void error(const std::string& message) = 0;
    
    /**
     * @brief 严重级别日志
     * @param message 日志消息
     */
    virtual void critical(const std::string& message) = 0;
    
    /**
     * @brief 格式化跟踪级别日志
     * @param format 格式化字符串
     * @param args 格式化参数
     */
    template<typename... Args>
    void trace(const std::string& format, Args&&... args) {
        trace(formatMessage(format, std::forward<Args>(args)...));
    }
    
    /**
     * @brief 格式化调试级别日志
     * @param format 格式化字符串
     * @param args 格式化参数
     */
    template<typename... Args>
    void debug(const std::string& format, Args&&... args) {
        debug(formatMessage(format, std::forward<Args>(args)...));
    }
    
    /**
     * @brief 格式化信息级别日志
     * @param format 格式化字符串
     * @param args 格式化参数
     */
    template<typename... Args>
    void info(const std::string& format, Args&&... args) {
        info(formatMessage(format, std::forward<Args>(args)...));
    }
    
    /**
     * @brief 格式化警告级别日志
     * @param format 格式化字符串
     * @param args 格式化参数
     */
    template<typename... Args>
    void warn(const std::string& format, Args&&... args) {
        warn(formatMessage(format, std::forward<Args>(args)...));
    }
    
    /**
     * @brief 格式化错误级别日志
     * @param format 格式化字符串
     * @param args 格式化参数
     */
    template<typename... Args>
    void error(const std::string& format, Args&&... args) {
        error(formatMessage(format, std::forward<Args>(args)...));
    }
    
    /**
     * @brief 格式化严重级别日志
     * @param format 格式化字符串
     * @param args 格式化参数
     */
    template<typename... Args>
    void critical(const std::string& format, Args&&... args) {
        critical(formatMessage(format, std::forward<Args>(args)...));
    }
    
    /**
     * @brief 刷新日志缓冲区
     */
    virtual void flush() = 0;
    
    /**
     * @brief 检查是否启用指定日志级别
     * @param level 要检查的日志级别
     * @return 是否启用
     */
    virtual bool shouldLog(LogLevel level) const = 0;
    
private:
    /**
     * @brief 格式化消息
     * @param format 格式化字符串
     * @param args 格式化参数
     * @return 格式化后的消息
     */
    template<typename... Args>
    std::string formatMessage(const std::string& format, Args&&... args) {
        // 简化实现：使用字符串拼接
        // TODO: 实现完整的格式化功能
        std::string result = format;
        // 这里可以集成fmt库或使用其他格式化方法
        return result;
    }
};

} // namespace elmer