/**
 * @file SpdlogAdapter.h
 * @brief spdlog适配器，实现ILogger接口
 * 
 * 将spdlog的具体实现隐藏在ILogger接口后面。
 */

#pragma once

// Windows头文件保护，避免ERROR宏冲突
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef ERROR
#endif

#include "LoggerInterface.h"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <memory>

namespace elmer {

/**
 * @brief spdlog适配器类
 * 
 * 实现ILogger接口，使用spdlog作为底层日志库。
 */
class SpdlogAdapter : public ILogger {
private:
    std::shared_ptr<spdlog::logger> logger_;
    
public:
    /**
     * @brief 构造函数
     * @param name 日志器名称
     */
    explicit SpdlogAdapter(const std::string& name) {
        // 创建控制台输出
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%n] %v");
        
        // 创建日志器
        logger_ = std::make_shared<spdlog::logger>(name, console_sink);
        
        // 设置默认级别
        logger_->set_level(spdlog::level::info);
        
        // 注册到spdlog
        spdlog::register_logger(logger_);
    }
    
    /**
     * @brief 析构函数
     */
    ~SpdlogAdapter() override = default;
    
    /**
     * @brief 设置日志级别
     * @param level 日志级别
     */
    void setLevel(LogLevel level) override {
        spdlog::level::level_enum spdlog_level;
        
        switch (level) {
            case LogLevel::TRACE:
                spdlog_level = spdlog::level::trace;
                break;
            case LogLevel::DEBUG:
                spdlog_level = spdlog::level::debug;
                break;
            case LogLevel::INFO:
                spdlog_level = spdlog::level::info;
                break;
            case LogLevel::WARN:
                spdlog_level = spdlog::level::warn;
                break;
            case LogLevel::ERROR:
                spdlog_level = spdlog::level::err;
                break;
            case LogLevel::CRITICAL:
                spdlog_level = spdlog::level::critical;
                break;
            default:
                spdlog_level = spdlog::level::info;
                break;
        }
        
        logger_->set_level(spdlog_level);
    }
    
    /**
     * @brief 获取当前日志级别
     * @return 当前日志级别
     */
    LogLevel getLevel() const override {
        auto spdlog_level = logger_->level();
        
        switch (spdlog_level) {
            case spdlog::level::trace:
                return LogLevel::TRACE;
            case spdlog::level::debug:
                return LogLevel::DEBUG;
            case spdlog::level::info:
                return LogLevel::INFO;
            case spdlog::level::warn:
                return LogLevel::WARN;
            case spdlog::level::err:
                return LogLevel::ERROR;
            case spdlog::level::critical:
                return LogLevel::CRITICAL;
            default:
                return LogLevel::INFO;
        }
    }
    
    /**
     * @brief 跟踪级别日志
     * @param message 日志消息
     */
    void trace(const std::string& message) override {
        logger_->trace(message);
    }
    
    /**
     * @brief 调试级别日志
     * @param message 日志消息
     */
    void debug(const std::string& message) override {
        logger_->debug(message);
    }
    
    /**
     * @brief 信息级别日志
     * @param message 日志消息
     */
    void info(const std::string& message) override {
        logger_->info(message);
    }
    
    /**
     * @brief 警告级别日志
     * @param message 日志消息
     */
    void warn(const std::string& message) override {
        logger_->warn(message);
    }
    
    /**
     * @brief 错误级别日志
     * @param message 日志消息
     */
    void error(const std::string& message) override {
        logger_->error(message);
    }
    
    /**
     * @brief 严重级别日志
     * @param message 日志消息
     */
    void critical(const std::string& message) override {
        logger_->critical(message);
    }
    
    /**
     * @brief 刷新日志缓冲区
     */
    void flush() override {
        logger_->flush();
    }
    
    /**
     * @brief 检查是否启用指定日志级别
     * @param level 要检查的日志级别
     * @return 是否启用
     */
    bool shouldLog(LogLevel level) const override {
        auto current_level = getLevel();
        
        // 日志级别从低到高：TRACE < DEBUG < INFO < WARN < ERROR < CRITICAL
        return static_cast<int>(level) >= static_cast<int>(current_level);
    }
    
    /**
     * @brief 添加文件输出
     * @param filename 日志文件名
     */
    void addFileSink(const std::string& filename) {
        try {
            auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename);
            file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] [%n] %v");
            logger_->sinks().push_back(file_sink);
        } catch (const spdlog::spdlog_ex& ex) {
            logger_->error("无法创建文件输出: {}", ex.what());
        }
    }
    
    /**
     * @brief 获取底层spdlog日志器
     * @return spdlog日志器
     */
    std::shared_ptr<spdlog::logger> getSpdlogLogger() const {
        return logger_;
    }
};

} // namespace elmer