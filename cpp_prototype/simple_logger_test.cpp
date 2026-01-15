/**
 * @file simple_logger_test.cpp
 * @brief 简化版日志接口测试
 * 
 * 使用最简化的方式测试日志接口功能，避免复杂的构建依赖。
 */

#include <iostream>
#include <string>

// 模拟日志接口的简化版本
class SimpleLogger {
public:
    enum class Level { TRACE, DEBUG, INFO, WARN, ERROR, CRITICAL };
    
    SimpleLogger(const std::string& name) : name_(name), level_(Level::INFO) {}
    
    void setLevel(Level level) { level_ = level; }
    Level getLevel() const { return level_; }
    
    void trace(const std::string& message) { log(Level::TRACE, message); }
    void debug(const std::string& message) { log(Level::DEBUG, message); }
    void info(const std::string& message) { log(Level::INFO, message); }
    void warn(const std::string& message) { log(Level::WARN, message); }
    void error(const std::string& message) { log(Level::ERROR, message); }
    void critical(const std::string& message) { log(Level::CRITICAL, message); }
    
    template<typename... Args>
    void trace(const std::string& format, Args&&... args) {
        trace(formatMessage(format, std::forward<Args>(args)...));
    }
    
    template<typename... Args>
    void debug(const std::string& format, Args&&... args) {
        debug(formatMessage(format, std::forward<Args>(args)...));
    }
    
    template<typename... Args>
    void info(const std::string& format, Args&&... args) {
        info(formatMessage(format, std::forward<Args>(args)...));
    }
    
    template<typename... Args>
    void warn(const std::string& format, Args&&... args) {
        warn(formatMessage(format, std::forward<Args>(args)...));
    }
    
    template<typename... Args>
    void error(const std::string& format, Args&&... args) {
        error(formatMessage(format, std::forward<Args>(args)...));
    }
    
    template<typename... Args>
    void critical(const std::string& format, Args&&... args) {
        critical(formatMessage(format, std::forward<Args>(args)...));
    }
    
    void flush() { std::cout.flush(); }
    
    bool shouldLog(Level level) const { return static_cast<int>(level) >= static_cast<int>(level_); }

private:
    void log(Level level, const std::string& message) {
        if (!shouldLog(level)) return;
        
        const char* levelStr = "";
        switch (level) {
            case Level::TRACE: levelStr = "TRACE"; break;
            case Level::DEBUG: levelStr = "DEBUG"; break;
            case Level::INFO: levelStr = "INFO"; break;
            case Level::WARN: levelStr = "WARN"; break;
            case Level::ERROR: levelStr = "ERROR"; break;
            case Level::CRITICAL: levelStr = "CRITICAL"; break;
        }
        
        std::cout << "[" << name_ << "] [" << levelStr << "] " << message << std::endl;
    }
    
    template<typename... Args>
    std::string formatMessage(const std::string& format, Args&&... args) {
        // 简化版格式化，实际实现中应该使用fmt库
        std::string result = format;
        
        // 简单的{}替换
        size_t pos = 0;
        int argIndex = 0;
        std::string argsArray[] = {toString(args)...};
        
        while ((pos = result.find("{}", pos)) != std::string::npos) {
            if (argIndex < sizeof...(args)) {
                result.replace(pos, 2, argsArray[argIndex]);
                argIndex++;
            } else {
                result.replace(pos, 2, "[MISSING_ARG]");
            }
            pos += 1;
        }
        
        return result;
    }
    
    template<typename T>
    std::string toString(const T& value) {
        if constexpr (std::is_same_v<T, std::string>) {
            return value;
        } else if constexpr (std::is_arithmetic_v<T>) {
            return std::to_string(value);
        } else {
            return "[UNKNOWN_TYPE]";
        }
    }
    
    std::string name_;
    Level level_;
};

// 模拟日志工厂
class SimpleLoggerFactory {
public:
    static SimpleLogger& getDefaultLogger() {
        static SimpleLogger logger("default");
        return logger;
    }
};

// 模拟全局日志宏
#define SIMPLE_TRACE(...)   SimpleLoggerFactory::getDefaultLogger().trace(__VA_ARGS__)
#define SIMPLE_DEBUG(...)   SimpleLoggerFactory::getDefaultLogger().debug(__VA_ARGS__)
#define SIMPLE_INFO(...)    SimpleLoggerFactory::getDefaultLogger().info(__VA_ARGS__)
#define SIMPLE_WARN(...)    SimpleLoggerFactory::getDefaultLogger().warn(__VA_ARGS__)
#define SIMPLE_ERROR(...)   SimpleLoggerFactory::getDefaultLogger().error(__VA_ARGS__)
#define SIMPLE_CRITICAL(...) SimpleLoggerFactory::getDefaultLogger().critical(__VA_ARGS__)

int main() {
    std::cout << "=== 简化版日志接口测试 ===" << std::endl;
    
    auto& logger = SimpleLoggerFactory::getDefaultLogger();
    
    // 测试不同级别的日志输出
    std::cout << "\n1. 测试不同级别的日志输出:" << std::endl;
    
    logger.trace("这是一条TRACE级别的日志消息");
    logger.debug("这是一条DEBUG级别的日志消息");
    logger.info("这是一条INFO级别的日志消息");
    logger.warn("这是一条WARN级别的日志消息");
    logger.error("这是一条ERROR级别的日志消息");
    logger.critical("这是一条CRITICAL级别的日志消息");
    
    // 测试格式化日志
    std::cout << "\n2. 测试格式化日志输出:" << std::endl;
    
    logger.trace("格式化TRACE日志: 参数1={}, 参数2={}", 123, "test");
    logger.debug("格式化DEBUG日志: 参数1={}, 参数2={}", 456, "debug");
    logger.info("格式化INFO日志: 参数1={}, 参数2={}", 789, "info");
    logger.warn("格式化WARN日志: 参数1={}, 参数2={}", 101, "warn");
    logger.error("格式化ERROR日志: 参数1={}, 参数2={}", 202, "error");
    logger.critical("格式化CRITICAL日志: 参数1={}, 参数2={}", 303, "critical");
    
    // 测试日志级别设置
    std::cout << "\n3. 测试日志级别设置:" << std::endl;
    
    auto currentLevel = logger.getLevel();
    std::cout << "当前日志级别: " << static_cast<int>(currentLevel) << std::endl;
    
    // 设置更严格的日志级别
    logger.setLevel(SimpleLogger::Level::WARN);
    std::cout << "设置日志级别为WARN后:" << std::endl;
    
    logger.trace("这条TRACE日志应该不会显示");
    logger.debug("这条DEBUG日志应该不会显示");
    logger.info("这条INFO日志应该不会显示");
    logger.warn("这条WARN日志应该显示");
    logger.error("这条ERROR日志应该显示");
    logger.critical("这条CRITICAL日志应该显示");
    
    // 恢复默认级别
    logger.setLevel(SimpleLogger::Level::INFO);
    
    // 测试全局日志宏
    std::cout << "\n4. 测试全局日志宏:" << std::endl;
    
    SIMPLE_TRACE("全局TRACE宏: 参数={}", "trace_macro");
    SIMPLE_DEBUG("全局DEBUG宏: 参数={}", "debug_macro");
    SIMPLE_INFO("全局INFO宏: 参数={}", "info_macro");
    SIMPLE_WARN("全局WARN宏: 参数={}", "warn_macro");
    SIMPLE_ERROR("全局ERROR宏: 参数={}", "error_macro");
    SIMPLE_CRITICAL("全局CRITICAL宏: 参数={}", "critical_macro");
    
    // 测试日志刷新
    logger.flush();
    
    std::cout << "\n=== 简化版日志接口测试完成 ===" << std::endl;
    
    return 0;
}