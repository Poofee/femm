/**
 * @file test_logger.cpp
 * @brief 测试日志接口功能
 * 
 * 验证日志接口和spdlog适配器的功能是否正常。
 */

#include "LoggerFactory.h"
#include <iostream>

int main() {
    std::cout << "开始测试日志接口..." << std::endl;
    
    // 测试默认日志器
    auto logger = elmer::LoggerFactory::getDefaultLogger();
    
    if (!logger) {
        std::cerr << "错误: 无法获取默认日志器" << std::endl;
        return 1;
    }
    
    std::cout << "默认日志器获取成功" << std::endl;
    
    // 测试不同级别的日志输出
    std::cout << "\n测试不同级别的日志输出:" << std::endl;
    
    logger->trace("这是一条TRACE级别的日志消息");
    logger->debug("这是一条DEBUG级别的日志消息");
    logger->info("这是一条INFO级别的日志消息");
    logger->warn("这是一条WARN级别的日志消息");
    logger->error("这是一条ERROR级别的日志消息");
    logger->critical("这是一条CRITICAL级别的日志消息");
    
    // 测试格式化日志
    std::cout << "\n测试格式化日志输出:" << std::endl;
    
    logger->trace("格式化TRACE日志: 参数1={}, 参数2={}", 123, "test");
    logger->debug("格式化DEBUG日志: 参数1={}, 参数2={}", 456, "debug");
    logger->info("格式化INFO日志: 参数1={}, 参数2={}", 789, "info");
    logger->warn("格式化WARN日志: 参数1={}, 参数2={}", 101, "warn");
    logger->error("格式化ERROR日志: 参数1={}, 参数2={}", 202, "error");
    logger->critical("格式化CRITICAL日志: 参数1={}, 参数2={}", 303, "critical");
    
    // 测试日志级别设置
    std::cout << "\n测试日志级别设置:" << std::endl;
    
    auto currentLevel = logger->getLevel();
    std::cout << "当前日志级别: " << static_cast<int>(currentLevel) << std::endl;
    
    // 设置更严格的日志级别
    logger->setLevel(elmer::LogLevel::WARN);
    std::cout << "设置日志级别为WARN后:" << std::endl;
    
    logger->trace("这条TRACE日志应该不会显示");
    logger->debug("这条DEBUG日志应该不会显示");
    logger->info("这条INFO日志应该不会显示");
    logger->warn("这条WARN日志应该显示");
    logger->error("这条ERROR日志应该显示");
    logger->critical("这条CRITICAL日志应该显示");
    
    // 恢复默认级别
    logger->setLevel(elmer::LogLevel::INFO);
    
    // 测试创建多个日志器
    std::cout << "\n测试创建多个日志器:" << std::endl;
    
    auto customLogger = elmer::LoggerFactory::createLogger("custom");
    if (customLogger) {
        customLogger->info("自定义日志器测试消息");
    }
    
    // 测试全局日志宏
    std::cout << "\n测试全局日志宏:" << std::endl;
    
    ELMER_TRACE("全局TRACE宏: 参数={}", "trace_macro");
    ELMER_DEBUG("全局DEBUG宏: 参数={}", "debug_macro");
    ELMER_INFO("全局INFO宏: 参数={}", "info_macro");
    ELMER_WARN("全局WARN宏: 参数={}", "warn_macro");
    ELMER_ERROR("全局ERROR宏: 参数={}", "error_macro");
    ELMER_CRITICAL("全局CRITICAL宏: 参数={}", "critical_macro");
    
    // 测试日志刷新
    logger->flush();
    
    std::cout << "\n日志接口测试完成" << std::endl;
    
    return 0;
}