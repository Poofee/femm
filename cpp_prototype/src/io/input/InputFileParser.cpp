/**
 * @file InputFileParser.cpp
 * @brief Elmer输入文件解析器实现
 * 
 * 移植自Fortran版本的输入文件解析功能，支持.sif文件格式
 */

#include "InputFileParser.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <regex>

namespace elmer {

// ===== 构造函数 =====

InputFileParser::InputFileParser() 
    : isParsed_(false) {
}

// ===== 私有辅助函数 =====

std::string InputFileParser::trim(const std::string& str) const {
    size_t start = str.find_first_not_of(" \t\n\r\f\v");
    if (start == std::string::npos) {
        return "";
    }
    size_t end = str.find_last_not_of(" \t\n\r\f\v");
    return str.substr(start, end - start + 1);
}

std::string InputFileParser::toLower(const std::string& str) const {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), 
                   [](unsigned char c) { return std::tolower(c); });
    return result;
}

bool InputFileParser::isSectionStart(const std::string& line, std::string& sectionName, 
                                    std::string& sectionType, int& sectionId) const {
    std::string trimmed = trim(line);
    
    // 检查是否是节段开始（如 "Header", "Simulation", "Body 1", "Solver 2"）
    if (trimmed.empty() || trimmed[0] == '!') {
        return false; // 空行或注释
    }
    
    // 使用正则表达式匹配节段格式
    std::regex sectionRegex(R"((\w+)(?:\s+(\d+))?)");
    std::smatch matches;
    
    if (std::regex_match(trimmed, matches, sectionRegex)) {
        sectionName = matches[1];
        sectionType = toLower(sectionName);
        
        if (matches.size() > 2 && matches[2].matched) {
            sectionId = std::stoi(matches[2]);
        } else {
            sectionId = 0; // 无ID的节段（如Header, Simulation）
        }
        
        return true;
    }
    
    return false;
}

bool InputFileParser::isParameterLine(const std::string& line, std::string& paramName, 
                                     std::string& paramValue) const {
    std::string trimmed = trim(line);
    
    if (trimmed.empty() || trimmed[0] == '!' || trimmed == "End") {
        return false; // 空行、注释或节段结束
    }
    
    // 查找等号位置
    size_t equalPos = trimmed.find('=');
    if (equalPos == std::string::npos) {
        return false; // 没有等号，不是参数行
    }
    
    paramName = trim(trimmed.substr(0, equalPos));
    paramValue = trim(trimmed.substr(equalPos + 1));
    
    return !paramName.empty() && !paramValue.empty();
}

bool InputFileParser::isEndSection(const std::string& line, const std::string& sectionName) const {
    std::string trimmed = trim(line);
    return trimmed == "End";
}

void InputFileParser::parseParameterValue(const std::string& valueStr, InputParameter& param) const {
    param.value = valueStr;
    
    // 尝试判断参数类型
    if (valueStr.empty()) {
        param.type = "String";
    } else if (std::regex_match(valueStr, std::regex(R"(\s*[+-]?\d+\.?\d*\s*)"))) {
        param.type = "Real";
    } else if (std::regex_match(valueStr, std::regex(R"(\s*[+-]?\d+\s*)"))) {
        param.type = "Integer";
    } else if (std::regex_match(valueStr, std::regex(R"(\s*(true|false|True|False|TRUE|FALSE)\s*)"))) {
        param.type = "Logical";
    } else if (valueStr.front() == '"' && valueStr.back() == '"') {
        param.type = "String";
        param.value = valueStr.substr(1, valueStr.length() - 2); // 去除引号
    } else {
        param.type = "String";
    }
}

// ===== 公共接口函数 =====

bool InputFileParser::parse(const std::string& filename) {
    filename_ = filename;
    sections_.clear();
    sectionMap_.clear();
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误: 无法打开输入文件: " << filename << std::endl;
        return false;
    }
    
    std::string line;
    int lineNumber = 0;
    InputSection* currentSection = nullptr;
    
    while (std::getline(file, line)) {
        lineNumber++;
        std::string trimmedLine = trim(line);
        
        // 跳过空行和注释
        if (trimmedLine.empty() || trimmedLine[0] == '!') {
            continue;
        }
        
        // 检查是否是节段开始
        std::string sectionName, sectionType;
        int sectionId;
        if (isSectionStart(trimmedLine, sectionName, sectionType, sectionId)) {
            // 创建新节段
            InputSection newSection(sectionName, sectionType, sectionId);
            sections_.push_back(newSection);
            currentSection = &sections_.back();
            sectionMap_[sectionType].push_back(newSection);
            continue;
        }
        
        // 检查是否是节段结束
        if (currentSection && isEndSection(trimmedLine, currentSection->name)) {
            currentSection = nullptr;
            continue;
        }
        
        // 检查是否是参数行
        if (currentSection) {
            std::string paramName, paramValue;
            if (isParameterLine(trimmedLine, paramName, paramValue)) {
                InputParameter param;
                param.name = paramName;
                param.section = currentSection->name;
                param.lineNumber = lineNumber;
                parseParameterValue(paramValue, param);
                
                currentSection->parameters[paramName] = param;
            }
        }
    }
    
    file.close();
    isParsed_ = true;
    
    std::cout << "输入文件解析完成: " << filename << std::endl;
    std::cout << "解析节段数量: " << sections_.size() << std::endl;
    
    return true;
}

std::vector<InputSection> InputFileParser::getSectionsByType(const std::string& type) const {
    std::string lowerType = toLower(type);
    auto it = sectionMap_.find(lowerType);
    if (it != sectionMap_.end()) {
        return it->second;
    }
    return std::vector<InputSection>();
}

InputSection InputFileParser::getSection(const std::string& type, int id) const {
    std::string lowerType = toLower(type);
    auto it = sectionMap_.find(lowerType);
    if (it != sectionMap_.end()) {
        for (const auto& section : it->second) {
            if (section.id == id) {
                return section;
            }
        }
    }
    
    // 返回空节段
    return InputSection();
}

bool InputFileParser::hasSection(const std::string& type, int id) const {
    std::string lowerType = toLower(type);
    auto it = sectionMap_.find(lowerType);
    if (it != sectionMap_.end()) {
        for (const auto& section : it->second) {
            if (section.id == id) {
                return true;
            }
        }
    }
    return false;
}

std::string InputFileParser::getParameterValue(const std::string& sectionType, int sectionId, 
                                              const std::string& paramName) const {
    InputSection section = getSection(sectionType, sectionId);
    auto it = section.parameters.find(paramName);
    if (it != section.parameters.end()) {
        return it->second.value;
    }
    return "";
}

double InputFileParser::getParameterReal(const std::string& sectionType, int sectionId, 
                                        const std::string& paramName, double defaultValue) const {
    std::string value = getParameterValue(sectionType, sectionId, paramName);
    if (value.empty()) {
        return defaultValue;
    }
    
    try {
        return std::stod(value);
    } catch (const std::exception& e) {
        std::cerr << "警告: 无法将参数 '" << paramName << "' 转换为实数: " << value << std::endl;
        return defaultValue;
    }
}

int InputFileParser::getParameterInteger(const std::string& sectionType, int sectionId, 
                                        const std::string& paramName, int defaultValue) const {
    std::string value = getParameterValue(sectionType, sectionId, paramName);
    if (value.empty()) {
        return defaultValue;
    }
    
    try {
        return std::stoi(value);
    } catch (const std::exception& e) {
        std::cerr << "警告: 无法将参数 '" << paramName << "' 转换为整数: " << value << std::endl;
        return defaultValue;
    }
}

bool InputFileParser::getParameterLogical(const std::string& sectionType, int sectionId, 
                                         const std::string& paramName, bool defaultValue) const {
    std::string value = getParameterValue(sectionType, sectionId, paramName);
    if (value.empty()) {
        return defaultValue;
    }
    
    std::string lowerValue = toLower(value);
    if (lowerValue == "true" || lowerValue == "1" || lowerValue == "yes") {
        return true;
    } else if (lowerValue == "false" || lowerValue == "0" || lowerValue == "no") {
        return false;
    } else {
        std::cerr << "警告: 无法将参数 '" << paramName << "' 转换为逻辑值: " << value << std::endl;
        return defaultValue;
    }
}

std::string InputFileParser::getSimulationType() const {
    return getParameterValue("simulation", 0, "Simulation Type");
}

std::string InputFileParser::getMeshFileName() const {
    // 从Header节段获取网格文件名
    InputSection header = getSection("header", 0);
    if (!header.name.empty()) {
        auto it = header.parameters.find("Mesh DB");
        if (it != header.parameters.end()) {
            // Mesh DB参数格式: "目录" "文件名"
            // 使用简单的字符串处理代替复杂的正则表达式
            std::string value = it->second.value;
            
            // 查找第一个引号
            size_t firstQuote = value.find('"');
            if (firstQuote != std::string::npos) {
                // 查找第二个引号
                size_t secondQuote = value.find('"', firstQuote + 1);
                if (secondQuote != std::string::npos) {
                    // 查找第三个引号
                    size_t thirdQuote = value.find('"', secondQuote + 1);
                    if (thirdQuote != std::string::npos) {
                        // 查找第四个引号
                        size_t fourthQuote = value.find('"', thirdQuote + 1);
                        if (fourthQuote != std::string::npos) {
                            // 提取文件名（第三个和第四个引号之间的内容）
                            return value.substr(thirdQuote + 1, fourthQuote - thirdQuote - 1);
                        }
                    }
                }
            }
        }
    }
    return "";
}

std::string InputFileParser::getMeshDirectory() const {
    // 从Header节段获取网格目录
    InputSection header = getSection("header", 0);
    if (!header.name.empty()) {
        auto it = header.parameters.find("Mesh DB");
        if (it != header.parameters.end()) {
            std::string value = it->second.value;
            
            // 使用简单的字符串处理代替复杂的正则表达式
            // 查找第一个引号
            size_t firstQuote = value.find('"');
            if (firstQuote != std::string::npos) {
                // 查找第二个引号
                size_t secondQuote = value.find('"', firstQuote + 1);
                if (secondQuote != std::string::npos) {
                    // 提取目录（第一个和第二个引号之间的内容）
                    return value.substr(firstQuote + 1, secondQuote - firstQuote - 1);
                }
            }
        }
    }
    return "."; // 默认当前目录
}

std::vector<int> InputFileParser::getActiveSolvers() const {
    std::vector<int> activeSolvers;
    
    // 查找所有Equation节段中的Active Solvers参数
    auto equationSections = getSectionsByType("equation");
    for (const auto& equation : equationSections) {
        std::string activeSolversStr = getParameterValue("equation", equation.id, "Active Solvers");
        
        if (!activeSolversStr.empty()) {
            // 解析格式如: "Active Solvers(1)=1" 或 "Active Solvers(2)=2 3"
            std::regex solverRegex(R"(\d+)");
            std::sregex_iterator it(activeSolversStr.begin(), activeSolversStr.end(), solverRegex);
            std::sregex_iterator end;
            
            while (it != end) {
                int solverId = std::stoi(it->str());
                activeSolvers.push_back(solverId);
                ++it;
            }
        }
    }
    
    return activeSolvers;
}

bool InputFileParser::validate() const {
    if (!isParsed_) {
        std::cerr << "错误: 输入文件未解析" << std::endl;
        return false;
    }
    
    // 检查必需节段
    if (!hasSection("simulation", 0)) {
        std::cerr << "错误: 缺少Simulation节段" << std::endl;
        return false;
    }
    
    // 检查必需参数
    if (getSimulationType().empty()) {
        std::cerr << "错误: 缺少Simulation Type参数" << std::endl;
        return false;
    }
    
    // 检查网格文件
    if (getMeshFileName().empty()) {
        std::cerr << "错误: 缺少网格文件名" << std::endl;
        return false;
    }
    
    return true;
}

void InputFileParser::printSummary() const {
    if (!isParsed_) {
        std::cout << "输入文件未解析" << std::endl;
        return;
    }
    
    std::cout << "=== 输入文件解析摘要 ===" << std::endl;
    std::cout << "文件: " << filename_ << std::endl;
    std::cout << "节段数量: " << sections_.size() << std::endl;
    
    for (const auto& section : sections_) {
        std::cout << "\n节段: " << section.name;
        if (section.id > 0) {
            std::cout << " " << section.id;
        }
        std::cout << " (类型: " << section.type << ")" << std::endl;
        
        for (const auto& paramPair : section.parameters) {
            const InputParameter& param = paramPair.second;
            std::cout << "  " << param.name << " = " << param.value 
                      << " (类型: " << param.type << ")" << std::endl;
        }
    }
    
    std::cout << "\n=== 仿真信息 ===" << std::endl;
    std::cout << "仿真类型: " << getSimulationType() << std::endl;
    std::cout << "网格目录: " << getMeshDirectory() << std::endl;
    std::cout << "网格文件: " << getMeshFileName() << std::endl;
    
    auto activeSolvers = getActiveSolvers();
    std::cout << "激活求解器: ";
    for (int solverId : activeSolvers) {
        std::cout << solverId << " ";
    }
    std::cout << std::endl;
}

} // namespace elmer