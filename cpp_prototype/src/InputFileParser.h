/**
 * @file InputFileParser.h
 * @brief Elmer输入文件解析器
 * 
 * 移植自Fortran版本的输入文件解析功能，支持.sif文件格式
 */

#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

namespace elmer {

/**
 * @brief 输入文件参数结构
 */
struct InputParameter {
    std::string name;           ///< 参数名称
    std::string value;          ///< 参数值（字符串形式）
    std::string type;           ///< 参数类型（Real, Integer, Logical, String）
    std::string section;        ///< 所属节段
    int lineNumber;             ///< 行号
    
    InputParameter() : lineNumber(0) {}
    InputParameter(const std::string& n, const std::string& v, const std::string& t, 
                   const std::string& s, int ln) 
        : name(n), value(v), type(t), section(s), lineNumber(ln) {}
};

/**
 * @brief 输入文件节段结构
 */
struct InputSection {
    std::string name;           ///< 节段名称
    std::string type;           ///< 节段类型（Header, Simulation, Body, Solver等）
    int id;                     ///< 节段ID
    std::map<std::string, InputParameter> parameters;  ///< 参数映射
    
    InputSection() : id(0) {}
    InputSection(const std::string& n, const std::string& t, int i) 
        : name(n), type(t), id(i) {}
};

/**
 * @brief 输入文件解析器类
 */
class InputFileParser {
private:
    std::string filename_;                      ///< 输入文件名
    std::vector<InputSection> sections_;        ///< 解析的节段列表
    std::map<std::string, std::vector<InputSection>> sectionMap_; ///< 按类型分组的节段
    bool isParsed_;                             ///< 是否已解析
    
    // 解析辅助函数
    std::string trim(const std::string& str) const;
    std::string toLower(const std::string& str) const;
    bool isSectionStart(const std::string& line, std::string& sectionName, std::string& sectionType, int& sectionId) const;
    bool isParameterLine(const std::string& line, std::string& paramName, std::string& paramValue) const;
    bool isEndSection(const std::string& line, const std::string& sectionName) const;
    void parseParameterValue(const std::string& valueStr, InputParameter& param) const;
    
public:
    InputFileParser();
    ~InputFileParser() = default;
    
    /**
     * @brief 解析输入文件
     * @param filename 输入文件名
     * @return 是否解析成功
     */
    bool parse(const std::string& filename);
    
    /**
     * @brief 获取所有节段
     */
    const std::vector<InputSection>& getSections() const { return sections_; }
    
    /**
     * @brief 根据类型获取节段
     */
    std::vector<InputSection> getSectionsByType(const std::string& type) const;
    
    /**
     * @brief 获取特定节段
     */
    InputSection getSection(const std::string& type, int id) const;
    
    /**
     * @brief 检查是否包含特定节段
     */
    bool hasSection(const std::string& type, int id) const;
    
    /**
     * @brief 获取参数值（字符串形式）
     */
    std::string getParameterValue(const std::string& sectionType, int sectionId, 
                                 const std::string& paramName) const;
    
    /**
     * @brief 获取参数值（实数形式）
     */
    double getParameterReal(const std::string& sectionType, int sectionId, 
                           const std::string& paramName, double defaultValue = 0.0) const;
    
    /**
     * @brief 获取参数值（整数形式）
     */
    int getParameterInteger(const std::string& sectionType, int sectionId, 
                           const std::string& paramName, int defaultValue = 0) const;
    
    /**
     * @brief 获取参数值（逻辑值形式）
     */
    bool getParameterLogical(const std::string& sectionType, int sectionId, 
                            const std::string& paramName, bool defaultValue = false) const;
    
    /**
     * @brief 获取仿真类型
     */
    std::string getSimulationType() const;
    
    /**
     * @brief 获取网格文件名
     */
    std::string getMeshFileName() const;
    
    /**
     * @brief 获取网格目录
     */
    std::string getMeshDirectory() const;
    
    /**
     * @brief 获取激活的求解器列表
     */
    std::vector<int> getActiveSolvers() const;
    
    /**
     * @brief 验证输入文件的完整性
     */
    bool validate() const;
    
    /**
     * @brief 打印解析结果（调试用）
     */
    void printSummary() const;
};

} // namespace elmer