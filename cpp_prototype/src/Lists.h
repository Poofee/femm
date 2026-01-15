/**
 * @file Lists.h
 * @brief 列表处理工具模块
 * 
 * 对应Fortran模块：Lists.F90
 * 提供Elmer FEM中的列表处理功能，包括关键字列表的创建、访问和管理。
 */

#pragma once

#include "Types.h"
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <variant>
#include <optional>
#include <stdexcept>
#include <iostream>

namespace elmer {

/**
 * @brief 列表数据类型枚举
 */
enum class ListType {
    LOGICAL = 1,                    // 逻辑类型
    STRING = 2,                     // 字符串类型
    INTEGER = 3,                    // 整数类型
    CONSTANT_SCALAR = 4,            // 常量标量
    VARIABLE_SCALAR = 5,            // 变量标量
    CONSTANT_SCALAR_STR = 6,        // 常量标量字符串
    VARIABLE_SCALAR_STR = 7,        // 变量标量字符串
    CONSTANT_SCALAR_PROC = 8,       // 常量标量过程
    CONSTANT_TENSOR = 9,            // 常量张量
    VARIABLE_TENSOR = 10,           // 变量张量
    CONSTANT_TENSOR_STR = 11,       // 常量张量字符串
    VARIABLE_TENSOR_STR = 12        // 变量张量字符串
};

/**
 * @brief 列表项数据结构
 */
struct ListEntry {
    std::string name;               // 条目名称
    ListType type;                  // 数据类型
    std::variant<bool, int, double, std::string> value; // 存储的值
    std::shared_ptr<ListEntry> next; // 下一个条目指针
    
    // 构造函数
    ListEntry(const std::string& entryName, ListType entryType) 
        : name(entryName), type(entryType), next(nullptr) {}
};

/**
 * @brief 列表数据结构
 */
class List {
private:
    std::shared_ptr<ListEntry> head; // 列表头指针
    std::string name;                // 列表名称
    
public:
    /**
     * @brief 构造函数
     */
    explicit List(const std::string& listName = "") : head(nullptr), name(listName) {}
    
    /**
     * @brief 析构函数
     */
    ~List() = default;
    
    /**
     * @brief 获取列表名称
     */
    const std::string& GetName() const { return name; }
    
    /**
     * @brief 设置列表名称
     */
    void SetName(const std::string& listName) { name = listName; }
    
    /**
     * @brief 获取列表头指针
     */
    std::shared_ptr<ListEntry> GetHead() const { return head; }
    
    /**
     * @brief 检查列表是否为空
     */
    bool IsEmpty() const { return head == nullptr; }
    
    /**
     * @brief 添加列表项
     */
    std::shared_ptr<ListEntry> Add(const std::string& entryName);
    
    /**
     * @brief 删除列表项
     */
    void Remove(const std::string& entryName);
    
    /**
     * @brief 获取列表项
     */
    std::shared_ptr<ListEntry> Get(const std::string& entryName) const;
    
    /**
     * @brief 检查列表项是否存在
     */
    bool Contains(const std::string& entryName) const;
    
    /**
     * @brief 清空列表
     */
    void Clear();
    
    /**
     * @brief 获取列表大小
     */
    size_t Size() const;
};

/**
 * @brief 列表工具类
 */
class Lists {
private:
    static std::map<std::string, std::shared_ptr<List>> lists; // 全局列表存储
    static std::string currentNamespace;                       // 当前命名空间
    
public:
    /**
     * @brief 创建新列表
     */
    static std::shared_ptr<List> CreateList(const std::string& listName);
    
    /**
     * @brief 获取列表
     */
    static std::shared_ptr<List> GetList(const std::string& listName);
    
    /**
     * @brief 删除列表
     */
    static void DeleteList(const std::string& listName);
    
    /**
     * @brief 设置当前命名空间
     */
    static void SetNamespace(const std::string& namespaceStr);
    
    /**
     * @brief 获取当前命名空间
     */
    static std::string GetNamespace();
    
    /**
     * @brief 添加逻辑值到列表
     */
    static bool AddLogical(std::shared_ptr<List> list, const std::string& name, bool value);
    
    /**
     * @brief 添加整数值到列表
     */
    static bool AddInteger(std::shared_ptr<List> list, const std::string& name, int value);
    
    /**
     * @brief 添加实数值到列表
     */
    static bool AddReal(std::shared_ptr<List> list, const std::string& name, double value);
    
    /**
     * @brief 添加字符串到列表
     */
    static bool AddString(std::shared_ptr<List> list, const std::string& name, const std::string& value);
    
    /**
     * @brief 从列表获取逻辑值
     */
    static std::optional<bool> GetLogical(std::shared_ptr<List> list, const std::string& name);
    
    /**
     * @brief 从列表获取整数值
     */
    static std::optional<int> GetInteger(std::shared_ptr<List> list, const std::string& name);
    
    /**
     * @brief 从列表获取实数值
     */
    static std::optional<double> GetReal(std::shared_ptr<List> list, const std::string& name);
    
    /**
     * @brief 从列表获取字符串值
     */
    static std::optional<std::string> GetString(std::shared_ptr<List> list, const std::string& name);
    
    /**
     * @brief 检查列表项是否存在
     */
    static bool Check(std::shared_ptr<List> list, const std::string& name);
    
    /**
     * @brief 获取列表中的所有条目名称
     */
    static std::vector<std::string> GetEntryNames(std::shared_ptr<List> list);
    
    /**
     * @brief 打印列表内容（用于调试）
     */
    static void PrintList(std::shared_ptr<List> list);
};

} // namespace elmer