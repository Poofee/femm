/**
 * @file Lists.cpp
 * @brief 列表处理工具模块实现
 * 
 * 对应Fortran模块：Lists.F90
 * 提供Elmer FEM中的列表处理功能实现。
 */

#include "Lists.h"
#include "ElmerCpp.h"
#include <algorithm>
#include <iostream>
#include <sstream>

namespace elmer {

// 静态成员变量定义
std::map<std::string, std::shared_ptr<List>> Lists::lists;
std::string Lists::currentNamespace;

// List类方法实现
std::shared_ptr<ListEntry> List::Add(const std::string& entryName) {
    // 检查是否已存在同名条目
    if (Contains(entryName)) {
        std::cout << "警告: 列表项 '" << entryName << "' 已存在" << std::endl;
        return Get(entryName);
    }
    
    // 创建新条目
    auto newEntry = std::make_shared<ListEntry>(entryName, ListType::STRING);
    
    // 添加到链表头部
    newEntry->next = head;
    head = newEntry;
    
    return newEntry;
}

void List::Remove(const std::string& entryName) {
    if (head == nullptr) {
        return;
    }
    
    // 如果头节点就是要删除的节点
    if (head->name == entryName) {
        head = head->next;
        return;
    }
    
    // 遍历链表查找要删除的节点
    auto current = head;
    while (current->next != nullptr) {
        if (current->next->name == entryName) {
            current->next = current->next->next;
            return;
        }
        current = current->next;
    }
}

std::shared_ptr<ListEntry> List::Get(const std::string& entryName) const {
    auto current = head;
    while (current != nullptr) {
        if (current->name == entryName) {
            return current;
        }
        current = current->next;
    }
    return nullptr;
}

bool List::Contains(const std::string& entryName) const {
    return Get(entryName) != nullptr;
}

void List::Clear() {
    head = nullptr;
}

size_t List::Size() const {
    size_t count = 0;
    auto current = head;
    while (current != nullptr) {
        count++;
        current = current->next;
    }
    return count;
}

// Lists类静态方法实现
std::shared_ptr<List> Lists::CreateList(const std::string& listName) {
    auto list = std::make_shared<List>(listName);
    lists[listName] = list;
    return list;
}

std::shared_ptr<List> Lists::GetList(const std::string& listName) {
    auto it = lists.find(listName);
    if (it != lists.end()) {
        return it->second;
    }
    return nullptr;
}

void Lists::DeleteList(const std::string& listName) {
    lists.erase(listName);
}

void Lists::SetNamespace(const std::string& namespaceStr) {
    currentNamespace = namespaceStr;
}

std::string Lists::GetNamespace() {
    return currentNamespace;
}

bool Lists::AddLogical(std::shared_ptr<List> list, const std::string& name, bool value) {
    if (!list) {
        std::cout << "错误: 列表为空" << std::endl;
        return false;
    }
    
    auto entry = list->Add(name);
    entry->type = ListType::LOGICAL;
    entry->value = value;
    return true;
}

bool Lists::AddInteger(std::shared_ptr<List> list, const std::string& name, int value) {
    if (!list) {
        std::cout << "错误: 列表为空" << std::endl;
        return false;
    }
    
    auto entry = list->Add(name);
    entry->type = ListType::INTEGER;
    entry->value = value;
    return true;
}

bool Lists::AddReal(std::shared_ptr<List> list, const std::string& name, double value) {
    if (!list) {
        std::cout << "错误: 列表为空" << std::endl;
        return false;
    }
    
    auto entry = list->Add(name);
    entry->type = ListType::CONSTANT_SCALAR;
    entry->value = value;
    return true;
}

bool Lists::AddString(std::shared_ptr<List> list, const std::string& name, const std::string& value) {
    if (!list) {
        std::cout << "错误: 列表为空" << std::endl;
        return false;
    }
    
    auto entry = list->Add(name);
    entry->type = ListType::STRING;
    entry->value = value;
    return true;
}

std::optional<bool> Lists::GetLogical(std::shared_ptr<List> list, const std::string& name) {
    if (!list) {
        return std::nullopt;
    }
    
    auto entry = list->Get(name);
    if (!entry || entry->type != ListType::LOGICAL) {
        return std::nullopt;
    }
    
    try {
        return std::get<bool>(entry->value);
    } catch (const std::bad_variant_access&) {
        return std::nullopt;
    }
}

std::optional<int> Lists::GetInteger(std::shared_ptr<List> list, const std::string& name) {
    if (!list) {
        return std::nullopt;
    }
    
    auto entry = list->Get(name);
    if (!entry || entry->type != ListType::INTEGER) {
        return std::nullopt;
    }
    
    try {
        return std::get<int>(entry->value);
    } catch (const std::bad_variant_access&) {
        return std::nullopt;
    }
}

std::optional<double> Lists::GetReal(std::shared_ptr<List> list, const std::string& name) {
    if (!list) {
        return std::nullopt;
    }
    
    auto entry = list->Get(name);
    if (!entry || (entry->type != ListType::CONSTANT_SCALAR && 
                   entry->type != ListType::VARIABLE_SCALAR)) {
        return std::nullopt;
    }
    
    try {
        return std::get<double>(entry->value);
    } catch (const std::bad_variant_access&) {
        return std::nullopt;
    }
}

std::optional<std::string> Lists::GetString(std::shared_ptr<List> list, const std::string& name) {
    if (!list) {
        return std::nullopt;
    }
    
    auto entry = list->Get(name);
    if (!entry || entry->type != ListType::STRING) {
        return std::nullopt;
    }
    
    try {
        return std::get<std::string>(entry->value);
    } catch (const std::bad_variant_access&) {
        return std::nullopt;
    }
}

bool Lists::Check(std::shared_ptr<List> list, const std::string& name) {
    return list && list->Contains(name);
}

std::vector<std::string> Lists::GetEntryNames(std::shared_ptr<List> list) {
    std::vector<std::string> names;
    if (!list) {
        return names;
    }
    
    auto current = list->GetHead();
    while (current != nullptr) {
        names.push_back(current->name);
        current = current->next;
    }
    
    return names;
}

void Lists::PrintList(std::shared_ptr<List> list) {
    if (!list) {
        std::cout << "列表为空" << std::endl;
        return;
    }
    
    std::cout << "列表 '" << list->GetName() << "' 内容:" << std::endl;
    
    auto current = list->GetHead();
    while (current != nullptr) {
        std::cout << "  " << current->name << " : ";
        
        switch (current->type) {
            case ListType::LOGICAL:
                std::cout << std::boolalpha << std::get<bool>(current->value);
                break;
            case ListType::INTEGER:
                std::cout << std::get<int>(current->value);
                break;
            case ListType::CONSTANT_SCALAR:
            case ListType::VARIABLE_SCALAR:
                std::cout << std::get<double>(current->value);
                break;
            case ListType::STRING:
                std::cout << "\"" << std::get<std::string>(current->value) << "\"";
                break;
            default:
                std::cout << "未知类型";
                break;
        }
        
        std::cout << std::endl;
        current = current->next;
    }
}

} // namespace elmer