/**
 * @file test_lists.cpp
 * @brief Lists模块测试文件
 * 
 * 测试Lists模块的核心功能，包括列表创建、添加、获取、删除等操作。
 */

#include "../src/Lists.h"
#include <iostream>
#include <cassert>

using namespace elmer;

/**
 * @brief 测试列表基本操作
 */
void testListBasicOperations() {
    std::cout << "=== 测试列表基本操作 ===" << std::endl;
    
    // 创建列表
    auto list = Lists::CreateList("TestList");
    assert(list != nullptr);
    assert(list->GetName() == "TestList");
    assert(list->IsEmpty());
    
    // 添加条目
    auto entry1 = list->Add("Item1");
    assert(entry1 != nullptr);
    assert(entry1->name == "Item1");
    assert(!list->IsEmpty());
    assert(list->Size() == 1);
    
    auto entry2 = list->Add("Item2");
    assert(entry2 != nullptr);
    assert(list->Size() == 2);
    
    // 检查条目存在性
    assert(list->Contains("Item1"));
    assert(list->Contains("Item2"));
    assert(!list->Contains("Item3"));
    
    // 获取条目
    auto retrievedEntry = list->Get("Item1");
    assert(retrievedEntry != nullptr);
    assert(retrievedEntry->name == "Item1");
    
    // 删除条目
    list->Remove("Item1");
    assert(!list->Contains("Item1"));
    assert(list->Size() == 1);
    
    // 清空列表
    list->Clear();
    assert(list->IsEmpty());
    assert(list->Size() == 0);
    
    std::cout << "PASSED: 列表基本操作测试" << std::endl;
}

/**
 * @brief 测试列表数据操作
 */
void testListDataOperations() {
    std::cout << "=== 测试列表数据操作 ===" << std::endl;
    
    auto list = Lists::CreateList("DataTestList");
    
    // 添加逻辑值
    assert(Lists::AddLogical(list, "BoolTrue", true));
    assert(Lists::AddLogical(list, "BoolFalse", false));
    
    // 添加整数值
    assert(Lists::AddInteger(list, "IntValue", 42));
    assert(Lists::AddInteger(list, "NegativeInt", -10));
    
    // 添加实数值
    assert(Lists::AddReal(list, "RealValue", 3.14159));
    assert(Lists::AddReal(list, "NegativeReal", -2.718));
    
    // 添加字符串
    assert(Lists::AddString(list, "StringValue", "Hello World"));
    assert(Lists::AddString(list, "EmptyString", ""));
    
    // 获取逻辑值
    auto boolTrue = Lists::GetLogical(list, "BoolTrue");
    assert(boolTrue.has_value() && boolTrue.value() == true);
    
    auto boolFalse = Lists::GetLogical(list, "BoolFalse");
    assert(boolFalse.has_value() && boolFalse.value() == false);
    
    // 获取整数值
    auto intValue = Lists::GetInteger(list, "IntValue");
    assert(intValue.has_value() && intValue.value() == 42);
    
    auto negativeInt = Lists::GetInteger(list, "NegativeInt");
    assert(negativeInt.has_value() && negativeInt.value() == -10);
    
    // 获取实数值
    auto realValue = Lists::GetReal(list, "RealValue");
    assert(realValue.has_value() && std::abs(realValue.value() - 3.14159) < 1e-10);
    
    auto negativeReal = Lists::GetReal(list, "NegativeReal");
    assert(negativeReal.has_value() && std::abs(negativeReal.value() - (-2.718)) < 1e-10);
    
    // 获取字符串
    auto stringValue = Lists::GetString(list, "StringValue");
    assert(stringValue.has_value() && stringValue.value() == "Hello World");
    
    auto emptyString = Lists::GetString(list, "EmptyString");
    assert(emptyString.has_value() && emptyString.value().empty());
    
    // 测试不存在的条目
    auto nonExistent = Lists::GetLogical(list, "NonExistent");
    assert(!nonExistent.has_value());
    
    std::cout << "PASSED: 列表数据操作测试" << std::endl;
}

/**
 * @brief 测试列表管理功能
 */
void testListManagement() {
    std::cout << "=== 测试列表管理功能 ===" << std::endl;
    
    // 创建多个列表
    auto list1 = Lists::CreateList("List1");
    auto list2 = Lists::CreateList("List2");
    auto list3 = Lists::CreateList("List3");
    
    // 验证列表创建
    assert(Lists::GetList("List1") == list1);
    assert(Lists::GetList("List2") == list2);
    assert(Lists::GetList("List3") == list3);
    
    // 添加数据到不同列表
    Lists::AddInteger(list1, "Value1", 100);
    Lists::AddInteger(list2, "Value2", 200);
    Lists::AddInteger(list3, "Value3", 300);
    
    // 验证数据隔离
    auto value1 = Lists::GetInteger(list1, "Value1");
    assert(value1.has_value() && value1.value() == 100);
    
    auto value2 = Lists::GetInteger(list2, "Value2");
    assert(value2.has_value() && value2.value() == 200);
    
    auto value3 = Lists::GetInteger(list3, "Value3");
    assert(value3.has_value() && value3.value() == 300);
    
    // 删除列表
    Lists::DeleteList("List2");
    assert(Lists::GetList("List2") == nullptr);
    
    // 验证其他列表不受影响
    assert(Lists::GetList("List1") == list1);
    assert(Lists::GetList("List3") == list3);
    
    std::cout << "PASSED: 列表管理功能测试" << std::endl;
}

/**
 * @brief 测试命名空间功能
 */
void testNamespaceFunctionality() {
    std::cout << "=== 测试命名空间功能 ===" << std::endl;
    
    // 设置命名空间
    Lists::SetNamespace("TestNamespace");
    assert(Lists::GetNamespace() == "TestNamespace");
    
    // 更改命名空间
    Lists::SetNamespace("AnotherNamespace");
    assert(Lists::GetNamespace() == "AnotherNamespace");
    
    // 清空命名空间
    Lists::SetNamespace("");
    assert(Lists::GetNamespace().empty());
    
    std::cout << "PASSED: 命名空间功能测试" << std::endl;
}

/**
 * @brief 测试条目名称获取
 */
void testEntryNames() {
    std::cout << "=== 测试条目名称获取 ===" << std::endl;
    
    auto list = Lists::CreateList("EntryNamesList");
    
    // 添加多个条目
    Lists::AddInteger(list, "First", 1);
    Lists::AddInteger(list, "Second", 2);
    Lists::AddInteger(list, "Third", 3);
    
    // 获取条目名称列表
    auto names = Lists::GetEntryNames(list);
    assert(names.size() == 3);
    
    // 检查名称是否存在（顺序可能不同）
    assert(std::find(names.begin(), names.end(), "First") != names.end());
    assert(std::find(names.begin(), names.end(), "Second") != names.end());
    assert(std::find(names.begin(), names.end(), "Third") != names.end());
    
    std::cout << "PASSED: 条目名称获取测试" << std::endl;
}

/**
 * @brief 测试检查功能
 */
void testCheckFunctionality() {
    std::cout << "=== 测试检查功能 ===" << std::endl;
    
    auto list = Lists::CreateList("CheckList");
    
    // 添加条目
    Lists::AddInteger(list, "Exists", 123);
    
    // 检查存在性
    assert(Lists::Check(list, "Exists"));
    assert(!Lists::Check(list, "NonExistent"));
    
    // 空列表检查
    auto emptyList = Lists::CreateList("EmptyList");
    assert(!Lists::Check(emptyList, "Anything"));
    
    // 空指针检查
    assert(!Lists::Check(nullptr, "Anything"));
    
    std::cout << "PASSED: 检查功能测试" << std::endl;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始Lists模块测试..." << std::endl;
    std::cout << "==================================" << std::endl;
    
    try {
        testListBasicOperations();
        testListDataOperations();
        testListManagement();
        testNamespaceFunctionality();
        testEntryNames();
        testCheckFunctionality();
        
        std::cout << "==================================" << std::endl;
        std::cout << "所有测试通过！Lists模块功能正常。" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cout << "测试失败: " << e.what() << std::endl;
        return 1;
    }
}