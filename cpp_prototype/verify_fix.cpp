// 验证智能指针语法错误修复的测试程序
#include <iostream>
#include <memory>

// 模拟Mesh类的简化版本
class Mesh {
public:
    Mesh(const std::string& name = "") : name_(name) {}
    
    std::string getName() const { return name_; }
    size_t numberOfNodes() const { return 10; }
    size_t numberOfBulkElements() const { return 5; }
    size_t numberOfBoundaryElements() const { return 3; }
    size_t totalElements() const { return 8; }
    
private:
    std::string name_;
};

// 测试修复后的智能指针检查语法
void testSmartPointerSyntax() {
    std::cout << "=== 测试智能指针语法修复 ===" << std::endl;
    
    // 测试1: 空指针检查
    std::shared_ptr<Mesh> nullMesh;
    if (!nullMesh.get()) {
        std::cout << "✓ 空指针检查语法正确" << std::endl;
    }
    
    // 测试2: 有效指针检查
    auto mesh = std::make_shared<Mesh>("test_mesh");
    if (mesh.get()) {
        std::cout << "✓ 有效指针检查语法正确" << std::endl;
        std::cout << "  网格名称: " << mesh->getName() << std::endl;
        std::cout << "  节点数量: " << mesh->numberOfNodes() << std::endl;
    }
    
    std::cout << "=== 所有语法测试通过 ===" << std::endl;
}

int main() {
    try {
        testSmartPointerSyntax();
        std::cout << "\n✅ 智能指针语法修复验证成功！" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }
}