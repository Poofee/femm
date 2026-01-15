#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <cassert>
#include "ElmerSolver.h"
#include "SolverRegistry.h"
#include "SolverBase.h"
#include "Mesh.h"
#include "Material.h"
#include "LinearAlgebra.h"
#include "BoundaryConditions.h"

using namespace elmer;

/**
 * @brief 简单的测试求解器类
 */
class TestSolver : public SolverBase {
private:
    std::vector<double> solution_;
    
public:
    TestSolver(const std::string& name = "TestSolver") 
        : SolverBase(name) {
        solution_ = {1.0, 2.0, 3.0}; // 简单的测试解
    }
    
    ~TestSolver() override = default;
    
    bool initialize() override {
        std::cout << "TestSolver " << name_ << " 初始化成功" << std::endl;
        status_ = SolverStatus::INITIALIZED;
        return true;
    }
    
    bool assemble() override {
        std::cout << "TestSolver " << name_ << " 组装完成" << std::endl;
        status_ = SolverStatus::ASSEMBLED;
        return true;
    }
    
    bool solve() override {
        std::cout << "TestSolver " << name_ << " 求解完成" << std::endl;
        status_ = SolverStatus::SOLVED;
        return true;
    }
    
    std::vector<double> getSolution() const override {
        return solution_;
    }
    
    void cleanup() override {
        std::cout << "TestSolver " << name_ << " 清理完成" << std::endl;
        status_ = SolverStatus::FINISHED;
    }
};

/**
 * @brief 线性测试求解器类
 */
class LinearTestSolver : public LinearSolverBase {
private:
    std::vector<double> solution_;
    
public:
    LinearTestSolver(const std::string& name = "LinearTestSolver") 
        : LinearSolverBase(name) {
        solution_ = {4.0, 5.0, 6.0}; // 简单的测试解
    }
    
    ~LinearTestSolver() override = default;
    
    bool initialize() override {
        std::cout << "LinearTestSolver " << name_ << " 初始化成功" << std::endl;
        status_ = SolverStatus::INITIALIZED;
        return true;
    }
    
    bool assemble() override {
        // 创建简单的刚度矩阵和右端向量
        stiffnessMatrix_ = std::make_shared<Matrix>(3, 3);
        rhsVector_ = std::make_shared<Vector>(3);
        solution_ = std::make_shared<Vector>(3);
        
        // 简单的对角矩阵
        for (int i = 0; i < 3; ++i) {
            stiffnessMatrix_->set(i, i, 2.0);
            rhsVector_->set(i, 1.0);
        }
        
        std::cout << "LinearTestSolver " << name_ << " 组装完成" << std::endl;
        status_ = SolverStatus::ASSEMBLED;
        return true;
    }
    
    bool solve() override {
        // 简单的求解：x = A^(-1) * b
        for (int i = 0; i < 3; ++i) {
            double value = rhsVector_->get(i) / stiffnessMatrix_->get(i, i);
            solution_->set(i, value);
        }
        
        std::cout << "LinearTestSolver " << name_ << " 求解完成" << std::endl;
        status_ = SolverStatus::SOLVED;
        return true;
    }
    
    std::vector<double> getSolution() const override {
        std::vector<double> result;
        for (int i = 0; i < 3; ++i) {
            result.push_back(solution_->get(i));
        }
        return result;
    }
    
    void cleanup() override {
        std::cout << "LinearTestSolver " << name_ << " 清理完成" << std::endl;
        status_ = SolverStatus::FINISHED;
    }
};

/**
 * @brief 创建简单的测试网格
 */
std::shared_ptr<Mesh> createSimpleTestMesh() {
    auto mesh = std::make_shared<Mesh>();
    
    // 创建4个节点
    std::vector<Node> nodes = {
        Node(0.0, 0.0, 0.0),
        Node(1.0, 0.0, 0.0),
        Node(1.0, 1.0, 0.0),
        Node(0.0, 1.0, 0.0)
    };
    
    // 创建单个四边形单元
    std::vector<size_t> elementNodes = {0, 1, 2, 3};
    Element element(ElementType::QUADRATIC, 0);
    element.setNodeIndices(elementNodes);
    
    // 添加节点到网格
    for (const auto& node : nodes) {
        mesh->getNodes().addNode(node);
    }
    
    // 添加单元到网格
    mesh->addBulkElement(element);
    
    return mesh;
}

/**
 * @brief 创建材料数据库
 */
elmer::MaterialDatabase createSimpleMaterialDatabase() {
    elmer::MaterialDatabase materialDB;
    
    // 添加测试材料
    auto material = std::make_shared<Material>("TestMaterial");
    material->setProperty("Density", 1000.0);
    material->setProperty("YoungsModulus", 2.1e11);
    material->setProperty("PoissonsRatio", 0.3);
    
    materialDB.addMaterial(material);
    
    return materialDB;
}

/**
 * @brief 测试求解器注册机制
 */
void testSolverRegistry() {
    std::cout << "\n=== 测试求解器注册机制 ===" << std::endl;
    
    auto& registry = SolverRegistry::getInstance();
    
    // 注册测试求解器
    registry.registerSolver("TestSolver", []() -> std::shared_ptr<SolverBase> {
        return std::make_shared<TestSolver>("RegisteredTestSolver");
    });
    
    registry.registerSolver("LinearTestSolver", []() -> std::shared_ptr<SolverBase> {
        return std::make_shared<LinearTestSolver>("RegisteredLinearTestSolver");
    });
    
    // 测试求解器创建
    auto testSolver = registry.createSolver("TestSolver");
    assert(testSolver != nullptr);
    assert(testSolver->getName() == "RegisteredTestSolver");
    
    auto linearSolver = registry.createSolver("LinearTestSolver");
    assert(linearSolver != nullptr);
    assert(linearSolver->getName() == "RegisteredLinearTestSolver");
    
    // 测试不存在的求解器
    auto invalidSolver = registry.createSolver("NonExistentSolver");
    assert(invalidSolver == nullptr);
    
    std::cout << "✓ 求解器注册机制测试通过" << std::endl;
}

/**
 * @brief 测试求解器管理器
 */
void testSolverManager() {
    std::cout << "\n=== 测试求解器管理器 ===" << std::endl;
    
    auto mesh = createSimpleTestMesh();
    auto materialDB = createSimpleMaterialDatabase();
    auto bc = std::make_shared<BoundaryConditions>();
    
    SolverManager manager;
    manager.setMesh(mesh);
    manager.setMaterialDatabase(materialDB);
    manager.setBoundaryConditions(bc);
    
    // 添加求解器
    auto solver1 = std::make_shared<TestSolver>("Solver1");
    auto solver2 = std::make_shared<LinearTestSolver>("Solver2");
    
    manager.addSolver(solver1);
    manager.addSolver(solver2);
    
    // 测试初始化
    assert(manager.initialize());
    
    // 测试执行
    assert(manager.executeAll());
    
    // 测试获取解
    auto solutions = manager.getSolutions();
    assert(solutions.size() == 2);
    assert(solutions.find("Solver1") != solutions.end());
    assert(solutions.find("Solver2") != solutions.end());
    
    // 验证解的正确性
    auto sol1 = solutions["Solver1"];
    auto sol2 = solutions["Solver2"];
    
    assert(sol1.size() == 3);
    assert(sol2.size() == 3);
    
    // 清理
    manager.cleanup();
    
    std::cout << "✓ 求解器管理器测试通过" << std::endl;
}

/**
 * @brief 测试主求解器稳态仿真
 */
void testElmerSolverSteadyState() {
    std::cout << "\n=== 测试主求解器稳态仿真 ===" << std::endl;
    
    ElmerSolver solver;
    
    // 设置仿真参数
    SimulationParameters params;
    params.type = SimulationType::STEADY_STATE;
    params.maxIterations = 10;
    params.tolerance = 1e-6;
    
    solver.setParameters(params);
    
    // 注册测试求解器
    auto& registry = SolverRegistry::getInstance();
    registry.registerSolver("TestSolver", []() -> std::shared_ptr<SolverBase> {
        return std::make_shared<TestSolver>("TestSolver");
    });
    
    // 添加求解器
    solver.addSolver("TestSolver");
    
    // 执行稳态仿真
    auto result = solver.execute();
    
    // 验证结果
    assert(result.success);
    assert(result.realTime > 0);
    assert(result.solutions.size() == 1);
    
    std::cout << "✓ 主求解器稳态仿真测试通过" << std::endl;
    std::cout << "  仿真耗时: " << result.realTime << " 秒" << std::endl;
}

/**
 * @brief 测试主求解器瞬态仿真
 */
void testElmerSolverTransient() {
    std::cout << "\n=== 测试主求解器瞬态仿真 ===" << std::endl;
    
    ElmerSolver solver;
    
    // 设置仿真参数
    SimulationParameters params;
    params.type = SimulationType::TRANSIENT;
    params.startTime = 0.0;
    params.endTime = 1.0;
    params.timeStep = 0.1;
    params.maxIterations = 5;
    params.tolerance = 1e-6;
    
    solver.setParameters(params);
    
    // 注册测试求解器
    auto& registry = SolverRegistry::getInstance();
    registry.registerSolver("LinearTestSolver", []() -> std::shared_ptr<SolverBase> {
        return std::make_shared<LinearTestSolver>("LinearTestSolver");
    });
    
    // 添加求解器
    solver.addSolver("LinearTestSolver");
    
    // 执行瞬态仿真
    auto result = solver.execute();
    
    // 验证结果
    assert(result.success);
    assert(result.realTime > 0);
    assert(result.solutions.size() == 1);
    
    std::cout << "✓ 主求解器瞬态仿真测试通过" << std::endl;
    std::cout << "  仿真耗时: " << result.realTime << " 秒" << std::endl;
}

/**
 * @brief 测试MPI通信功能（串行模式）
 */
void testMPICommunicatorSerial() {
    std::cout << "\n=== 测试MPI通信功能（串行模式） ===" << std::endl;
    
    MPICommunicator comm;
    
    // 测试串行模式初始化
    bool initSuccess = comm.initialize();
    assert(initSuccess);
    
    // 测试进程信息
    assert(comm.getRank() == 0);
    assert(comm.getSize() == 1);
    
    // 测试数据交换（串行模式）
    std::vector<double> sendData = {1.0, 2.0, 3.0};
    std::vector<double> recvData(3);
    
    comm.allReduce(sendData.data(), recvData.data(), 3, "sum");
    
    // 在串行模式下，数据应该保持不变
    for (int i = 0; i < 3; ++i) {
        assert(recvData[i] == sendData[i]);
    }
    
    // 测试广播（串行模式）
    std::vector<double> broadcastData = {4.0, 5.0, 6.0};
    comm.broadcast(broadcastData.data(), 3, 0);
    
    // 数据应该保持不变
    for (int i = 0; i < 3; ++i) {
        assert(broadcastData[i] == 4.0 + i);
    }
    
    std::cout << "✓ MPI通信功能串行模式测试通过" << std::endl;
}

/**
 * @brief 测试求解器基类功能
 */
void testSolverBaseFunctionality() {
    std::cout << "\n=== 测试求解器基类功能 ===" << std::endl;
    
    TestSolver solver("BaseTestSolver");
    
    // 测试基本属性
    assert(solver.getName() == "BaseTestSolver");
    assert(solver.getStatus() == SolverStatus::INITIALIZED);
    
    // 测试参数设置
    SolverParameters params;
    params.maxIterations = 100;
    params.tolerance = 1e-8;
    params.linearSolverType = "CG";
    
    solver.setParameters(params);
    auto retrievedParams = solver.getParameters();
    
    assert(retrievedParams.maxIterations == 100);
    assert(retrievedParams.tolerance == 1e-8);
    assert(retrievedParams.linearSolverType == "CG");
    
    // 测试求解器生命周期
    assert(solver.initialize());
    assert(solver.getStatus() == SolverStatus::INITIALIZED);
    
    assert(solver.assemble());
    assert(solver.getStatus() == SolverStatus::ASSEMBLED);
    
    assert(solver.solve());
    assert(solver.getStatus() == SolverStatus::SOLVED);
    
    auto solution = solver.getSolution();
    assert(solution.size() == 3);
    assert(solution[0] == 1.0);
    assert(solution[1] == 2.0);
    assert(solution[2] == 3.0);
    
    solver.cleanup();
    assert(solver.getStatus() == SolverStatus::FINISHED);
    
    std::cout << "✓ 求解器基类功能测试通过" << std::endl;
}

/**
 * @brief 主测试函数
 */
int main() {
    std::cout << "开始测试Elmer求解器功能..." << std::endl;
    
    try {
        testSolverRegistry();
        testSolverManager();
        testSolverBaseFunctionality();
        testMPICommunicatorSerial();
        testElmerSolverSteadyState();
        testElmerSolverTransient();
        
        std::cout << "\n🎉 所有测试通过！Elmer求解器功能验证成功。" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "❌ 测试失败: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "❌ 未知错误导致测试失败" << std::endl;
        return 1;
    }
}