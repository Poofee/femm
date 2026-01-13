/**
 * @file test_mpi_magnetodynamics2d.cpp
 * @brief MPI并行2D磁动力学求解器测试
 * 
 * 测试MPI并行功能在2D磁动力学求解器中的正确性，包括：
 * - 域分解
 * - 并行矩阵组装
 * - 并行求解
 * - 幽灵数据交换
 * - 负载均衡
 */

#include <gtest/gtest.h>
#include <mpi.h>
#include <memory>
#include <vector>
#include <array>

#include "../src/MagnetoDynamics2DSolver.h"
#include "../src/Mesh.h"
#include "../src/MPIUtils.h"
#include "../src/DomainDecomposition.h"
#include "../src/ParallelMatrixAssembly.h"
#include "../src/ParallelLinearSolver.h"

namespace elmer {

class MPIMagnetoDynamics2DSolverTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 初始化MPI（如果尚未初始化）
        MPIUtils::initializeMPI(&argc, &argv);
        
        // 获取MPI通信器
        comm = MPIUtils::getDefaultComm();
        
        // 创建测试网格
        createTestMesh();
    }
    
    void TearDown() override {
        // 清理MPI资源
        MPIUtils::finalizeMPI();
    }
    
    void createTestMesh() {
        // 创建简单的矩形网格用于测试
        mesh = std::make_shared<Mesh>();
        
        // 添加节点
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                double x = static_cast<double>(i) * 0.1;
                double y = static_cast<double>(j) * 0.1;
                mesh->getNodes().addNode({x, y, 0.0});
            }
        }
        
        // 添加四边形单元
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                std::vector<size_t> nodeIndices = {
                    static_cast<size_t>(i * 5 + j),
                    static_cast<size_t>(i * 5 + j + 1),
                    static_cast<size_t>((i + 1) * 5 + j + 1),
                    static_cast<size_t>((i + 1) * 5 + j)
                };
                
                Element element(ElementType::QUADRILATERAL, nodeIndices);
                element.setMaterialName("copper");
                mesh->getBulkElements().push_back(element);
            }
        }
        
        // 设置边界条件
        BoundaryCondition bc;
        bc.setType(BoundaryConditionType::DIRICHLET);
        bc.setValue(0.0);
        
        // 添加边界条件
        mesh->getBoundaryConditions().push_back(bc);
    }
    
    int argc = 0;
    char** argv = nullptr;
    std::shared_ptr<MPICommunicator> comm;
    std::shared_ptr<Mesh> mesh;
};

// 测试域分解功能
TEST_F(MPIMagnetoDynamics2DSolverTest, DomainDecomposition) {
    if (comm->getSize() < 2) {
        GTEST_SKIP() << "需要至少2个MPI进程进行域分解测试";
    }
    
    MagnetoDynamics2DSolver solver(mesh, comm);
    
    // 执行域分解
    auto decomposition = solver.performDomainDecomposition();
    
    // 验证域分解结果
    EXPECT_EQ(decomposition.numPartitions, comm->getSize());
    EXPECT_GT(decomposition.loadBalance, 0.5); // 负载均衡度应大于0.5
    
    // 验证每个分区都有元素
    for (int i = 0; i < decomposition.numPartitions; ++i) {
        EXPECT_GT(decomposition.partitionSizes[i], 0);
    }
    
    // 验证元素分区映射的一致性
    int totalElements = 0;
    for (int i = 0; i < decomposition.numPartitions; ++i) {
        totalElements += decomposition.partitionSizes[i];
    }
    
    EXPECT_EQ(totalElements, 16); // 4x4网格有16个元素
}

// 测试并行矩阵组装
TEST_F(MPIMagnetoDynamics2DSolverTest, ParallelMatrixAssembly) {
    if (comm->getSize() < 2) {
        GTEST_SKIP() << "需要至少2个MPI进程进行并行矩阵组装测试";
    }
    
    MagnetoDynamics2DSolver solver(mesh, comm);
    
    // 设置求解器参数
    MagnetoDynamics2DParameters params;
    params.isTransient = true;
    params.tolerance = 1e-6;
    params.maxIterations = 1000;
    solver.setParameters(params);
    
    // 执行域分解
    auto decomposition = solver.performDomainDecomposition();
    
    // 组装系统矩阵
    solver.assembleSystem();
    
    // 验证系统矩阵是否成功组装
    EXPECT_TRUE(solver.isSystemAssembled());
    
    // 验证分布式矩阵的维度
    auto stiffnessMatrix = solver.getStiffnessMatrix();
    if (stiffnessMatrix) {
        EXPECT_EQ(stiffnessMatrix->NumberOfRows(), 25); // 5x5网格有25个节点
        EXPECT_EQ(stiffnessMatrix->NumberOfColumns(), 25);
    }
}

// 测试并行求解
TEST_F(MPIMagnetoDynamics2DSolverTest, ParallelSolving) {
    if (comm->getSize() < 2) {
        GTEST_SKIP() << "需要至少2个MPI进程进行并行求解测试";
    }
    
    MagnetoDynamics2DSolver solver(mesh, comm);
    
    // 设置求解器参数
    MagnetoDynamics2DParameters params;
    params.isTransient = false; // 稳态分析
    params.tolerance = 1e-6;
    params.maxIterations = 1000;
    solver.setParameters(params);
    
    // 组装系统矩阵
    solver.assembleSystem();
    
    // 求解线性系统
    auto solution = solver.solveLinearSystemParallel();
    
    // 验证求解结果
    EXPECT_FALSE(solution.empty());
    EXPECT_EQ(solution.size(), 25); // 5x5网格有25个节点
    
    // 验证解向量的合理性（不应全为零）
    bool hasNonZero = false;
    for (double val : solution) {
        if (std::abs(val) > 1e-10) {
            hasNonZero = true;
            break;
        }
    }
    EXPECT_TRUE(hasNonZero);
}

// 测试负载均衡
TEST_F(MPIMagnetoDynamics2DSolverTest, LoadBalancing) {
    if (comm->getSize() < 2) {
        GTEST_SKIP() << "需要至少2个MPI进程进行负载均衡测试";
    }
    
    MagnetoDynamics2DSolver solver(mesh, comm);
    
    // 执行域分解
    auto decomposition = solver.performDomainDecomposition();
    
    // 计算负载均衡度
    double loadBalance = solver.computeLoadBalance(decomposition);
    
    // 验证负载均衡度在合理范围内
    EXPECT_GT(loadBalance, 0.0);
    EXPECT_LE(loadBalance, 1.0);
    
    // 对于均匀网格，负载均衡度应接近1.0
    if (comm->getSize() == 2 || comm->getSize() == 4) {
        EXPECT_GT(loadBalance, 0.8);
    }
}

// 测试串行与并行结果的一致性
TEST_F(MPIMagnetoDynamics2DSolverTest, SerialParallelConsistency) {
    // 仅在有多个进程时测试
    if (comm->getSize() < 2) {
        GTEST_SKIP() << "需要至少2个MPI进程进行一致性测试";
    }
    
    // 创建串行求解器
    auto serialComm = MPIUtils::createSerialComm();
    MagnetoDynamics2DSolver serialSolver(mesh, serialComm);
    
    // 创建并行求解器
    MagnetoDynamics2DSolver parallelSolver(mesh, comm);
    
    // 设置相同的参数
    MagnetoDynamics2DParameters params;
    params.isTransient = false;
    params.tolerance = 1e-6;
    params.maxIterations = 1000;
    
    serialSolver.setParameters(params);
    parallelSolver.setParameters(params);
    
    // 串行求解
    serialSolver.assembleSystem();
    auto serialSolution = serialSolver.solveLinearSystem();
    
    // 并行求解
    parallelSolver.assembleSystem();
    auto parallelSolution = parallelSolver.solveLinearSystemParallel();
    
    // 验证解向量维度一致
    EXPECT_EQ(serialSolution.size(), parallelSolution.size());
    
    // 验证解向量的数值一致性（在容差范围内）
    if (comm->getRank() == 0) {
        for (size_t i = 0; i < serialSolution.size(); ++i) {
            double diff = std::abs(serialSolution[i] - parallelSolution[i]);
            EXPECT_LT(diff, 1e-8); // 串行和并行解应非常接近
        }
    }
}

// 测试幽灵数据交换
TEST_F(MPIMagnetoDynamics2DSolverTest, GhostDataExchange) {
    if (comm->getSize() < 2) {
        GTEST_SKIP() << "需要至少2个MPI进程进行幽灵数据交换测试";
    }
    
    MagnetoDynamics2DSolver solver(mesh, comm);
    
    // 执行域分解
    auto decomposition = solver.performDomainDecomposition();
    
    // 获取本地元素和幽灵元素
    auto localElements = solver.getLocalElements(decomposition);
    auto ghostElements = solver.getGhostBoundaryElements(decomposition);
    
    // 验证本地元素不为空
    EXPECT_FALSE(localElements.empty());
    
    // 验证幽灵元素的存在（在边界进程上）
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
        // 边界进程应有幽灵元素
        EXPECT_FALSE(ghostElements.empty());
    }
    
    // 执行幽灵数据交换
    solver.exchangeGhostData();
    
    // 验证幽灵数据交换成功完成（无异常）
    SUCCEED();
}

// 测试MPI通信器集成
TEST_F(MPIMagnetoDynamics2DSolverTest, MPICommunicatorIntegration) {
    MagnetoDynamics2DSolver solver(mesh, comm);
    
    // 验证MPI通信器正确集成
    EXPECT_TRUE(solver.isParallel());
    
    // 验证并行组件正确初始化
    if (comm->getSize() > 1) {
        EXPECT_NE(solver.getDecompositionManager(), nullptr);
        EXPECT_NE(solver.getParallelAssembler(), nullptr);
        EXPECT_NE(solver.getParallelSolver(), nullptr);
    }
}

// 测试性能基准
TEST_F(MPIMagnetoDynamics2DSolverTest, PerformanceBenchmark) {
    if (comm->getSize() < 2) {
        GTEST_SKIP() << "需要至少2个MPI进程进行性能基准测试";
    }
    
    MagnetoDynamics2DSolver solver(mesh, comm);
    
    // 设置求解器参数
    MagnetoDynamics2DParameters params;
    params.isTransient = false;
    params.tolerance = 1e-6;
    params.maxIterations = 1000;
    solver.setParameters(params);
    
    // 测量并行组装时间
    auto startTime = MPI_Wtime();
    solver.assembleSystem();
    auto assemblyTime = MPI_Wtime() - startTime;
    
    // 测量并行求解时间
    startTime = MPI_Wtime();
    auto solution = solver.solveLinearSystemParallel();
    auto solvingTime = MPI_Wtime() - startTime;
    
    // 输出性能信息（仅在主进程）
    if (comm->getRank() == 0) {
        std::cout << "MPI并行性能基准:" << std::endl;
        std::cout << "  进程数: " << comm->getSize() << std::endl;
        std::cout << "  组装时间: " << assemblyTime << " 秒" << std::endl;
        std::cout << "  求解时间: " << solvingTime << " 秒" << std::endl;
        std::cout << "  总时间: " << (assemblyTime + solvingTime) << " 秒" << std::endl;
    }
    
    // 验证性能指标（时间应为正数）
    EXPECT_GT(assemblyTime, 0.0);
    EXPECT_GT(solvingTime, 0.0);
}

} // namespace elmer

// MPI测试主函数
int main(int argc, char** argv) {
    // 初始化Google Test
    ::testing::InitGoogleTest(&argc, argv);
    
    // 初始化MPI
    MPI_Init(&argc, &argv);
    
    // 运行测试
    int result = RUN_ALL_TESTS();
    
    // 清理MPI
    MPI_Finalize();
    
    return result;
}