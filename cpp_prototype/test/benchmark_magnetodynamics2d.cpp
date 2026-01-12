/**
 * @file benchmark_magnetodynamics2d.cpp
 * @brief 2D磁动力学求解器性能基准测试
 * @author Elmer C++ Porting Team
 * @date 2026-01-12
 * 
 * 用于对比性能优化前后的计算时间差异
 */

#include <iostream>
#include <chrono>
#include <vector>
#include <memory>
#include "../src/MagnetoDynamics2DSolver.h"
#include "../src/Mesh.h"
#include "../src/LinearAlgebra.h"

namespace elmer {

/**
 * @brief 性能基准测试类
 */
class PerformanceBenchmark {
private:
    std::shared_ptr<MagnetoDynamics2DSolver> solver;
    
public:
    PerformanceBenchmark() {
        solver = std::make_shared<MagnetoDynamics2DSolver>();
    }
    
    /**
     * @brief 创建测试网格
     */
    std::shared_ptr<Mesh> createTestMesh(int numElements = 100) {
        auto mesh = std::make_shared<Mesh>("BenchmarkTestMesh");
        
        // 创建简化的测试网格 - 2x2 网格
        // 添加4个节点
        mesh->getNodes().addNode(0.0, 0.0, 0.0); // 节点0
        mesh->getNodes().addNode(1.0, 0.0, 0.0); // 节点1
        mesh->getNodes().addNode(0.0, 1.0, 0.0); // 节点2
        mesh->getNodes().addNode(1.0, 1.0, 0.0); // 节点3
        
        // 创建2个三角形单元
        Element elem1(ElementType::LINEAR, 1);
        elem1.addNodeIndex(0); // 节点0
        elem1.addNodeIndex(1); // 节点1
        elem1.addNodeIndex(2); // 节点2
        mesh->addBulkElement(elem1);
        
        Element elem2(ElementType::LINEAR, 1);
        elem2.addNodeIndex(1); // 节点1
        elem2.addNodeIndex(3); // 节点3
        elem2.addNodeIndex(2); // 节点2
        mesh->addBulkElement(elem2);
        
        // 初始化并行信息（简化版本）
        mesh->getParallelInfo().initialize(mesh->numberOfNodes());
        
        return mesh;
    }
    
    /**
     * @brief 基准测试：传统矩阵组装
     */
    double benchmarkTraditionalAssembly(int iterations = 10) {
        auto mesh = createTestMesh();
        solver->setMesh(mesh);
        
        // 禁用缓存
        solver->enableBasisFunctionCache(false);
        
        auto start = std::chrono::high_resolution_clock::now();
        
        for (int i = 0; i < iterations; ++i) {
            solver->assembleSystem();
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        return duration.count() / 1000.0; // 返回毫秒
    }
    
    /**
     * @brief 基准测试：优化矩阵组装
     */
    double benchmarkOptimizedAssembly(int iterations = 10) {
        auto mesh = createTestMesh();
        solver->setMesh(mesh);
        
        // 启用缓存并预计算
        solver->enableBasisFunctionCache(true);
        solver->precomputeBasisFunctionCache();
        
        auto start = std::chrono::high_resolution_clock::now();
        
        for (int i = 0; i < iterations; ++i) {
            solver->assembleSystemOptimized();
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        return duration.count() / 1000.0; // 返回毫秒
    }
    
    /**
     * @brief 基准测试：并行矩阵组装
     */
    double benchmarkParallelAssembly(int iterations = 10, int numThreads = 4) {
        auto mesh = createTestMesh();
        solver->setMesh(mesh);
        
        // 设置并行线程数
        solver->setParallelThreads(numThreads);
        
        auto start = std::chrono::high_resolution_clock::now();
        
        for (int i = 0; i < iterations; ++i) {
            solver->assembleSystemParallel();
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        return duration.count() / 1000.0; // 返回毫秒
    }
    
    /**
     * @brief 基准测试：不同线程数的并行性能
     */
    void benchmarkParallelScaling(int iterations = 5) {
        std::cout << "\n7. 并行扩展性基准测试" << std::endl;
        
        // 测试不同线程数的性能
        std::vector<int> threadCounts = {1, 2, 4, 8};
        
        for (int threads : threadCounts) {
            double parallelTime = benchmarkParallelAssembly(iterations, threads);
            double avgTime = parallelTime / iterations;
            
            std::cout << "   " << threads << " 线程 - 平均时间: " << avgTime << " ms/次" << std::endl;
        }
    }
    
    /**
     * @brief 基准测试：缓存预计算时间
     */
    double benchmarkCachePrecomputation() {
        auto mesh = createTestMesh();
        solver->setMesh(mesh);
        
        auto start = std::chrono::high_resolution_clock::now();
        
        solver->precomputeBasisFunctionCache();
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        return duration.count() / 1000.0; // 返回毫秒
    }
    
    /**
     * @brief 运行完整的性能基准测试
     */
    void runComprehensiveBenchmark() {
        std::cout << "========================================" << std::endl;
        std::cout << "2D磁动力学求解器性能基准测试" << std::endl;
        std::cout << "========================================" << std::endl;
        
        const int iterations = 5;
        const int warmupRuns = 2;
        
        // 预热运行
        std::cout << "\n预热运行..." << std::endl;
        for (int i = 0; i < warmupRuns; ++i) {
            benchmarkTraditionalAssembly(1);
            benchmarkOptimizedAssembly(1);
            benchmarkParallelAssembly(1, 4);
        }
        
        // 缓存预计算基准测试
        std::cout << "\n1. 缓存预计算基准测试" << std::endl;
        double cacheTime = benchmarkCachePrecomputation();
        std::cout << "   缓存预计算时间: " << cacheTime << " ms" << std::endl;
        
        // 传统组装基准测试
        std::cout << "\n2. 传统矩阵组装基准测试 (" << iterations << "次迭代)" << std::endl;
        double traditionalTime = benchmarkTraditionalAssembly(iterations);
        std::cout << "   总时间: " << traditionalTime << " ms" << std::endl;
        std::cout << "   平均时间: " << traditionalTime / iterations << " ms/次" << std::endl;
        
        // 优化组装基准测试
        std::cout << "\n3. 优化矩阵组装基准测试 (" << iterations << "次迭代)" << std::endl;
        double optimizedTime = benchmarkOptimizedAssembly(iterations);
        std::cout << "   总时间: " << optimizedTime << " ms" << std::endl;
        std::cout << "   平均时间: " << optimizedTime / iterations << " ms/次" << std::endl;
        
        // 并行组装基准测试
        std::cout << "\n4. 并行矩阵组装基准测试 (" << iterations << "次迭代, 4线程)" << std::endl;
        double parallelTime = benchmarkParallelAssembly(iterations, 4);
        std::cout << "   总时间: " << parallelTime << " ms" << std::endl;
        std::cout << "   平均时间: " << parallelTime / iterations << " ms/次" << std::endl;
        
        // 性能对比分析
        std::cout << "\n5. 性能对比分析" << std::endl;
        double speedupOptimized = traditionalTime / optimizedTime;
        double speedupParallel = traditionalTime / parallelTime;
        
        std::cout << "   优化加速比: " << speedupOptimized << "x" << std::endl;
        std::cout << "   并行加速比: " << speedupParallel << "x" << std::endl;
        
        if (speedupOptimized > 1.0) {
            std::cout << "   优化性能提升: +" << ((speedupOptimized - 1.0) * 100) << "%" << std::endl;
        }
        
        if (speedupParallel > 1.0) {
            std::cout << "   并行性能提升: +" << ((speedupParallel - 1.0) * 100) << "%" << std::endl;
        }
        
        // 缓存效率分析
        std::cout << "\n6. 缓存效率分析" << std::endl;
        double amortizedCacheTime = cacheTime / iterations;
        double effectiveOptimizedTime = optimizedTime + cacheTime;
        double effectiveSpeedup = traditionalTime / effectiveOptimizedTime;
        
        std::cout << "   分摊缓存时间: " << amortizedCacheTime << " ms/次" << std::endl;
        std::cout << "   有效优化时间: " << effectiveOptimizedTime << " ms" << std::endl;
        std::cout << "   有效加速比: " << effectiveSpeedup << "x" << std::endl;
        
        std::cout << "\n========================================" << std::endl;
        std::cout << "基准测试完成" << std::endl;
        std::cout << "========================================" << std::endl;
    }
    
    /**
     * @brief 内存使用基准测试
     */
    void benchmarkMemoryUsage() {
        std::cout << "\n6. 内存使用基准测试" << std::endl;
        
        auto mesh = createTestMesh();
        solver->setMesh(mesh);
        
        // 测试无缓存状态
        solver->enableBasisFunctionCache(false);
        solver->clearCache();
        std::cout << "   无缓存状态: 缓存已清除" << std::endl;
        
        // 测试有缓存状态
        solver->enableBasisFunctionCache(true);
        solver->precomputeBasisFunctionCache();
        std::cout << "   有缓存状态: 缓存已预计算" << std::endl;
        
        // 获取缓存统计信息
        solver->getCacheStatistics();
    }
};

} // namespace elmer

/**
 * @brief 主函数 - 运行性能基准测试
 */
int main() {
    try {
        elmer::PerformanceBenchmark benchmark;
        
        // 运行综合性能基准测试
        benchmark.runComprehensiveBenchmark();
        
        // 运行并行扩展性测试
        benchmark.benchmarkParallelScaling();
        
        // 运行内存使用基准测试
        benchmark.benchmarkMemoryUsage();
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cout << "基准测试过程中出现异常: " << e.what() << std::endl;
        return 1;
    }
}