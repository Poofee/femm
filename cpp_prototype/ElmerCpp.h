// ElmerCpp - C++移植版Elmer FEM主头文件
// 基于Elmer FEM架构设计的现代C++接口

#ifndef ELMER_CPP_H
#define ELMER_CPP_H

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <functional>

namespace elmer {

// 基础类型定义
using Real = double;
using Integer = int;

// 前向声明
class Mesh;
class Solver;
class Matrix;
class Vector;
class Model;
class Variable;

// 核心数据结构
struct Element {
    Integer id;
    std::vector<Integer> nodes;
    std::vector<Real> coordinates;
    // 其他元素属性
};

struct Node {
    Integer id;
    Real x, y, z;
    // 其他节点属性
};

// 并行环境
class ParallelEnvironment {
public:
    static std::shared_ptr<ParallelEnvironment> Init();
    Integer GetMyPE() const { return my_pe_; }
    Integer GetNumPE() const { return num_pe_; }
    
private:
    Integer my_pe_ = 0;
    Integer num_pe_ = 1;
};

// 矩阵接口（对应CRSMatrix.F90）
class Matrix {
public:
    enum class Format { CRS, LIST, BAND };
    
    virtual ~Matrix() = default;
    virtual void Zero() = 0;
    virtual void SetElement(Integer i, Integer j, Real value) = 0;
    virtual Real GetElement(Integer i, Integer j) const = 0;
    virtual void AddToElement(Integer i, Integer j, Real value) = 0;
    
    // 工厂方法
    static std::unique_ptr<Matrix> CreateCRS(Integer nrows, Integer ncols);
    static std::unique_ptr<Matrix> CreateBand(Integer nrows, Integer bandwidth);
};

// 向量接口
class Vector {
public:
    virtual ~Vector() = default;
    virtual Integer Size() const = 0;
    virtual Real& operator[](Integer i) = 0;
    virtual const Real& operator[](Integer i) const = 0;
    virtual void Zero() = 0;
    
    static std::unique_ptr<Vector> Create(Integer size);
};

// 网格类（对应MeshUtils.F90）
class Mesh {
public:
    virtual ~Mesh() = default;
    virtual Integer GetNumberOfNodes() const = 0;
    virtual Integer GetNumberOfElements() const = 0;
    virtual const Node& GetNode(Integer id) const = 0;
    virtual const Element& GetElement(Integer id) const = 0;
    
    // 工厂方法
    static std::unique_ptr<Mesh> LoadFromFile(const std::string& filename);
};

// 求解器基类（对应各种物理场求解器）
class Solver {
public:
    virtual ~Solver() = default;
    virtual void Initialize() = 0;
    virtual void Execute() = 0;
    virtual void Finalize() = 0;
    
    // 求解器注册机制
    static void RegisterSolver(const std::string& name, 
                              std::function<std::unique_ptr<Solver>()> factory);
    static std::unique_ptr<Solver> CreateSolver(const std::string& name);
};

// 主模型类（对应Model_t结构）
class Model {
public:
    Model();
    ~Model();
    
    void LoadSifFile(const std::string& filename);
    void AddSolver(std::unique_ptr<Solver> solver);
    void Solve();
    
    std::shared_ptr<Mesh> GetMesh() const { return mesh_; }
    std::shared_ptr<ParallelEnvironment> GetParallelEnv() const { return parallel_env_; }
    
private:
    std::shared_ptr<Mesh> mesh_;
    std::shared_ptr<ParallelEnvironment> parallel_env_;
    std::vector<std::unique_ptr<Solver>> solvers_;
    std::map<std::string, std::shared_ptr<Variable>> variables_;
};

// 主求解器类（对应ElmerSolver.F90）
class ElmerSolver {
public:
    ElmerSolver();
    
    void Initialize(Integer argc, char** argv);
    void Run();
    void Finalize();
    
    // Fortran互操作性接口
    extern "C" {
        void elmer_solver_init_cpp(int* initialize);
        void elmer_solver_run_cpp();
        void elmer_solver_finalize_cpp();
    }
    
private:
    std::unique_ptr<Model> model_;
    std::shared_ptr<ParallelEnvironment> parallel_env_;
};

} // namespace elmer

#endif // ELMER_CPP_H