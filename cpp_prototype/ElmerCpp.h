// ElmerCpp - C++ port of Elmer FEM main header file
// Modern C++ interface based on Elmer FEM architecture

#ifndef ELMER_CPP_H
#define ELMER_CPP_H

#include <memory>
#include <vector>
#include <string>
#include <map>
#include <functional>

namespace elmer {

// Basic type definitions
using Real = double;
using Integer = int;

// Forward declarations
class Mesh;
class Solver;
class Matrix;
class Vector;
class Model;
class Variable;

// Core data structures
struct Element {
    Integer id;
    std::vector<Integer> nodes;
    std::vector<Real> coordinates;
};

struct Node {
    Integer id;
    Real x, y, z;
};

// Parallel environment
class ParallelEnvironment {
public:
    static std::shared_ptr<ParallelEnvironment> Init();
    Integer GetMyPE() const { return my_pe_; }
    Integer GetNumPE() const { return num_pe_; }
    
private:
    Integer my_pe_ = 0;
    Integer num_pe_ = 1;
};

// Matrix interface (corresponds to CRSMatrix.F90)
class Matrix {
public:
    enum class Format { CRS, LIST, BAND };
    
    virtual ~Matrix() = default;
    virtual void Zero() = 0;
    virtual void SetElement(Integer i, Integer j, Real value) = 0;
    virtual Real GetElement(Integer i, Integer j) const = 0;
    virtual void AddToElement(Integer i, Integer j, Real value) = 0;
    
    // Factory methods
    static std::unique_ptr<Matrix> CreateCRS(Integer nrows, Integer ncols);
    static std::unique_ptr<Matrix> CreateBand(Integer nrows, Integer bandwidth);
};

// Vector interface
class Vector {
public:
    virtual ~Vector() = default;
    virtual Integer Size() const = 0;
    virtual Real& operator[](Integer i) = 0;
    virtual const Real& operator[](Integer i) const = 0;
    virtual void Zero() = 0;
    
    static std::unique_ptr<Vector> Create(Integer size);
};

// Mesh class (corresponds to MeshUtils.F90)
class Mesh {
public:
    virtual ~Mesh() = default;
    virtual Integer GetNumberOfNodes() const = 0;
    virtual Integer GetNumberOfElements() const = 0;
    virtual const Node& GetNode(Integer id) const = 0;
    virtual const Element& GetElement(Integer id) const = 0;
    
    static std::unique_ptr<Mesh> LoadFromFile(const std::string& filename);
};

// Solver base class (corresponds to various physics solvers)
class Solver {
public:
    virtual ~Solver() = default;
    virtual void Initialize() = 0;
    virtual void Execute() = 0;
    virtual void Finalize() = 0;
    
    virtual void SetName(const std::string& name) = 0;
    virtual std::string GetName() const = 0;
};

// Model class (main simulation container)
class Model {
public:
    virtual ~Model() = default;
    virtual void AddSolver(std::unique_ptr<Solver> solver) = 0;
    virtual void SetMesh(std::unique_ptr<Mesh> mesh) = 0;
    virtual void Run() = 0;
    
    static std::unique_ptr<Model> Create();
};

// Variable class (field variables)
class Variable {
public:
    virtual ~Variable() = default;
    virtual void SetName(const std::string& name) = 0;
    virtual std::string GetName() const = 0;
    virtual void SetValues(const std::vector<Real>& values) = 0;
    virtual const std::vector<Real>& GetValues() const = 0;
};

// Main solver class (entry point)
class ElmerSolver {
public:
    ElmerSolver();
    
    void Initialize(Integer argc, char** argv);
    void Run();
    void Finalize();
    
private:
    std::unique_ptr<Model> model_;
    std::shared_ptr<ParallelEnvironment> parallel_env_;
};

} // namespace elmer

#endif // ELMER_CPP_H