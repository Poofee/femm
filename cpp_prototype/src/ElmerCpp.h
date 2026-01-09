#pragma once

#include <memory>
#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>

namespace elmer {

// Forward declarations for mesh-related types
struct Node;
enum class ElementType;

// Type aliases for numerical types
using Integer = int;
using Real = double;

/**
 * @brief Vector interface for linear algebra operations
 */
class Vector {
public:
    virtual ~Vector() = default;
    
    virtual Integer Size() const = 0;
    virtual Real& operator[](Integer i) = 0;
    virtual const Real& operator[](Integer i) const = 0;
    virtual void Zero() = 0;
    
    // Factory method for creating vectors
    static std::unique_ptr<Vector> Create(Integer size);
};

/**
 * @brief Matrix interface for linear algebra operations
 */
class Matrix {
public:
    virtual ~Matrix() = default;
    
    virtual void Zero() = 0;
    virtual void SetElement(Integer i, Integer j, Real value) = 0;
    virtual Real GetElement(Integer i, Integer j) const = 0;
    virtual void AddToElement(Integer i, Integer j, Real value) = 0;
    
    // Matrix size information
    virtual Integer GetNumRows() const = 0;
    virtual Integer GetNumCols() const = 0;
    
    // Matrix-vector multiplication
    virtual void Multiply(const Vector& x, Vector& y) const = 0;
    
    // Matrix-vector multiplication with scaling: y = alpha * A * x + beta * y
    virtual void Multiply(Real alpha, const Vector& x, Real beta, Vector& y) const = 0;
    
    // Factory methods for creating different matrix types
    static std::unique_ptr<Matrix> CreateCRS(Integer nrows, Integer ncols);
    static std::unique_ptr<Matrix> CreateBand(Integer nrows, Integer bandwidth);
    static std::unique_ptr<Matrix> CreateDense(Integer nrows, Integer ncols);
};

} // namespace elmer