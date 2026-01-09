// IterativeSolver.h - Iterative Linear System Solver
// Corresponds to Elmer FEM's IterSolve.F90

#pragma once

#include "ElmerCpp.h"
#include "CRSMatrix.h"
#include <memory>
#include <vector>
#include <cmath>
#include <iostream>

namespace elmer {

/**
 * @brief Base class for iterative solvers
 */
class IterativeSolver {
protected:
    Integer max_iterations_;      ///< Maximum number of iterations
    Real tolerance_;              ///< Convergence tolerance
    Integer iteration_count_;     ///< Current iteration count
    Real residual_norm_;          ///< Current residual norm
    bool converged_;              ///< Convergence status
    
public:
    IterativeSolver(Integer max_iter = 1000, Real tol = 1e-8)
        : max_iterations_(max_iter), tolerance_(tol), 
          iteration_count_(0), residual_norm_(0.0), converged_(false) {}
    
    virtual ~IterativeSolver() = default;
    
    // Solver configuration
    void SetMaxIterations(Integer max_iter) { max_iterations_ = max_iter; }
    void SetTolerance(Real tol) { tolerance_ = tol; }
    
    // Solver status
    Integer GetIterationCount() const { return iteration_count_; }
    Real GetResidualNorm() const { return residual_norm_; }
    bool IsConverged() const { return converged_; }
    
    // Main solve method
    virtual bool Solve(const Matrix& A, Vector& x, const Vector& b) = 0;
    
    // Preconditioned solve
    virtual bool Solve(const Matrix& A, Vector& x, const Vector& b, 
                      const Matrix& preconditioner) {
        // Default implementation: ignore preconditioner
        return Solve(A, x, b);
    }
    
    // Reset solver state
    virtual void Reset() {
        iteration_count_ = 0;
        residual_norm_ = 0.0;
        converged_ = false;
    }
};

/**
 * @brief Conjugate Gradient solver for symmetric positive definite systems
 */
class ConjugateGradientSolver : public IterativeSolver {
private:
    std::vector<Real> r_;  ///< Residual vector
    std::vector<Real> p_;  ///< Search direction
    std::vector<Real> Ap_; ///< Matrix-vector product A*p
    
public:
    ConjugateGradientSolver(Integer max_iter = 1000, Real tol = 1e-8)
        : IterativeSolver(max_iter, tol) {}
    
    bool Solve(const Matrix& A, Vector& x, const Vector& b) override;
    
    // Preconditioned CG (not implemented yet)
    bool Solve(const Matrix& A, Vector& x, const Vector& b, 
               const Matrix& preconditioner) override;
};

/**
 * @brief Generalized Minimal Residual solver for non-symmetric systems
 */
class GMRESSolver : public IterativeSolver {
private:
    Integer restart_;  ///< Restart parameter
    
public:
    GMRESSolver(Integer max_iter = 1000, Real tol = 1e-8, Integer restart = 30)
        : IterativeSolver(max_iter, tol), restart_(restart) {}
    
    bool Solve(const Matrix& A, Vector& x, const Vector& b) override;
    
    void SetRestart(Integer restart) { restart_ = restart; }
};

/**
 * @brief Bi-Conjugate Gradient Stabilized solver
 */
class BiCGSTABSolver : public IterativeSolver {
public:
    BiCGSTABSolver(Integer max_iter = 1000, Real tol = 1e-8)
        : IterativeSolver(max_iter, tol) {}
    
    bool Solve(const Matrix& A, Vector& x, const Vector& b) override;
};

/**
 * @brief Jacobi iterative solver (simple but robust)
 */
class JacobiSolver : public IterativeSolver {
private:
    Real relaxation_;  ///< Relaxation parameter
    
public:
    JacobiSolver(Integer max_iter = 1000, Real tol = 1e-8, Real relax = 1.0)
        : IterativeSolver(max_iter, tol), relaxation_(relax) {}
    
    bool Solve(const Matrix& A, Vector& x, const Vector& b) override;
    
    void SetRelaxation(Real relax) { relaxation_ = relax; }
};

/**
 * @brief Gauss-Seidel iterative solver
 */
class GaussSeidelSolver : public IterativeSolver {
private:
    Real relaxation_;  ///< Relaxation parameter
    
public:
    GaussSeidelSolver(Integer max_iter = 1000, Real tol = 1e-8, Real relax = 1.0)
        : IterativeSolver(max_iter, tol), relaxation_(relax) {}
    
    bool Solve(const Matrix& A, Vector& x, const Vector& b) override;
    
    void SetRelaxation(Real relax) { relaxation_ = relax; }
};

/**
 * @brief Utility functions for iterative solvers
 */
namespace IterativeSolverUtils {
    
    /**
     * @brief Compute residual vector r = b - A*x
     */
    void ComputeResidual(const Matrix& A, const Vector& x, const Vector& b, Vector& r);
    
    /**
     * @brief Compute L2 norm of vector
     */
    Real ComputeNorm(const Vector& v);
    
    /**
     * @brief Compute dot product of two vectors
     */
    Real DotProduct(const Vector& a, const Vector& b);
    
    /**
     * @brief Vector addition: y = alpha*x + y
     */
    void VectorAdd(Real alpha, const Vector& x, Vector& y);
    
    /**
     * @brief Vector scaling and addition: y = alpha*x + beta*y
     */
    void VectorAdd(Real alpha, const Vector& x, Real beta, Vector& y);
    
    /**
     * @brief Vector copy: y = x
     */
    void VectorCopy(const Vector& x, Vector& y);
    
    /**
     * @brief Check if matrix is symmetric positive definite (for CG)
     */
    bool IsSymmetricPositiveDefinite(const Matrix& A, Real tolerance = 1e-12);
    
} // namespace IterativeSolverUtils

} // namespace elmer