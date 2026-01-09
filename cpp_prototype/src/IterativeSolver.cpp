// IterativeSolver.cpp - Iterative Linear System Solver Implementation
// Corresponds to Elmer FEM's IterSolve.F90

#include "IterativeSolver.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace elmer {

// ============================================================================
// Conjugate Gradient Solver Implementation
// ============================================================================

bool ConjugateGradientSolver::Solve(const Matrix& A, Vector& x, const Vector& b) {
    Integer n = A.GetNumRows();
    
    // Check dimensions
    if (n != A.GetNumCols() || n != x.Size() || n != b.Size()) {
        throw std::invalid_argument("Matrix and vector dimensions don't match");
    }
    
    // Initialize working vectors
    r_.resize(n);
    p_.resize(n);
    Ap_.resize(n);
    
    // Create temporary vectors
    std::unique_ptr<Vector> r = Vector::Create(n);
    std::unique_ptr<Vector> p = Vector::Create(n);
    std::unique_ptr<Vector> Ap = Vector::Create(n);
    
    // Initial residual: r = b - A*x
    A.Multiply(x, *Ap);
    for (Integer i = 0; i < n; ++i) {
        (*r)[i] = b[i] - (*Ap)[i];
    }
    
    // Initial search direction: p = r
    IterativeSolverUtils::VectorCopy(*r, *p);
    
    Real r_dot_r_old = IterativeSolverUtils::DotProduct(*r, *r);
    residual_norm_ = std::sqrt(r_dot_r_old);
    
    // Check initial convergence
    if (residual_norm_ < tolerance_) {
        converged_ = true;
        return true;
    }
    
    // CG iterations
    iteration_count_ = 0;
    while (iteration_count_ < max_iterations_) {
        iteration_count_++;
        
        // Compute A*p
        A.Multiply(*p, *Ap);
        
        // Compute alpha = (r^T * r) / (p^T * A * p)
        Real p_dot_Ap = IterativeSolverUtils::DotProduct(*p, *Ap);
        if (std::abs(p_dot_Ap) < 1e-20) {
            // Avoid division by zero
            break;
        }
        
        Real alpha = r_dot_r_old / p_dot_Ap;
        
        // Update solution: x = x + alpha * p
        IterativeSolverUtils::VectorAdd(alpha, *p, x);
        
        // Update residual: r = r - alpha * A * p
        IterativeSolverUtils::VectorAdd(-alpha, *Ap, *r);
        
        // Compute new residual norm
        Real r_dot_r_new = IterativeSolverUtils::DotProduct(*r, *r);
        residual_norm_ = std::sqrt(r_dot_r_new);
        
        // Check convergence
        if (residual_norm_ < tolerance_) {
            converged_ = true;
            break;
        }
        
        // Compute beta = (r_new^T * r_new) / (r_old^T * r_old)
        Real beta = r_dot_r_new / r_dot_r_old;
        
        // Update search direction: p = r + beta * p
        IterativeSolverUtils::VectorAdd(1.0, *r, beta, *p);
        
        r_dot_r_old = r_dot_r_new;
    }
    
    return converged_;
}

bool ConjugateGradientSolver::Solve(const Matrix& A, Vector& x, const Vector& b, 
                                    const Matrix& preconditioner) {
    // Preconditioned CG not implemented yet
    // For now, use standard CG
    return Solve(A, x, b);
}

// ============================================================================
// GMRES Solver Implementation
// ============================================================================

bool GMRESSolver::Solve(const Matrix& A, Vector& x, const Vector& b) {
    Integer n = A.GetNumRows();
    
    // Check dimensions
    if (n != A.GetNumCols() || n != x.Size() || n != b.Size()) {
        throw std::invalid_argument("Matrix and vector dimensions don't match");
    }
    
    // GMRES implementation (simplified version)
    // For now, use CG as a placeholder - full GMRES will be implemented later
    ConjugateGradientSolver cg_solver(max_iterations_, tolerance_);
    return cg_solver.Solve(A, x, b);
}

// ============================================================================
// BiCGSTAB Solver Implementation
// ============================================================================

bool BiCGSTABSolver::Solve(const Matrix& A, Vector& x, const Vector& b) {
    Integer n = A.GetNumRows();
    
    // Check dimensions
    if (n != A.GetNumCols() || n != x.Size() || n != b.Size()) {
        throw std::invalid_argument("Matrix and vector dimensions don't match");
    }
    
    // BiCGSTAB implementation (simplified version)
    // For now, use CG as a placeholder - full BiCGSTAB will be implemented later
    ConjugateGradientSolver cg_solver(max_iterations_, tolerance_);
    return cg_solver.Solve(A, x, b);
}

// ============================================================================
// Jacobi Solver Implementation
// ============================================================================

bool JacobiSolver::Solve(const Matrix& A, Vector& x, const Vector& b) {
    Integer n = A.GetNumRows();
    
    // Check dimensions
    if (n != A.GetNumCols() || n != x.Size() || n != b.Size()) {
        throw std::invalid_argument("Matrix and vector dimensions don't match");
    }
    
    // Create temporary vectors
    std::unique_ptr<Vector> x_new = Vector::Create(n);
    std::unique_ptr<Vector> residual = Vector::Create(n);
    
    iteration_count_ = 0;
    converged_ = false;
    
    while (iteration_count_ < max_iterations_) {
        iteration_count_++;
        
        // Jacobi iteration: x_new[i] = (b[i] - sum_{j≠i} A[i,j]*x[j]) / A[i,i]
        for (Integer i = 0; i < n; ++i) {
            Real sum = 0.0;
            Real diag = A.GetElement(i, i);
            
            if (std::abs(diag) < 1e-20) {
                // Avoid division by zero
                (*x_new)[i] = x[i];
                continue;
            }
            
            for (Integer j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A.GetElement(i, j) * x[j];
                }
            }
            
            (*x_new)[i] = (b[i] - sum) / diag;
        }
        
        // Apply relaxation: x = (1 - ω)*x + ω*x_new
        for (Integer i = 0; i < n; ++i) {
            x[i] = (1.0 - relaxation_) * x[i] + relaxation_ * (*x_new)[i];
        }
        
        // Compute residual
        IterativeSolverUtils::ComputeResidual(A, x, b, *residual);
        residual_norm_ = IterativeSolverUtils::ComputeNorm(*residual);
        
        // Check convergence
        if (residual_norm_ < tolerance_) {
            converged_ = true;
            break;
        }
    }
    
    return converged_;
}

// ============================================================================
// Gauss-Seidel Solver Implementation
// ============================================================================

bool GaussSeidelSolver::Solve(const Matrix& A, Vector& x, const Vector& b) {
    Integer n = A.GetNumRows();
    
    // Check dimensions
    if (n != A.GetNumCols() || n != x.Size() || n != b.Size()) {
        throw std::invalid_argument("Matrix and vector dimensions don't match");
    }
    
    // Create temporary vector for residual
    std::unique_ptr<Vector> residual = Vector::Create(n);
    
    iteration_count_ = 0;
    converged_ = false;
    
    while (iteration_count_ < max_iterations_) {
        iteration_count_++;
        
        // Gauss-Seidel iteration
        for (Integer i = 0; i < n; ++i) {
            Real sum = 0.0;
            Real diag = A.GetElement(i, i);
            
            if (std::abs(diag) < 1e-20) {
                // Avoid division by zero
                continue;
            }
            
            for (Integer j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A.GetElement(i, j) * x[j];
                }
            }
            
            // Update x[i] immediately (uses latest values)
            Real x_new = (b[i] - sum) / diag;
            x[i] = (1.0 - relaxation_) * x[i] + relaxation_ * x_new;
        }
        
        // Compute residual
        IterativeSolverUtils::ComputeResidual(A, x, b, *residual);
        residual_norm_ = IterativeSolverUtils::ComputeNorm(*residual);
        
        // Check convergence
        if (residual_norm_ < tolerance_) {
            converged_ = true;
            break;
        }
    }
    
    return converged_;
}

// ============================================================================
// Iterative Solver Utility Functions
// ============================================================================

namespace IterativeSolverUtils {

void ComputeResidual(const Matrix& A, const Vector& x, const Vector& b, Vector& r) {
    Integer n = A.GetNumRows();
    
    if (n != A.GetNumCols() || n != x.Size() || n != b.Size() || n != r.Size()) {
        throw std::invalid_argument("Matrix and vector dimensions don't match");
    }
    
    // r = b - A*x
    A.Multiply(x, r);
    for (Integer i = 0; i < n; ++i) {
        r[i] = b[i] - r[i];
    }
}

Real ComputeNorm(const Vector& v) {
    Real sum = 0.0;
    Integer n = v.Size();
    
    for (Integer i = 0; i < n; ++i) {
        sum += v[i] * v[i];
    }
    
    return std::sqrt(sum);
}

Real DotProduct(const Vector& a, const Vector& b) {
    Integer n = a.Size();
    
    if (n != b.Size()) {
        throw std::invalid_argument("Vector dimensions don't match");
    }
    
    Real sum = 0.0;
    for (Integer i = 0; i < n; ++i) {
        sum += a[i] * b[i];
    }
    
    return sum;
}

void VectorAdd(Real alpha, const Vector& x, Vector& y) {
    Integer n = x.Size();
    
    if (n != y.Size()) {
        throw std::invalid_argument("Vector dimensions don't match");
    }
    
    for (Integer i = 0; i < n; ++i) {
        y[i] += alpha * x[i];
    }
}

void VectorAdd(Real alpha, const Vector& x, Real beta, Vector& y) {
    Integer n = x.Size();
    
    if (n != y.Size()) {
        throw std::invalid_argument("Vector dimensions don't match");
    }
    
    for (Integer i = 0; i < n; ++i) {
        y[i] = alpha * x[i] + beta * y[i];
    }
}

void IterativeSolverUtils::VectorCopy(const Vector& x, Vector& y) {
    Integer n = x.Size();
    
    if (n != y.Size()) {
        throw std::invalid_argument("Vector dimensions don't match");
    }
    
    for (Integer i = 0; i < n; ++i) {
        y[i] = x[i];
    }
}

bool IsSymmetricPositiveDefinite(const Matrix& A, Real tolerance) {
    Integer n = A.GetNumRows();
    
    if (n != A.GetNumCols()) {
        return false;
    }
    
    // Check symmetry
    for (Integer i = 0; i < n; ++i) {
        for (Integer j = i + 1; j < n; ++j) {
            Real a_ij = A.GetElement(i, j);
            Real a_ji = A.GetElement(j, i);
            
            if (std::abs(a_ij - a_ji) > tolerance) {
                return false;
            }
        }
    }
    
    // Check positive definiteness (simplified check)
    // For a proper check, we would need to compute eigenvalues
    // Here we just check that diagonal elements are positive
    for (Integer i = 0; i < n; ++i) {
        if (A.GetElement(i, i) <= 0.0) {
            return false;
        }
    }
    
    return true;
}

} // namespace IterativeSolverUtils

} // namespace elmer