// test_linear_algebra.cpp - Test for linear algebra modules

#include "CRSMatrix.h"
#include "IterativeSolver.h"
#include <iostream>
#include <cmath>
#include <memory>

using namespace elmer;

// Test function for CRS matrix operations
void testCRSMatrix() {
    std::cout << "=== Testing CRS Matrix Operations ===" << std::endl;
    
    // Create a 3x3 matrix
    auto matrix = Matrix::CreateCRS(3, 3);
    
    // Set some values (creating a symmetric positive definite matrix)
    matrix->SetElement(0, 0, 4.0);
    matrix->SetElement(0, 1, 1.0);
    matrix->SetElement(0, 2, 0.0);
    
    matrix->SetElement(1, 0, 1.0);
    matrix->SetElement(1, 1, 3.0);
    matrix->SetElement(1, 2, 1.0);
    
    matrix->SetElement(2, 0, 0.0);
    matrix->SetElement(2, 1, 1.0);
    matrix->SetElement(2, 2, 2.0);
    
    // Test matrix-vector multiplication
    auto x = Vector::Create(3);
    auto y = Vector::Create(3);
    
    (*x)[0] = 1.0;
    (*x)[1] = 2.0;
    (*x)[2] = 3.0;
    
    matrix->Multiply(*x, *y);
    
    std::cout << "Matrix-vector multiplication: A*x = [" 
              << (*y)[0] << ", " << (*y)[1] << ", " << (*y)[2] << "]" << std::endl;
    
    // Expected result: [4*1 + 1*2 + 0*3, 1*1 + 3*2 + 1*3, 0*1 + 1*2 + 2*3] = [6, 10, 8]
    
    // Test symmetry
    bool symmetric = IterativeSolverUtils::IsSymmetricPositiveDefinite(*matrix);
    std::cout << "Matrix is symmetric positive definite: " << (symmetric ? "Yes" : "No") << std::endl;
    
    std::cout << "CRS Matrix test completed successfully!" << std::endl;
}

// Test function for Conjugate Gradient solver
void testConjugateGradient() {
    std::cout << "\n=== Testing Conjugate Gradient Solver ===" << std::endl;
    
    // Create a symmetric positive definite matrix
    auto matrix = Matrix::CreateCRS(3, 3);
    
    matrix->SetElement(0, 0, 4.0);
    matrix->SetElement(0, 1, 1.0);
    matrix->SetElement(0, 2, 0.0);
    
    matrix->SetElement(1, 0, 1.0);
    matrix->SetElement(1, 1, 3.0);
    matrix->SetElement(1, 2, 1.0);
    
    matrix->SetElement(2, 0, 0.0);
    matrix->SetElement(2, 1, 1.0);
    matrix->SetElement(2, 2, 2.0);
    
    // Create right-hand side vector
    auto b = Vector::Create(3);
    (*b)[0] = 6.0;
    (*b)[1] = 10.0;
    (*b)[2] = 8.0;
    
    // Solution should be [1, 2, 3]
    auto x = Vector::Create(3);
    x->Zero();  // Start with zero initial guess
    
    // Create and run CG solver
    ConjugateGradientSolver cg_solver(100, 1e-10);
    bool success = cg_solver.Solve(*matrix, *x, *b);
    
    std::cout << "CG solver converged: " << (success ? "Yes" : "No") << std::endl;
    std::cout << "Iterations: " << cg_solver.GetIterationCount() << std::endl;
    std::cout << "Final residual: " << cg_solver.GetResidualNorm() << std::endl;
    std::cout << "Solution: [" << (*x)[0] << ", " << (*x)[1] << ", " << (*x)[2] << "]" << std::endl;
    
    // Check solution accuracy
    auto residual = Vector::Create(3);
    IterativeSolverUtils::ComputeResidual(*matrix, *x, *b, *residual);
    Real residual_norm = IterativeSolverUtils::ComputeNorm(*residual);
    
    std::cout << "Residual norm: " << residual_norm << std::endl;
    std::cout << "Conjugate Gradient test completed successfully!" << std::endl;
}

// Test function for Jacobi solver
void testJacobiSolver() {
    std::cout << "\n=== Testing Jacobi Solver ===" << std::endl;
    
    // Create a diagonally dominant matrix
    auto matrix = Matrix::CreateCRS(3, 3);
    
    matrix->SetElement(0, 0, 4.0);
    matrix->SetElement(0, 1, 1.0);
    matrix->SetElement(0, 2, 0.0);
    
    matrix->SetElement(1, 0, 1.0);
    matrix->SetElement(1, 1, 5.0);
    matrix->SetElement(1, 2, 1.0);
    
    matrix->SetElement(2, 0, 0.0);
    matrix->SetElement(2, 1, 1.0);
    matrix->SetElement(2, 2, 6.0);
    
    // Create right-hand side vector
    auto b = Vector::Create(3);
    (*b)[0] = 1.0;
    (*b)[1] = 2.0;
    (*b)[2] = 3.0;
    
    auto x = Vector::Create(3);
    x->Zero();  // Start with zero initial guess
    
    // Create and run Jacobi solver
    JacobiSolver jacobi_solver(1000, 1e-8, 1.0);
    bool success = jacobi_solver.Solve(*matrix, *x, *b);
    
    std::cout << "Jacobi solver converged: " << (success ? "Yes" : "No") << std::endl;
    std::cout << "Iterations: " << jacobi_solver.GetIterationCount() << std::endl;
    std::cout << "Final residual: " << jacobi_solver.GetResidualNorm() << std::endl;
    std::cout << "Solution: [" << (*x)[0] << ", " << (*x)[1] << ", " << (*x)[2] << "]" << std::endl;
    
    std::cout << "Jacobi solver test completed successfully!" << std::endl;
}

// Test function for Gauss-Seidel solver
void testGaussSeidelSolver() {
    std::cout << "\n=== Testing Gauss-Seidel Solver ===" << std::endl;
    
    // Create a diagonally dominant matrix
    auto matrix = Matrix::CreateCRS(3, 3);
    
    matrix->SetElement(0, 0, 4.0);
    matrix->SetElement(0, 1, 1.0);
    matrix->SetElement(0, 2, 0.0);
    
    matrix->SetElement(1, 0, 1.0);
    matrix->SetElement(1, 1, 5.0);
    matrix->SetElement(1, 2, 1.0);
    
    matrix->SetElement(2, 0, 0.0);
    matrix->SetElement(2, 1, 1.0);
    matrix->SetElement(2, 2, 6.0);
    
    // Create right-hand side vector
    auto b = Vector::Create(3);
    (*b)[0] = 1.0;
    (*b)[1] = 2.0;
    (*b)[2] = 3.0;
    
    auto x = Vector::Create(3);
    x->Zero();  // Start with zero initial guess
    
    // Create and run Gauss-Seidel solver
    GaussSeidelSolver gs_solver(1000, 1e-8, 1.0);
    bool success = gs_solver.Solve(*matrix, *x, *b);
    
    std::cout << "Gauss-Seidel solver converged: " << (success ? "Yes" : "No") << std::endl;
    std::cout << "Iterations: " << gs_solver.GetIterationCount() << std::endl;
    std::cout << "Final residual: " << gs_solver.GetResidualNorm() << std::endl;
    std::cout << "Solution: [" << (*x)[0] << ", " << (*x)[1] << ", " << (*x)[2] << "]" << std::endl;
    
    std::cout << "Gauss-Seidel solver test completed successfully!" << std::endl;
}

int main() {
    try {
        testCRSMatrix();
        testConjugateGradient();
        testJacobiSolver();
        testGaussSeidelSolver();
        
        std::cout << "\n=== All linear algebra tests passed! ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}