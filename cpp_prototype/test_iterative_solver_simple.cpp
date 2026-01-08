// test_iterative_solver_simple.cpp - Simple iterative solver test
// Tests basic iterative solver functionality

#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <cmath>

// Simple matrix interface for testing
class SimpleMatrix {
public:
    virtual ~SimpleMatrix() = default;
    virtual void zero() = 0;
    virtual void setElement(int i, int j, double value) = 0;
    virtual double getElement(int i, int j) const = 0;
    virtual void addToElement(int i, int j, double value) = 0;
    virtual std::vector<double> multiply(const std::vector<double>& x) const = 0;
};

// CRS matrix implementation
class CRSMatrix : public SimpleMatrix {
public:
    CRSMatrix(int nrows, int ncols) 
        : nrows_(nrows), ncols_(ncols) {
        rows_.resize(nrows + 1, 0);
    }
    
    ~CRSMatrix() override = default;
    
    void zero() override {
        std::fill(values_.begin(), values_.end(), 0.0);
    }
    
    void setElement(int i, int j, double value) override {
        // Find or insert element
        auto pos = findElementPosition(i, j);
        if (pos != values_.end()) {
            *pos = value;
        } else {
            insertElement(i, j, value);
        }
    }
    
    double getElement(int i, int j) const override {
        auto pos = findElementPosition(i, j);
        return (pos != values_.end()) ? *pos : 0.0;
    }
    
    void addToElement(int i, int j, double value) override {
        auto pos = findElementPosition(i, j);
        if (pos != values_.end()) {
            *pos += value;
        } else {
            insertElement(i, j, value);
        }
    }
    
    // Matrix-vector multiplication
    std::vector<double> multiply(const std::vector<double>& x) const override {
        std::vector<double> result(nrows_, 0.0);
        for (int i = 0; i < nrows_; ++i) {
            double sum = 0.0;
            for (int k = rows_[i]; k < rows_[i + 1]; ++k) {
                sum += values_[k] * x[cols_[k]];
            }
            result[i] = sum;
        }
        return result;
    }
    
private:
    int nrows_, ncols_;
    std::vector<int> rows_;  // Row pointers
    std::vector<int> cols_;  // Column indices
    std::vector<double> values_;   // Non-zero values
    
    // Find element position
    std::vector<double>::iterator findElementPosition(int i, int j) {
        int start = rows_[i];
        int end = rows_[i + 1];
        
        for (int k = start; k < end; ++k) {
            if (cols_[k] == j) {
                return values_.begin() + k;
            }
        }
        return values_.end();
    }
    
    std::vector<double>::const_iterator findElementPosition(int i, int j) const {
        int start = rows_[i];
        int end = rows_[i + 1];
        
        for (int k = start; k < end; ++k) {
            if (cols_[k] == j) {
                return values_.begin() + k;
            }
        }
        return values_.end();
    }
    
    // Insert new element
    void insertElement(int i, int j, double value) {
        int start = rows_[i];
        int end = rows_[i + 1];
        
        // Find insertion position
        int insert_pos = start;
        while (insert_pos < end && cols_[insert_pos] < j) {
            ++insert_pos;
        }
        
        // Insert new element
        cols_.insert(cols_.begin() + insert_pos, j);
        values_.insert(values_.begin() + insert_pos, value);
        
        // Update row pointers
        for (int k = i + 1; k <= nrows_; ++k) {
            rows_[k]++;
        }
    }
};

// Simple iterative solver interface
class IterativeSolver {
public:
    virtual ~IterativeSolver() = default;
    virtual bool solve(const SimpleMatrix& A, std::vector<double>& x, const std::vector<double>& b) = 0;
    
    void setMaxIterations(int max_iter) { max_iterations_ = max_iter; }
    void setTolerance(double tol) { tolerance_ = tol; }
    
    int getIterations() const { return iterations_; }
    double getResidual() const { return residual_; }
    
protected:
    int max_iterations_ = 1000;
    double tolerance_ = 1e-8;
    int iterations_ = 0;
    double residual_ = 0.0;
    
    // Compute vector norm
    double norm(const std::vector<double>& v) const {
        double sum = 0.0;
        for (double val : v) {
            sum += val * val;
        }
        return std::sqrt(sum);
    }
    
    // Compute residual
    double computeResidual(const SimpleMatrix& A, const std::vector<double>& x, const std::vector<double>& b) {
        std::vector<double> r = A.multiply(x);
        
        for (size_t i = 0; i < b.size(); ++i) {
            r[i] = b[i] - r[i];
        }
        
        return norm(r);
    }
};

// Conjugate Gradient solver
class CGSolver : public IterativeSolver {
public:
    bool solve(const SimpleMatrix& A, std::vector<double>& x, const std::vector<double>& b) override {
        int n = static_cast<int>(b.size());
        
        // Initialize vectors
        std::vector<double> r(n), p(n), Ap(n);
        
        // Compute initial residual
        Ap = A.multiply(x);
        for (int i = 0; i < n; ++i) {
            r[i] = b[i] - Ap[i];
        }
        
        p = r;  // p0 = r0
        
        double r_norm_sq = 0.0;
        for (double val : r) {
            r_norm_sq += val * val;
        }
        
        iterations_ = 0;
        
        while (iterations_ < max_iterations_) {
            // Compute A*p
            Ap = A.multiply(p);
            
            // Compute alpha
            double pAp = 0.0;
            for (int i = 0; i < n; ++i) {
                pAp += p[i] * Ap[i];
            }
            
            if (std::abs(pAp) < 1e-15) {
                break;  // Avoid division by zero
            }
            
            double alpha = r_norm_sq / pAp;
            
            // Update solution and residual
            for (int i = 0; i < n; ++i) {
                x[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
            }
            
            // Check convergence
            double new_r_norm_sq = 0.0;
            for (double val : r) {
                new_r_norm_sq += val * val;
            }
            
            residual_ = std::sqrt(new_r_norm_sq);
            
            if (residual_ < tolerance_) {
                break;
            }
            
            // Compute beta and update p
            double beta = new_r_norm_sq / r_norm_sq;
            
            for (int i = 0; i < n; ++i) {
                p[i] = r[i] + beta * p[i];
            }
            
            r_norm_sq = new_r_norm_sq;
            iterations_++;
        }
        
        return residual_ < tolerance_;
    }
};

int main() {
    std::cout << "Testing Iterative Solver Implementation" << std::endl;
    std::cout << "=======================================" << std::endl;
    
    // Create a 3x3 symmetric positive definite matrix
    CRSMatrix A(3, 3);
    
    // Set matrix elements
    A.setElement(0, 0, 4.0);
    A.setElement(0, 1, -1.0);
    A.setElement(0, 2, 0.0);
    
    A.setElement(1, 0, -1.0);
    A.setElement(1, 1, 4.0);
    A.setElement(1, 2, -1.0);
    
    A.setElement(2, 0, 0.0);
    A.setElement(2, 1, -1.0);
    A.setElement(2, 2, 4.0);
    
    // Create right-hand side vector
    std::vector<double> b = {1.0, 2.0, 3.0};
    
    // Initial guess (zero vector)
    std::vector<double> x(3, 0.0);
    
    // Create and configure CG solver
    CGSolver solver;
    solver.setMaxIterations(100);
    solver.setTolerance(1e-10);
    
    std::cout << "Solving linear system Ax = b using Conjugate Gradient method" << std::endl;
    std::cout << "Matrix A:" << std::endl;
    std::cout << "[4, -1,  0]" << std::endl;
    std::cout << "[-1, 4, -1]" << std::endl;
    std::cout << "[0, -1,  4]" << std::endl;
    std::cout << "Right-hand side b = [" << b[0] << ", " << b[1] << ", " << b[2] << "]" << std::endl;
    
    // Solve the linear system
    bool converged = solver.solve(A, x, b);
    
    std::cout << "\nSolution:" << std::endl;
    std::cout << "x = [" << x[0] << ", " << x[1] << ", " << x[2] << "]" << std::endl;
    
    std::cout << "\nSolver statistics:" << std::endl;
    std::cout << "Iterations: " << solver.getIterations() << std::endl;
    std::cout << "Final residual: " << solver.getResidual() << std::endl;
    std::cout << "Converged: " << (converged ? "YES" : "NO") << std::endl;
    
    // Verify solution by computing A*x
    std::vector<double> Ax = A.multiply(x);
    std::cout << "\nVerification:" << std::endl;
    std::cout << "A*x = [" << Ax[0] << ", " << Ax[1] << ", " << Ax[2] << "]" << std::endl;
    std::cout << "b = [" << b[0] << ", " << b[1] << ", " << b[2] << "]" << std::endl;
    
    // Check accuracy
    double error = 0.0;
    for (size_t i = 0; i < b.size(); ++i) {
        error += std::abs(Ax[i] - b[i]);
    }
    
    std::cout << "Absolute error: " << error << std::endl;
    std::cout << "Test " << (error < 1e-10 ? "PASSED" : "FAILED") << std::endl;
    
    return error < 1e-10 ? 0 : 1;
}