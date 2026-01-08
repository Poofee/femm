// test_matrix_simple.cpp - Simple matrix compilation test
// Tests basic matrix functionality without complex interfaces

#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>

// Simple matrix interface for testing
class SimpleMatrix {
public:
    virtual ~SimpleMatrix() = default;
    virtual void zero() = 0;
    virtual void setElement(int i, int j, double value) = 0;
    virtual double getElement(int i, int j) const = 0;
    virtual void addToElement(int i, int j, double value) = 0;
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
    std::vector<double> multiply(const std::vector<double>& x) const {
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

int main() {
    std::cout << "Testing CRS Matrix Implementation" << std::endl;
    std::cout << "=================================" << std::endl;
    
    // Create a 3x3 matrix
    CRSMatrix matrix(3, 3);
    
    // Set matrix elements (symmetric positive definite matrix)
    matrix.setElement(0, 0, 4.0);
    matrix.setElement(0, 1, -1.0);
    matrix.setElement(0, 2, 0.0);
    
    matrix.setElement(1, 0, -1.0);
    matrix.setElement(1, 1, 4.0);
    matrix.setElement(1, 2, -1.0);
    
    matrix.setElement(2, 0, 0.0);
    matrix.setElement(2, 1, -1.0);
    matrix.setElement(2, 2, 4.0);
    
    // Verify matrix elements
    std::cout << "Matrix element verification:" << std::endl;
    std::cout << "A[0,0] = " << matrix.getElement(0, 0) << " (expected: 4.0)" << std::endl;
    std::cout << "A[1,1] = " << matrix.getElement(1, 1) << " (expected: 4.0)" << std::endl;
    std::cout << "A[2,2] = " << matrix.getElement(2, 2) << " (expected: 4.0)" << std::endl;
    
    // Test matrix-vector multiplication
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> b = matrix.multiply(x);
    
    std::cout << "Matrix-vector multiplication test:" << std::endl;
    std::cout << "x = [" << x[0] << ", " << x[1] << ", " << x[2] << "]" << std::endl;
    std::cout << "A*x = [" << b[0] << ", " << b[1] << ", " << b[2] << "]" << std::endl;
    
    // Verify results
    double expected_b0 = 4.0*1.0 - 1.0*2.0 + 0.0*3.0;
    double expected_b1 = -1.0*1.0 + 4.0*2.0 - 1.0*3.0;
    double expected_b2 = 0.0*1.0 - 1.0*2.0 + 4.0*3.0;
    
    std::cout << "Expected result: [" << expected_b0 << ", " << expected_b1 << ", " << expected_b2 << "]" << std::endl;
    
    bool test_passed = std::abs(b[0] - expected_b0) < 1e-10 &&
                      std::abs(b[1] - expected_b1) < 1e-10 &&
                      std::abs(b[2] - expected_b2) < 1e-10;
    
    std::cout << "Test " << (test_passed ? "PASSED" : "FAILED") << std::endl;
    
    return test_passed ? 0 : 1;
}