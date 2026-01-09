// test_def_utils.cpp - Simple tests for DefUtils module

#include "DefUtils.h"
#include "CRSMatrix.h"
#include "ElmerCpp.h"
#include <iostream>
#include <cassert>

namespace elmer {

// Test setup function
void setup_test() {
    // Reset thread-local storage before each test
    DefUtils::GetIndexStore().clear();
    DefUtils::GetPermIndexStore().clear();
    DefUtils::GetValueStore(1).clear();
}

void test_VersionInformation() {
    std::string version = DefUtils::GetVersion();
    assert(!version.empty());
    
    bool found = false;
    std::string sif_name = DefUtils::GetSifName(&found);
    assert(!found); // Should not find SIF name without model context
    
    std::string revision = DefUtils::GetRevision(&found);
    // Revision may or may not be available depending on compilation
    
    std::string comp_date = DefUtils::GetCompilationDate(&found);
    // Compilation date may or may not be available
}

void test_MemoryStoreManagement() {
    // Test index store
    auto& index_store = DefUtils::GetIndexStore();
    assert(index_store.size() == DefUtils::ISTORE_MAX_SIZE);
    assert(index_store[0] == 0);
    
    // Test perm index store
    auto& perm_store = DefUtils::GetPermIndexStore();
    assert(perm_store.size() == DefUtils::ISTORE_MAX_SIZE);
    assert(perm_store[0] == 0);
    
    // Test value store
    Integer test_size = 100;
    auto& value_store = DefUtils::GetValueStore(test_size);
    assert(value_store.size() >= test_size);
    assert(value_store[0] == 0.0);
    
    // Test store size validation
    try {
        DefUtils::GetValueStore(DefUtils::VSTORE_MAX_SIZE + 1);
        assert(false && "Expected exception not thrown");
    } catch (const std::runtime_error&) {
        // Expected
    } catch (...) {
        assert(false && "Wrong exception type");
    }
}

void test_StringUtilities() {
    bool found = false;
    std::string result = DefUtils::GetString("nonexistent_key", &found);
    assert(!found);
    assert(result.empty());
    
    // Test with default values
    Real real_val = DefUtils::GetConstReal("nonexistent", 42.0);
    assert(real_val == 42.0);
    
    Integer int_val = DefUtils::GetConstInteger("nonexistent", 123);
    assert(int_val == 123);
    
    bool bool_val = DefUtils::GetLogical("nonexistent", true);
    assert(bool_val);
}

void test_ErrorHandling() {
    // Test info message (should not throw)
    try {
        DefUtils::Info("TestSolver", "Test message", 1);
    } catch (...) {
        assert(false && "Unexpected exception thrown");
    }
    
    // Test warning message (should not throw)
    try {
        DefUtils::Warning("TestSolver", "Test warning");
    } catch (...) {
        assert(false && "Unexpected exception thrown");
    }
    
    // Test fatal error (should throw)
    try {
        DefUtils::Fatal("TestRoutine", "Test fatal error");
        assert(false && "Expected exception not thrown");
    } catch (const std::runtime_error&) {
        // Expected
    } catch (...) {
        assert(false && "Wrong exception type");
    }
}

void test_TimeIntegration() {
    std::vector<Real> old_values = {1.0, 2.0, 3.0};
    std::vector<Real> new_values = {4.0, 5.0, 6.0};
    std::vector<Real> old_derivatives = {0.1, 0.2, 0.3};
    std::vector<Real> new_derivatives(3, 0.0);
    
    Real dt = 0.1;
    Real alpha = 0.5;
    Real beta = 0.25;
    
    // Test first order time integration
    auto test_new_values = new_values; // Copy for testing
    try {
        DefUtils::Default1stOrderTime(dt, alpha, old_values, test_new_values);
    } catch (...) {
        assert(false && "Unexpected exception thrown");
    }
    
    // Verify size remains the same
    assert(test_new_values.size() == old_values.size());
    
    // Test second order time integration
    try {
        DefUtils::Default2ndOrderTime(dt, alpha, beta, 
                                      old_values, old_derivatives,
                                      new_values, new_derivatives);
    } catch (...) {
        assert(false && "Unexpected exception thrown");
    }
    
    // Test size mismatch errors
    std::vector<Real> wrong_size = {1.0, 2.0};
    try {
        DefUtils::Default1stOrderTime(dt, alpha, wrong_size, new_values);
        assert(false && "Expected exception not thrown");
    } catch (const std::invalid_argument&) {
        // Expected
    } catch (...) {
        assert(false && "Wrong exception type");
    }
    
    try {
        DefUtils::Default2ndOrderTime(dt, alpha, beta, 
                                      wrong_size, old_derivatives,
                                      new_values, new_derivatives);
        assert(false && "Expected exception not thrown");
    } catch (const std::invalid_argument&) {
        // Expected
    } catch (...) {
        assert(false && "Wrong exception type");
    }
}

void test_ThreadLocalStorage() {
    // Test that each thread gets its own storage
    auto& store1 = DefUtils::GetIndexStore();
    store1[0] = 42; // Modify in this thread
    
    // The modification should persist within the same thread
    auto& store2 = DefUtils::GetIndexStore();
    assert(store2[0] == 42);
    
    // Test value store with different sizes
    auto& small_store = DefUtils::GetValueStore(10);
    assert(small_store.size() >= 10);
    
    auto& large_store = DefUtils::GetValueStore(100);
    assert(large_store.size() >= 100);
}

} // namespace elmer

// Test runner
int main() {
    std::cout << "Running DefUtils tests..." << std::endl;
    
    try {
        elmer::setup_test();
        elmer::test_VersionInformation();
        std::cout << "✓ VersionInformation test passed" << std::endl;
        
        elmer::setup_test();
        elmer::test_MemoryStoreManagement();
        std::cout << "✓ MemoryStoreManagement test passed" << std::endl;
        
        elmer::setup_test();
        elmer::test_StringUtilities();
        std::cout << "✓ StringUtilities test passed" << std::endl;
        
        elmer::setup_test();
        elmer::test_ErrorHandling();
        std::cout << "✓ ErrorHandling test passed" << std::endl;
        
        elmer::setup_test();
        elmer::test_TimeIntegration();
        std::cout << "✓ TimeIntegration test passed" << std::endl;
        
        elmer::setup_test();
        elmer::test_ThreadLocalStorage();
        std::cout << "✓ ThreadLocalStorage test passed" << std::endl;
        
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown test failure" << std::endl;
        return 1;
    }
}