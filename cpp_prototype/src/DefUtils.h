// DefUtils.h - Elmer FEM C++ Default Utilities Module
// Corresponds to Fortran module: DefUtils.F90

#pragma once

#include "ElmerCpp.h"
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <functional>
#include <stdexcept>
#include <iostream>

namespace elmer {

// Forward declarations - these classes will be implemented later
// class Model;
// class Solver;
// class Element;
// class Variable;

/**
 * @brief Default utilities module for Elmer FEM
 * 
 * This class provides default implementations for various utility functions
 * used throughout Elmer FEM. It corresponds to the functionality in DefUtils.F90.
 */
class DefUtils {
public:
    // Memory management constants
    static constexpr Integer ISTORE_MAX_SIZE = 1024;
    static constexpr Integer VSTORE_MAX_SIZE = 1024;

    // Version information
    static std::string GetVersion();
    static std::string GetSifName(bool* found = nullptr);
    static std::string GetRevision(bool* found = nullptr);
    static std::string GetCompilationDate(bool* found = nullptr);

    // Memory store management (thread-local)
    static std::vector<Integer>& GetIndexStore();
    static std::vector<Integer>& GetPermIndexStore();
    static std::vector<Real>& GetValueStore(Integer n);

    // Solver and model access (placeholder implementations)
    // static std::shared_ptr<Solver> GetSolver();
    // static std::shared_ptr<Model> GetCurrentModel();
    // static std::shared_ptr<Element> GetCurrentElement();

    // String utilities
    static std::string GetString(const std::string& key, bool* found = nullptr);
    // static std::string GetString(const std::shared_ptr<Model>& model, 
    //                             const std::string& key, bool* found = nullptr);
    
    // Numerical utilities
    static Real GetConstReal(const std::string& key, Real default_value = 0.0);
    static Integer GetConstInteger(const std::string& key, Integer default_value = 0);
    static bool GetLogical(const std::string& key, bool default_value = false);

    // Error handling and logging
    static void Info(const std::string& solver_name, const std::string& message, 
                    Integer level = 1);
    static void Warning(const std::string& solver_name, const std::string& message);
    static void Fatal(const std::string& routine_name, const std::string& message);

    // Element utilities (placeholder implementations)
    // static std::vector<Real> GetLocalSolution(const std::shared_ptr<Element>& element,
    //                                         const std::shared_ptr<Variable>& variable);
    // static std::vector<Real> GetVectorLocalSolution(const std::shared_ptr<Element>& element,
    //                                               const std::shared_ptr<Variable>& variable);

    // Time integration utilities
    static void Default1stOrderTime(Real dt, Real alpha, 
                                   const std::vector<Real>& old_values,
                                   std::vector<Real>& new_values);
    static void Default2ndOrderTime(Real dt, Real alpha, Real beta,
                                   const std::vector<Real>& old_values,
                                   const std::vector<Real>& old_derivatives,
                                   std::vector<Real>& new_values,
                                   std::vector<Real>& new_derivatives);

    // Matrix update utilities (placeholder implementations)
    // static void DefaultUpdateEquations(const std::shared_ptr<Matrix>& matrix,
    //                                   const std::vector<Real>& local_matrix,
    //                                   const std::shared_ptr<Element>& element);
    // static void DefaultUpdateForce(const std::shared_ptr<Vector>& force_vector,
    //                               const std::vector<Real>& local_force,
    //                               const std::shared_ptr<Element>& element);

private:
    // Thread-local storage
    static thread_local std::vector<Integer> index_store_;
    static thread_local std::vector<Integer> perm_index_store_;
    static thread_local std::vector<Real> value_store_;
    
    // Current context (thread-local)
    // static thread_local std::shared_ptr<Element> current_element_;
    
    // Initialize thread-local storage
    static void InitializeThreadLocalStorage();
    
    // Internal helper functions
    static void CheckStoreSize(Integer required_size, Integer max_size, 
                              const std::string& routine_name);
};

} // namespace elmer