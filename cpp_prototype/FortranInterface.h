// Fortran-C++互操作性接口
// 提供Elmer FEM Fortran代码与C++代码的桥梁

#ifndef FORTRAN_INTERFACE_H
#define FORTRAN_INTERFACE_H

#include <cstddef>

#ifdef __cplusplus
extern "C" {
#endif

// 基础数据类型映射
using FortranInteger = int;
using FortranReal = double;
using FortranLogical = int; // Fortran LOGICAL通常映射为int

// Fortran字符串结构（Fortran字符串有长度信息）
struct FortranString {
    char* data;
    std::size_t length;
};

// 矩阵操作接口（对应CRSMatrix.F90）
void crs_zero_matrix_fortran(void* matrix_ptr);
void crs_set_matrix_element_fortran(void* matrix_ptr, FortranInteger* i, 
                                   FortranInteger* j, FortranReal* value);
FortranReal crs_get_matrix_element_fortran(void* matrix_ptr, FortranInteger* i, 
                                          FortranInteger* j);
void crs_add_to_matrix_element_fortran(void* matrix_ptr, FortranInteger* i, 
                                      FortranInteger* j, FortranReal* value);

// 迭代求解器接口（对应IterSolve.F90）
void iter_solver_fortran(void* matrix_ptr, FortranReal* x, FortranReal* b, 
                        void* solver_ptr, FortranInteger* ndim);

// 矩阵组装接口（对应MatrixAssembly.F90）
void set_matrix_element_fortran(void* matrix_ptr, FortranInteger* i, 
                               FortranInteger* j, FortranReal* value);
void add_to_matrix_element_fortran(void* matrix_ptr, FortranInteger* i, 
                                  FortranInteger* j, FortranReal* value);

// 并行计算接口（对应ParallelUtils.F90）
void* parallel_init_fortran();
FortranInteger parallel_get_my_pe_fortran(void* parallel_env);
FortranInteger parallel_get_num_pe_fortran(void* parallel_env);

// 主求解器接口（对应ElmerSolver.F90）
void elmer_solver_fortran(FortranInteger* initialize);

// 工具函数接口（对应DefUtils.F90）
FortranReal get_control_value_fortran(FortranString* name, FortranInteger* len);

// 内存管理接口
void* allocate_fortran_memory(FortranInteger* size, FortranInteger* type_size);
void deallocate_fortran_memory(void* ptr);

#ifdef __cplusplus
}

// C++包装器类
namespace elmer {
namespace fortran {

class MatrixWrapper {
public:
    MatrixWrapper(void* fortran_matrix) : fortran_matrix_(fortran_matrix) {}
    
    void Zero() { crs_zero_matrix_fortran(fortran_matrix_); }
    void SetElement(int i, int j, double value) { 
        crs_set_matrix_element_fortran(fortran_matrix_, &i, &j, &value); 
    }
    double GetElement(int i, int j) const { 
        return crs_get_matrix_element_fortran(fortran_matrix_, &i, &j); 
    }
    void AddToElement(int i, int j, double value) { 
        crs_add_to_matrix_element_fortran(fortran_matrix_, &i, &j, &value); 
    }
    
private:
    void* fortran_matrix_;
};

class ParallelEnvironmentWrapper {
public:
    ParallelEnvironmentWrapper() {
        fortran_env_ = parallel_init_fortran();
    }
    
    int GetMyPE() const { 
        return parallel_get_my_pe_fortran(fortran_env_); 
    }
    int GetNumPE() const { 
        return parallel_get_num_pe_fortran(fortran_env_); 
    }
    
private:
    void* fortran_env_;
};

} // namespace fortran
} // namespace elmer

#endif // __cplusplus

#endif // FORTRAN_INTERFACE_H