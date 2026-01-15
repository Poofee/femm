/**
 * @file ElementUtils.h
 * @brief 有限元元素工具模块
 * 
 * 对应Fortran模块：ElementUtils.F90
 * 提供有限元计算中的元素级别实用函数，包括矩阵创建、初始化、释放等操作。
 */

#pragma once

#include "Types.h"
#include <memory>
#include <vector>
#include <stdexcept>

namespace elmer {

/**
 * @brief 元素工具类
 * 
 * 提供有限元计算中的元素级别实用函数，包括矩阵操作、元素索引计算等。
 */
class ElementUtils {
public:
    /**
     * @brief 矩阵格式枚举
     */
    enum class MatrixFormat {
        CRS,        ///< 压缩行存储格式
        LIST,       ///< 列表存储格式  
        BAND,       ///< 带状矩阵格式
        SBAND       ///< 对称带状矩阵格式
    };

    /**
     * @brief 释放矩阵结构
     * @param matrix 要释放的矩阵指针
     * 
     * 对应Fortran子程序：FreeMatrix
     */
    static void freeMatrix(std::shared_ptr<Matrix> matrix);

    /**
     * @brief 创建列表矩阵
     * @param model 模型对象
     * @param solver 求解器对象
     * @param mesh 网格对象
     * @param list 列表结构
     * @param reorder 是否重新排序
     * @param dofs 自由度数量
     * @param matrixFormat 矩阵格式
     * @return 创建的矩阵指针
     * 
     * 对应Fortran子程序：MakeListMatrix
     */
    static std::shared_ptr<Matrix> makeListMatrix(
        std::shared_ptr<void> model,
        std::shared_ptr<void> solver,
        std::shared_ptr<void> mesh,
        std::shared_ptr<void> list,
        bool reorder,
        int dofs,
        MatrixFormat matrixFormat);

    /**
     * @brief 创建列表矩阵数组
     * @param model 模型对象
     * @param solver 求解器对象
     * @param mesh 网格对象
     * @param list 列表结构
     * @param reorder 是否重新排序
     * @param dofs 自由度数量
     * @param matrixFormat 矩阵格式
     * @return 创建的矩阵指针数组
     * 
     * 对应Fortran子程序：MakeListMatrixArray
     */
    static std::vector<std::shared_ptr<Matrix>> makeListMatrixArray(
        std::shared_ptr<void> model,
        std::shared_ptr<void> solver,
        std::shared_ptr<void> mesh,
        std::shared_ptr<void> list,
        bool reorder,
        int dofs,
        MatrixFormat matrixFormat);

    /**
     * @brief 初始化矩阵
     * @param matrix 矩阵指针
     * @param n 矩阵大小
     * @param list 列表结构
     * @param dofs 自由度数量
     * @param reorder 是否重新排序
     * @param invInitialReorder 初始逆序
     * 
     * 对应Fortran子程序：InitializeMatrix
     */
    static void initializeMatrix(
        std::shared_ptr<Matrix> matrix,
        int n,
        std::shared_ptr<void> list,
        int dofs,
        bool reorder,
        std::shared_ptr<void> invInitialReorder);

    /**
     * @brief 创建矩阵
     * @param model 模型对象
     * @param solver 求解器对象
     * @param mesh 网格对象
     * @param perm 排列数组
     * @param dofs 自由度数量
     * @param matrixFormat 矩阵格式
     * @param reorder 是否重新排序
     * @return 创建的矩阵指针
     * 
     * 对应Fortran子程序：CreateMatrix
     */
    static std::shared_ptr<Matrix> createMatrix(
        std::shared_ptr<void> model,
        std::shared_ptr<void> solver,
        std::shared_ptr<void> mesh,
        const std::vector<int>& perm,
        int dofs,
        MatrixFormat matrixFormat,
        bool reorder);

    /**
     * @brief 计算DG辐射索引
     * @param element 元素对象
     * @param n 元素数量
     * @param elemInds 元素索引数组
     * @param diffuseGray 是否使用漫射灰度
     * 
     * 对应Fortran子程序：DgRadiationIndexes
     */
    static void dgRadiationIndexes(
        std::shared_ptr<void> element,
        int n,
        std::vector<int>& elemInds,
        bool diffuseGray = false);

    /**
     * @brief 打包排序添加操作
     * @param n 元素数量
     * @param ind 索引数组
     * @param perm 排列数组
     * 
     * 对应Fortran子程序：PackSortAdd
     */
    static void packSortAdd(int n, std::vector<int>& ind, std::vector<int>& perm);

    /**
     * @brief 检查矩阵是否为空
     * @param matrix 矩阵指针
     * @return 如果矩阵为空返回true，否则返回false
     */
    static bool isMatrixEmpty(std::shared_ptr<Matrix> matrix);

    /**
     * @brief 获取矩阵大小
     * @param matrix 矩阵指针
     * @return 矩阵大小
     */
    static int getMatrixSize(std::shared_ptr<Matrix> matrix);

    /**
     * @brief 验证矩阵格式
     * @param format 矩阵格式
     * @return 如果格式有效返回true，否则返回false
     */
    static bool validateMatrixFormat(MatrixFormat format);

private:
    /**
     * @brief 内部初始化矩阵实现
     */
    static void initializeMatrixInternal(std::shared_ptr<Matrix> matrix, int n);

    /**
     * @brief 创建CRS矩阵
     */
    static std::shared_ptr<Matrix> createCRSMatrix(int n, int dofs);

    /**
     * @brief 创建带状矩阵
     */
    static std::shared_ptr<Matrix> createBandMatrix(int n, int dofs);

    /**
     * @brief 创建对称带状矩阵
     */
    static std::shared_ptr<Matrix> createSBandMatrix(int n, int dofs);
};

} // namespace elmer