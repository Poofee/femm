/**
 * @file test_wedge_functions.cpp
 * @brief 棱柱单元形函数测试程序
 * 
 * 测试棱柱单元P型和L型基函数及其导数的正确性。
 */

#include "src/ElementDescription.h"
#include <iostream>
#include <vector>
#include <iomanip>

using namespace ElmerCpp;

/**
 * @brief 测试棱柱单元基函数
 */
void testWedgeBasisFunctions() {
    std::cout << "=== 棱柱单元基函数测试 ===" << std::endl;
    
    // 测试点：棱柱单元中心
    double u = 0.0, v = 1.0/3.0, w = 0.0;
    
    std::vector<double> phi_p(6), phi_l(6);
    std::vector<std::vector<double>> gradphi_p(6, std::vector<double>(3));
    std::vector<std::vector<double>> gradphi_l(6, std::vector<double>(3));
    
    // 测试P型基函数
    WedgeNodalPBasisAll(u, v, w, phi_p);
    std::cout << "P型基函数在(" << u << ", " << v << ", " << w << ")处的值：" << std::endl;
    for (int i = 0; i < 6; ++i) {
        std::cout << "  phi_p[" << i << "] = " << std::setprecision(10) << phi_p[i] << std::endl;
    }
    
    // 测试L型基函数
    WedgeNodalLBasisAll(u, v, w, phi_l);
    std::cout << "\nL型基函数在(" << u << ", " << v << ", " << w << ")处的值：" << std::endl;
    for (int i = 0; i < 6; ++i) {
        std::cout << "  phi_l[" << i << "] = " << std::setprecision(10) << phi_l[i] << std::endl;
    }
    
    // 测试P型基函数导数
    dWedgeNodalPBasisAll(u, v, w, gradphi_p);
    std::cout << "\nP型基函数导数：" << std::endl;
    for (int i = 0; i < 6; ++i) {
        std::cout << "  节点" << i << ": du=" << gradphi_p[i][0] 
                  << ", dv=" << gradphi_p[i][1] 
                  << ", dw=" << gradphi_p[i][2] << std::endl;
    }
    
    // 测试L型基函数导数
    dWedgeNodalLBasisAll(u, v, w, gradphi_l);
    std::cout << "\nL型基函数导数：" << std::endl;
    for (int i = 0; i < 6; ++i) {
        std::cout << "  节点" << i << ": du=" << gradphi_l[i][0] 
                  << ", dv=" << gradphi_l[i][1] 
                  << ", dw=" << gradphi_l[i][2] << std::endl;
    }
    
    // 验证基函数性质：所有基函数之和应为1
    double sum_p = 0.0, sum_l = 0.0;
    for (int i = 0; i < 6; ++i) {
        sum_p += phi_p[i];
        sum_l += phi_l[i];
    }
    
    std::cout << "\n基函数性质验证：" << std::endl;
    std::cout << "P型基函数和 = " << sum_p << " (应为1)" << std::endl;
    std::cout << "L型基函数和 = " << sum_l << " (应为1)" << std::endl;
    
    // 验证导数性质：所有基函数导数之和应为0
    double sum_du_p = 0.0, sum_dv_p = 0.0, sum_dw_p = 0.0;
    double sum_du_l = 0.0, sum_dv_l = 0.0, sum_dw_l = 0.0;
    
    for (int i = 0; i < 6; ++i) {
        sum_du_p += gradphi_p[i][0];
        sum_dv_p += gradphi_p[i][1];
        sum_dw_p += gradphi_p[i][2];
        
        sum_du_l += gradphi_l[i][0];
        sum_dv_l += gradphi_l[i][1];
        sum_dw_l += gradphi_l[i][2];
    }
    
    std::cout << "P型基函数导数之和： du=" << sum_du_p << " (应为0), dv=" << sum_dv_p 
              << " (应为0), dw=" << sum_dw_p << " (应为0)" << std::endl;
    std::cout << "L型基函数导数之和： du=" << sum_du_l << " (应为0), dv=" << sum_dv_l 
              << " (应为0), dw=" << sum_dw_l << " (应为0)" << std::endl;
}

/**
 * @brief 测试多个点
 */
void testMultiplePoints() {
    std::cout << "\n=== 多点测试 ===" << std::endl;
    
    // 测试多个点
    std::vector<std::vector<double>> test_points = {
        {0.0, 0.0, 0.0},    // 节点0
        {1.0, 0.0, 0.0},    // 节点1
        {0.0, 1.0, 0.0},    // 节点2
        {0.0, 0.0, 1.0},    // 节点3
        {1.0, 0.0, 1.0},    // 节点4
        {0.0, 1.0, 1.0}     // 节点5
    };
    
    for (int node = 0; node < 6; ++node) {
        double u = test_points[node][0];
        double v = test_points[node][1];
        double w = test_points[node][2];
        
        std::vector<double> phi_l(6);
        WedgeNodalLBasisAll(u, v, w, phi_l);
        
        std::cout << "在节点" << node << " (" << u << ", " << v << ", " << w << ")：" << std::endl;
        for (int i = 0; i < 6; ++i) {
            std::cout << "  phi_l[" << i << "] = " << phi_l[i];
            if (i == node) {
                std::cout << " (应为1)";
            } else {
                std::cout << " (应为0)";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

int main() {
    try {
        testWedgeBasisFunctions();
        testMultiplePoints();
        
        std::cout << "\n=== 测试完成 ===" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误：" << e.what() << std::endl;
        return 1;
    }
}