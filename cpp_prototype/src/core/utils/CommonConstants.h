#ifndef ELMER_CPP_COMMON_CONSTANTS_H
#define ELMER_CPP_COMMON_CONSTANTS_H

#include <complex>

namespace elmer {

// =============================================================================
// 数学常数
// =============================================================================

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI_4
#define M_PI_4 0.78539816339744830962
#endif

#ifndef M_2_PI
#define M_2_PI 0.63661977236758134308
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#ifndef M_LOG2E
#define M_LOG2E 1.44269504088896340736
#endif

#ifndef M_LOG10E
#define M_LOG10E 0.43429448190325182765
#endif

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_LN10
#define M_LN10 2.30258509299404568402
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

// =============================================================================
// 物理常数
// =============================================================================

// 真空介电常数 (F/m)
constexpr double VACUUM_PERMITTIVITY = 8.854187817e-12;

// 真空磁导率 (H/m)
constexpr double VACUUM_PERMEABILITY = 4.0 * M_PI * 1e-7;

// 光速 (m/s)
constexpr double SPEED_OF_LIGHT = 299792458.0;

// 玻尔兹曼常数 (J/K)
constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;

// 普朗克常数 (J·s)
constexpr double PLANCK_CONSTANT = 6.62607015e-34;

// 阿伏伽德罗常数 (1/mol)
constexpr double AVOGADRO_CONSTANT = 6.02214076e23;

// 电子电荷 (C)
constexpr double ELECTRON_CHARGE = 1.602176634e-19;

// 电子质量 (kg)
constexpr double ELECTRON_MASS = 9.10938356e-31;

// =============================================================================
// 数值计算常数
// =============================================================================

// 机器精度
constexpr double MACHINE_EPSILON = 1e-15;

// 浮点数比较容差
constexpr double FLOATING_POINT_TOLERANCE = 1e-12;

// 最大迭代次数
constexpr int MAX_ITERATIONS = 1000;

// 收敛容差
constexpr double CONVERGENCE_TOLERANCE = 1e-8;

// =============================================================================
// 电磁学常数
// =============================================================================

// 自由空间波阻抗 (Ω)
constexpr double FREE_SPACE_IMPEDANCE = 376.730313668;

// 磁通量量子 (Wb)
constexpr double MAGNETIC_FLUX_QUANTUM = 2.067833848e-15;

// 约瑟夫森常数 (Hz/V)
constexpr double JOSEPHSON_CONSTANT = 483597.8484e9;

// =============================================================================
// 单位转换常数
// =============================================================================

// 角度转换
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// 温度转换
constexpr double CELSIUS_TO_KELVIN = 273.15;

// 长度转换
constexpr double INCH_TO_METER = 0.0254;
constexpr double FOOT_TO_METER = 0.3048;
constexpr double MILE_TO_METER = 1609.344;

// =============================================================================
// 工程常数
// =============================================================================

// 标准大气压 (Pa)
constexpr double STANDARD_ATMOSPHERE = 101325.0;

// 标准重力加速度 (m/s²)
constexpr double STANDARD_GRAVITY = 9.80665;

// 水的密度 (kg/m³)
constexpr double WATER_DENSITY = 1000.0;

// 空气密度 (kg/m³)
constexpr double AIR_DENSITY = 1.225;

// =============================================================================
// 宏定义工具
// =============================================================================

// 数组长度
#define ARRAY_LENGTH(arr) (sizeof(arr) / sizeof((arr)[0]))

// 最小值/最大值
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// 绝对值
#define ABS(x) (((x) < 0) ? -(x) : (x))

// 平方和立方
#define SQUARE(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))

// 符号函数
#define SIGN(x) (((x) > 0) ? 1 : (((x) < 0) ? -1 : 0))

// 角度转换宏
#define DEG2RAD(x) ((x) * DEG_TO_RAD)
#define RAD2DEG(x) ((x) * RAD_TO_DEG)

// =============================================================================
// 类型定义
// =============================================================================

// 浮点数类型别名
using Real = double;
using Float = float;

// 复数类型别名
using Complex = std::complex<double>;
using ComplexF = std::complex<float>;

// 向量和矩阵类型别名
using RealVector = std::vector<double>;
using FloatVector = std::vector<float>;
using ComplexVector = std::vector<Complex>;
using ComplexVectorF = std::vector<ComplexF>;

} // namespace elmer

#endif // ELMER_CPP_COMMON_CONSTANTS_H