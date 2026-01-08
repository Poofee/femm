@echo off
echo ================================================
echo Elmer FEM C++移植版 - 构建脚本
echo ================================================

setlocal

:: 检查编译器
where cl >nul 2>&1
if %errorlevel% neq 0 (
    echo 错误: 未找到Visual Studio编译器
    echo 请安装Visual Studio或设置VCVARS环境
    exit /b 1
)

:: 设置构建目录
set BUILD_DIR=build
if not exist %BUILD_DIR% mkdir %BUILD_DIR%
cd %BUILD_DIR%

:: 编译线性代数模块
echo 编译线性代数模块...
cl /c /EHsc /std:c++17 /I..\src ..\src\LinearAlgebra.cpp
if %errorlevel% neq 0 (
    echo 编译LinearAlgebra.cpp失败
    exit /b 1
)

:: 编译测试程序
echo 编译测试程序...
cl /c /EHsc /std:c++17 /I..\src ..\test\test_linear_algebra.cpp
if %errorlevel% neq 0 (
    echo 编译test_linear_algebra.cpp失败
    exit /b 1
)

:: 链接可执行文件
echo 链接可执行文件...
cl /EHsc /Fe:test_linear_algebra.exe test_linear_algebra.obj LinearAlgebra.obj
if %errorlevel% neq 0 (
    echo 链接失败
    exit /b 1
)

:: 运行测试
echo.
echo 运行线性代数模块测试...
test_linear_algebra.exe

if %errorlevel% equ 0 (
    echo.
    echo ================================================
    echo 构建成功！C++线性代数模块测试通过。
    echo ================================================
) else (
    echo.
    echo ================================================
    echo 测试失败，请检查代码。
    echo ================================================
)

cd ..
endlocal