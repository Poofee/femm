@echo off
chcp 65001 >nul
echo Cleaning build directory and reconfiguring...

REM 清理build目录
if exist "build" (
    echo Removing build directory...
    rmdir /s /q build
)

REM 创建新的build目录
mkdir build
cd build

REM 重新配置CMake
echo Configuring CMake...
cmake -G "Visual Studio 17 2022" -A x64 ..

REM 编译项目
echo Building project...
cmake --build . --config Release

echo Build process completed!
pause