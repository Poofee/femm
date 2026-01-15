@echo off
chcp 65001 >nul
echo Fixing compilation issues...

REM 1. 使用正确的CMakeLists.txt
echo Applying correct CMake configuration...
copy CMakeLists_correct.txt CMakeLists.txt >nul

REM 2. 清理build目录
echo Cleaning build directory...
if exist "build" (
    rmdir /s /q build
)

REM 3. 创建新的build目录并配置
mkdir build
cd build

echo Configuring CMake with proper include paths...
cmake -G "Visual Studio 17 2022" -A x64 ..

REM 4. 编译项目
echo Building project...
cmake --build . --config Release

REM 5. 检查编译结果
if %errorlevel% equ 0 (
    echo Build successful!
    echo.
    echo Checking generated executables...
    dir Release\*.exe
) else (
    echo Build failed with error code %errorlevel%
)

echo Build process completed!
pause