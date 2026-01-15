@echo off
chcp 65001 >nul
echo Rebuilding project...

REM 清理build目录
if exist "build" (
    echo Removing build directory...
    rmdir /s /q build
)

REM 创建新的build目录并配置
mkdir build
cd build

echo Configuring CMake...
cmake -G "Visual Studio 17 2022" -A x64 ..

echo Building project...
cmake --build . --config Release

REM 检查编译结果
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