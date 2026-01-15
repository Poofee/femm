@echo off
chcp 65001 >nul
echo Applying correct CMake configuration...

REM 使用正确的CMakeLists.txt
copy CMakeLists_correct.txt CMakeLists.txt >nul

echo Using correct CMakeLists.txt with proper include directories...

REM 清理并重新构建
if exist "build" (
    echo Removing build directory...
    rmdir /s /q build
)

mkdir build
cd build

echo Configuring CMake with proper include paths...
cmake -G "Visual Studio 17 2022" -A x64 ..

echo Building project...
cmake --build . --config Release

echo Build process completed!
pause