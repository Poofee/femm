@echo off
chcp 65001 >nul
echo Fixing compilation issues...

REM 1. 使用简化的CMakeLists.txt
copy CMakeLists_simple.txt CMakeLists.txt >nul

echo Using simplified CMakeLists.txt without MPI modules...

REM 2. 清理并重新构建
if exist "build" (
    echo Removing build directory...
    rmdir /s /q build
)

mkdir build
cd build

echo Configuring CMake...
cmake -G "Visual Studio 17 2022" -A x64 ..

echo Building project...
cmake --build . --config Release

echo Build process completed!
pause