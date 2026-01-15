@echo off
chcp 65001 >nul
echo Final compilation fix...

REM 1. 修复头文件包含路径
cd src

echo Fixing header include paths...

REM 修复Mesh.h中的Types.h引用
powershell -Command "(Get-Content 'core\mesh\Mesh.h') -replace '#include \"Types.h\"', '#include \"../base/Types.h\"' | Set-Content 'core\mesh\Mesh.h'"

REM 修复BoundaryConditions.h中的Mesh.h引用
powershell -Command "(Get-Content 'boundary\BoundaryConditions.h') -replace '#include \"Mesh.h\"', '#include \"../core/mesh/Mesh.h\"' | Set-Content 'boundary\BoundaryConditions.h'"

REM 修复ElementMatrix.h中的头文件引用
powershell -Command "(Get-Content 'ElementMatrix.h') -replace '#include \"ShapeFunctions.h\"', '#include \"core\utils\ShapeFunctions.h\"' | Set-Content 'ElementMatrix.h'"
powershell -Command "(Get-Content 'ElementMatrix.h') -replace '#include \"GaussIntegration.h\"', '#include \"core\utils\GaussIntegration.h\"' | Set-Content 'ElementMatrix.h'"

REM 修复MatrixAssembly.h中的头文件引用
powershell -Command "(Get-Content 'MatrixAssembly.h') -replace '#include \"Types.h\"', '#include \"core\base\Types.h\"' | Set-Content 'MatrixAssembly.h'"

REM 修复其他核心头文件
powershell -Command "(Get-Content 'core\math\IterativeSolver.h') -replace '#include \"Types.h\"', '#include \"../base\Types.h\"' | Set-Content 'core\math\IterativeSolver.h'"
powershell -Command "(Get-Content 'core\utils\ShapeFunctions.h') -replace '#include \"Types.h\"', '#include \"../base\Types.h\"' | Set-Content 'core\utils\ShapeFunctions.h'"

cd ..

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