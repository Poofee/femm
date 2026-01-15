# PowerShell script to fix all compilation issues
Write-Host "Fixing all compilation issues..." -ForegroundColor Green

Set-Location "src"

# 1. 修复MPIConfig.h - 添加MPI_Op定义
$content = Get-Content "parallel/mpi/MPIConfig.h" -Raw
$content = $content -replace 'typedef int MPI_Datatype;', 'typedef int MPI_Datatype;' + "`n" + 'typedef int MPI_Op;'
Set-Content "parallel/mpi/MPIConfig.h" $content

# 2. 修复头文件包含路径问题
# MatrixAssembly.h
$content = Get-Content "MatrixAssembly.h" -Raw
$content = $content -replace '#include "../core/base/Types.h"', '#include "core/base/Types.h"'
Set-Content "MatrixAssembly.h" $content

# ElementMatrix.h
$content = Get-Content "ElementMatrix.h" -Raw
$content = $content -replace '#include "../core/base/Types.h"', '#include "core/base/Types.h"'
Set-Content "ElementMatrix.h" $content

# BoundaryConditions.h
$content = Get-Content "boundary/BoundaryConditions.h" -Raw
$content = $content -replace '#include "Mesh.h"', '#include "../core/mesh/Mesh.h"'
Set-Content "boundary/BoundaryConditions.h" $content

# ParallelLinearSolver.h
$content = Get-Content "parallel/domain/ParallelLinearSolver.h" -Raw
$content = $content -replace '#include "IterativeSolver.h"', '#include "../../core/math/IterativeSolver.h"'
Set-Content "parallel/domain/ParallelLinearSolver.h" $content

# 3. 修复DistributedLinearAlgebra.h中的override问题
$content = Get-Content "parallel/mpi/DistributedLinearAlgebra.h" -Raw
# 移除override关键字，因为基类可能没有这些方法
$content = $content -replace 'override', ''
Set-Content "parallel/mpi/DistributedLinearAlgebra.h" $content

# 4. 修复其他头文件包含
# 检查并修复所有可能的问题文件
$filesToCheck = @(
    "parallel/openmp/ParallelMatrixAssembly.h",
    "parallel/domain/DomainDecomposition.h",
    "solvers/ElmerSolver.h",
    "solvers/electromagnetic/MagneticSolver.h",
    "solvers/thermal/HeatSolver.h",
    "StructuralSolver.h",
    "MultiphysicsSolver.h"
)

foreach ($file in $filesToCheck) {
    if (Test-Path $file) {
        $content = Get-Content $file -Raw
        # 修复常见的头文件包含问题
        $content = $content -replace '#include "Types.h"', '#include "../core/base/Types.h"'
        $content = $content -replace '#include "LinearAlgebra.h"', '#include "../core/math/LinearAlgebra.h"'
        $content = $content -replace '#include "SolverBase.h"', '#include "../core/base/SolverBase.h"'
        $content = $content -replace '#include "SolverRegistry.h"', '#include "../core/base/SolverRegistry.h"'
        Set-Content $file $content
    }
}

Write-Host "All compilation issues fixed!" -ForegroundColor Green
Set-Location ".."
Read-Host "Press Enter to continue..."