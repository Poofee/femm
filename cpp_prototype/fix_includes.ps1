# PowerShell script to fix header include paths
Write-Host "Fixing header include paths..." -ForegroundColor Green

Set-Location "src"

# 修复并行计算模块的头文件包含
# MPI模块
$content = Get-Content "parallel/mpi/DistributedLinearAlgebra.h" -Raw
$content = $content -replace '#include "MPIConfig.h"', '#include "../../parallel/mpi/MPIConfig.h"'
$content = $content -replace '#include "LinearAlgebra.h"', '#include "../../core/math/LinearAlgebra.h"'
Set-Content "parallel/mpi/DistributedLinearAlgebra.h" $content

$content = Get-Content "parallel/openmp/ParallelMatrixAssembly.h" -Raw
$content = $content -replace '#include "DistributedLinearAlgebra.h"', '#include "../../parallel/mpi/DistributedLinearAlgebra.h"'
Set-Content "parallel/openmp/ParallelMatrixAssembly.h" $content

$content = Get-Content "parallel/domain/DomainDecomposition.h" -Raw
$content = $content -replace '#include "MPIConfig.h"', '#include "../../parallel/mpi/MPIConfig.h"'
Set-Content "parallel/domain/DomainDecomposition.h" $content

$content = Get-Content "parallel/domain/ParallelLinearSolver.h" -Raw
$content = $content -replace '#include "DistributedLinearAlgebra.h"', '#include "../../parallel/mpi/DistributedLinearAlgebra.h"'
Set-Content "parallel/domain/ParallelLinearSolver.h" $content

# 修复核心模块的头文件包含
$content = Get-Content "MatrixAssembly.h" -Raw
$content = $content -replace '#include "Types.h"', '#include "../core/base/Types.h"'
Set-Content "MatrixAssembly.h" $content

$content = Get-Content "ElementMatrix.h" -Raw
$content = $content -replace '#include "Types.h"', '#include "../core/base/Types.h"'
Set-Content "ElementMatrix.h" $content

# 修复边界条件模块的头文件包含
$content = Get-Content "boundary/BoundaryConditions.h" -Raw
$content = $content -replace '#include "LinearAlgebra.h"', '#include "../core/math/LinearAlgebra.h"'
Set-Content "boundary/BoundaryConditions.h" $content

# 修复求解器模块的头文件包含
$content = Get-Content "solvers/ElmerSolver.h" -Raw
$content = $content -replace '#include "SolverBase.h"', '#include "../core/base/SolverBase.h"'
$content = $content -replace '#include "SolverRegistry.h"', '#include "../core/base/SolverRegistry.h"'
$content = $content -replace '#include "Types.h"', '#include "../core/base/Types.h"'
Set-Content "solvers/ElmerSolver.h" $content

# 修复电磁求解器的头文件包含
$content = Get-Content "solvers/electromagnetic/MagneticSolver.h" -Raw
$content = $content -replace '#include "SolverBase.h"', '#include "../../core/base/SolverBase.h"'
Set-Content "solvers/electromagnetic/MagneticSolver.h" $content

# 修复热求解器的头文件包含
$content = Get-Content "solvers/thermal/HeatSolver.h" -Raw
$content = $content -replace '#include "SolverBase.h"', '#include "../../core/base/SolverBase.h"'
Set-Content "solvers/thermal/HeatSolver.h" $content

# 修复结构求解器的头文件包含
$content = Get-Content "StructuralSolver.h" -Raw
$content = $content -replace '#include "SolverBase.h"', '#include "../core/base/SolverBase.h"'
Set-Content "StructuralSolver.h" $content

# 修复多物理场求解器的头文件包含
$content = Get-Content "MultiphysicsSolver.h" -Raw
$content = $content -replace '#include "SolverBase.h"', '#include "../core/base/SolverBase.h"'
Set-Content "MultiphysicsSolver.h" $content

Write-Host "Header include paths fixed!" -ForegroundColor Green
Set-Location ".."
Read-Host "Press Enter to continue..."