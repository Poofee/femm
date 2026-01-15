@echo off
chcp 65001 >nul
echo Starting file organization by functional modules...

cd src

REM Move core infrastructure module files
echo Moving core infrastructure module files...

REM Base classes
move SolverBase.h core\base\
move SolverRegistry.h core\base\
move Types.h core\base\

REM Math computation
move LinearAlgebra.h core\math\
move LinearAlgebra.cpp core\math\
move CRSMatrix.h core\math\
move CRSMatrix.cpp core\math\
move Matrix.cpp core\math\
move IterativeSolver.h core\math\
move IterativeSolver.cpp core\math\

REM Mesh processing
move Mesh.h core\mesh\
move MeshIO.h core\mesh\
move MeshIO.cpp core\mesh\
move MeshUtils.h core\mesh\
move MeshUtils.cpp core\mesh\
move ElementDescription.h core\mesh\
move ElementDescription.cpp core\mesh\

REM Utility classes
move DefUtils.h core\utils\
move DefUtils.cpp core\utils\
move Lists.h core\utils\
move Lists.cpp core\utils\
move Integration.h core\utils\
move Integration.cpp core\utils\
move Interpolation.h core\utils\
move Interpolation.cpp core\utils\
move ShapeFunctions.h core\utils\
move ShapeFunctions.cpp core\utils\
move ElementUtils.h core\utils\
move ElementUtils.cpp core\utils\
move GaussIntegration.h core\utils\
move GaussIntegration.cpp core\utils\

REM Logging system
move LoggerInterface.h core\logging\
move SpdlogAdapter.h core\logging\
move LoggerFactory.h core\logging\
move LoggerFactory.cpp core\logging\

REM Move solver module files
echo Moving solver module files...

REM Electromagnetic solvers
move MagneticSolver.h solvers\electromagnetic\
move MagneticSolver.cpp solvers\electromagnetic\
move MagnetoDynamics2DSolver.h solvers\electromagnetic\
move MagnetoDynamics2DSolver.cpp solvers\electromagnetic\
move MagnetodynamicsSolver.h solvers\electromagnetic\
move ElectromagneticMaterial.h solvers\electromagnetic\

REM Thermal solvers
move HeatSolver.h solvers\thermal\
move HeatSolver.cpp solvers\thermal\
move HeatTransferSolver.h solvers\thermal\
move HeatTransferSolver.cpp solvers\thermal\

REM Structural solvers
move StructuralSolver.h solvers\structural\
move StructuralSolver.cpp solvers\structural\

REM Multiphysics solvers
move MultiphysicsSolver.h solvers\multiphysics\
move MultiphysicsSolver.cpp solvers\multiphysics\

REM Solver base classes
move NonlinearSolver.h solvers\base\

REM Move parallel computing module files
echo Moving parallel computing module files...

REM MPI parallel
move MPICommunicator.h parallel\mpi\
move MPIConfig.h parallel\mpi\
move MPIConfig.cpp parallel\mpi\
move MPIUtils.h parallel\mpi\
move DistributedLinearAlgebra.h parallel\mpi\
move DistributedLinearAlgebra.cpp parallel\mpi\

REM OpenMP parallel
move ParallelMatrixAssembly.h parallel\openmp\
move ParallelMatrixAssembly.cpp parallel\openmp\

REM Domain decomposition
move DomainDecomposition.h parallel\domain\
move DomainDecomposition.cpp parallel\domain\
move ParallelLinearSolver.h parallel\domain\
move ParallelLinearSolver.cpp parallel\domain\

REM Move I/O module files
echo Moving I/O module files...

REM Input processing
move InputFileParser.h io\input\
move InputFileParser.cpp io\input\

REM Material database
move Material.h io\material\
move MaterialDatabase.h io\material\

REM Matrix assembly
move MatrixAssembly.h io\matrix\
move MatrixAssembly.cpp io\matrix\

REM Move boundary conditions module files
echo Moving boundary conditions module files...

move BoundaryConditions.h boundary\
move BoundaryConditions.cpp boundary\

REM Move main solver files
echo Moving main solver files...

move ElmerSolver.h solvers\
move ElmerSolver.cpp solvers\

echo File organization completed!
cd ..
pause