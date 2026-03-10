# CPP-SIMPLE Project Context

## Project Overview

**CPP-SIMPLE** is a 2D CFD (Computational Fluid Dynamics) solver implementing the SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm using the Finite Volume Method. It simulates fluid flow problems like lid-driven cavity flow on structured grids.

- **Language**: C++17
- **Build System**: CMake 3.14+
- **Project Type**: Scientific computing / CFD simulation

## Architecture

```
CPP-SIMPLE/
├── src/
│   ├── config/              # Configuration parameters
│   │   └── simulationConfig.hpp
│   ├── grid/                # Mesh generation and management
│   │   ├── structuredMesh.hpp
│   │   └── structuredMesh.cpp
│   ├── field/               # Field data containers
│   │   ├── scalarField.hpp/cpp      # Scalar data (temperature, pressure)
│   │   ├── vectorField.hpp/cpp      # Vector data (velocity)
│   │   ├── boundaryField.hpp/cpp    # Boundary conditions
│   │   └── fluidProperties.hpp/cpp  # Fluid properties
│   ├── equation/            # Equation discretization
│   │   ├── scalarEquation.hpp/cpp   # Scalar transport equation
│   ├── main.cpp             # Entry point
│   ├── test.cpp/hpp         # Test suite
├── utils/
│   └── formatter4eigen.h    # fmt formatting for Eigen matrices
├── data/                    # Output data directory
└── CMakeLists.txt
```

## Dependencies

| Package  | Purpose                          | CMake Target        |
|----------|----------------------------------|---------------------|
| Eigen3   | Linear algebra (matrix solvers)  | `Eigen3::Eigen`     |
| Boost    | multi_array, math constants      | `Boost::multi_array`, `Boost::math` |
| fmt      | Formatted output                 | `fmt::fmt`          |

**Package Manager**: vcpkg (recommended)

## Build Commands

```bash
# Create build directory
mkdir -p build && cd build

# Configure (ensure vcpkg toolchain is available)
cmake ..

# Build
make

# Run
./test_deps
```

## Key Configuration (simulationConfig.hpp)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Lx`, `Ly` | 1.0 | Cavity dimensions (m) |
| `ncx`, `ncy` | 10 | Grid cells in x/y |
| `Re` | 100.0 | Reynolds number |
| `lid_velocity` | 1.0 | Lid velocity (m/s) |
| `max_iterations` | 10000 | Solver max iterations |
| `continuity_tolerance` | 1e-6 | Convergence criterion |

## Core Classes

### StructuredMesh
Manages 2D structured grid geometry and topology.
- `createVolumeMesh()` - Generate cell centers and node coordinates
- `createBoundaryMesh()` - Identify boundary cell faces
- `saveMeshInfo()` - Export mesh to `data/mesh_info.json`

### ScalarField
Container for scalar data using `boost::multi_array<float, 2>`.
- Indexing: `field(i, j)` maps to internal `data_[j][i]` (x-direction contiguous)
- No boundary condition handling (pure data storage)

### VectorField
Velocity field composed of two ScalarField instances (u, v components).

### BoundaryField
Stores boundary conditions for all four walls (west, east, south, north).

### FluidPropertyField
Stores fluid properties (rho, mu, k, cp) with computed values (nu, alpha).

### ScalarEquation
Discretizes scalar transport equations (momentum, pressure correction, energy).
- Coefficients: aE, aW, aN, aS, aP, bsrc
- Methods: `resetCoefficients()`, `addDiffusionTerm()`, etc.

## Data Indexing Convention

- **Global coordinates**: x (0 to ncx-1, left to right), y (0 to ncy-1, bottom to top)
- **Array storage**: `data_[j][i]` - row-major, y-index outer, x-index inner
- **Access**: `field(i, j)` - intuitive (x, y) order
- **Boundaries**: West/East use j-index (0 to ncy-1), South/North use i-index (0 to ncx-1)

## Development Notes

1. **Boundary handling**: ScalarField/VectorField don't handle boundaries; BoundaryField provides BC data for equation assembly
2. **Constant properties**: FluidPropertyField currently uses constant properties; variable property update is TODO
3. **SIMPLE loop**: Main algorithm loop is not yet fully implemented
4. **Memory layout**: boost::multi_array uses row-major storage; x-direction access is cache-friendly

## Testing

Tests are in `src/test.cpp`:
- `test_basic_fields()` - Field initialization and access
- `test_linear_algebra()` - Eigen matrix operations
- `test_scalar_equation()` - Equation discretization tests

## Output

- Console: Simulation info, test results, timing
- File: `data/mesh_info.json` - Grid geometry data

## TODO / Extension Points

- Implement SIMPLE algorithm main loop
- Add pressure-velocity coupling (SIMPLE/SIMPLER/PISO)
- Integrate efficient linear solvers (BiCGSTAB, GMRES)
- Add time-stepping (explicit/implicit)
- Output to VTK/Tecplot formats
- Parallelization (OpenMP/MPI)
- Support non-uniform grids
