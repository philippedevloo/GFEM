# GFEM Project - AI Agent Instructions

## Project Overview
This project implements the **Generalized Finite Element Method (GFEM)** with enrichment functions for fracture mechanics simulations. It builds on the NeoPZ finite element library and extends it with specialized elements for modeling discontinuities and singularities near fracture tips in 2D and 3D domains.

## Build System & Dependencies

### Required NeoPZ Installation
This project depends on the NeoPZ library, which must be installed **first**:
```bash
# Build and install NeoPZ from neopz_GFEM workspace
cd ../neopz_GFEM
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=../neopz_install ..
make -j && make install
```

### GFEM Build Process
```bash
cd GFEM
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j
```

**Critical:** CMakeLists.txt searches for NeoPZ in `${CMAKE_SOURCE_DIR}/../neopz_install/` or `${CMAKE_SOURCE_DIR}/neopz_install/` (see GFEM/CMakeLists.txt line 25).

### Key Dependencies
- **NeoPZ::pz**: Core FEM library (must be installed separately)
- **Eigen3**: Linear algebra (matrix operations)
- **nlohmann_json**: JSON parsing for configuration files (e.g., `CrackDef.json`)
- **Catch2**: Unit testing framework
- **autodiff**: Automatic differentiation (fetched by CMake)

## Architecture

### Three-Level Mesh Hierarchy
GFEM uses a sophisticated multi-mesh architecture:

1. **Geometric Mesh (TPZGeoMesh)**: Contains geometry and topology
   - Read from `.msh` files using `TPZGmshReader`
   - Material IDs: `volmat=1`, `BCb=2`, `BCr=3`, `BCt=4`, `BCl=5`, `cutmat=6`, `fracedge=7`

2. **Computational Meshes**: Multiple meshes for different approximation spaces
   - **H1 Mesh**: Standard continuous FEM elements
   - **GFem Mesh (TPZGFemCompMesh)**: Enriched elements near fractures
   - **H1 Mirror Mesh**: For orthogonal enrichment schemes

3. **Multiphysics Mesh (TPZMultiphysicsCompMesh)**: Combines multiple computational meshes
   - Coordinates different approximation spaces
   - Manages connect restraints between enriched and standard elements

### Element Enrichment System

**Color-Based Element Classification:**
Elements are classified by their position relative to fractures:
- `Black`: Enriched elements (use custom shape functions)
- `White`: Standard elements (no enrichment)
- `Grey`: Transition elements

**Key Classes:**
- `TPZGFemCompElH1<TSHAPE>`: H1 element with GFEM enrichment
- `TPZGFemEnrichment`: Computes enrichment functions using SBFem eigenvalue analysis
- `TPZFrac2d`: Defines 2D fracture geometry and displacement jump functions
- `TPZGFemOrthogonal`: Orthogonal enrichment for specific problem types

### Material Models
GFEM extends NeoPZ materials with multi-space support:
- `TPZGFemElasticity2D`: 2D elasticity for enriched elements (inherits from `TPZMatCombinedSpacesT`)
- `TPZGFemElasticity3D`: 3D elasticity for enriched elements
- `TPZGFemDarcyFlow`: Darcy flow with discontinuities

All GFEM materials inherit from `TPZMatCombinedSpacesT` to handle multiple meshes.

## Workflow Patterns

### Typical Simulation Flow (see `2D/GFem2D.cpp::main`)
```cpp
1. Initialize refinement patterns: gRefDBase.InitializeRefPatterns(3)
2. Load geometric mesh from Gmsh
3. Uniform refinement
4. Directional refinement toward fracture: RefineTowardsFrac(gmesh, nref)
5. Build element coloring: BuildBlueRedElements(gmesh, blue, red)
6. Create H1 mesh: CreateH1CompMesh(gmesh)
7. Create GFem enriched mesh: CreateGFemCompMesh(gmesh)
8. Create multiphysics mesh: CreateMultiphysicsMesh(cmeshH1, cmeshGFem)
9. Compute connect restraints: ComputeConnectRestraints(cmesh_m)
10. Solve and visualize: Simulate(cmesh_m, plotfile)
```

### JSON Configuration Files
Enrichment and material properties are configured via JSON (e.g., `2D/CrackDef.json`):
```json
{
  "radius": 1.0,
  "angles": [0.0, 1.57, 3.141592653589793],
  "sector material ids": [1, 2],
  "problem_type": "Elastic",
  "materials": [...]
}
```

### Visualization Output
- Geometric meshes: `.vtk` files via `TPZVTKGeoMesh::PrintGMeshVTK()`
- Solution fields: `TPZVTKGenerator` with `res` parameter for subdivision
- Debug outputs: `.txt` files with `Print()` methods

## Project-Specific Conventions

### Directory Structure
- `2D/`: 2D GFEM simulations (GFem2D, SBFem2D, TestGFemEnrichment)
- `3D/`: 3D GFEM simulations (GFem3D)
- `NeoPZ_extensions/`: Core GFEM library (`GFEM_library` target)
  - Custom element types, materials, and enrichment classes
- `Meshes/`: Input mesh files (`.msh` format)
- `build/`: Out-of-source build directory

### Executables
Each subdirectory creates executables:
- `GFem2D`, `SBFem2D`: Main 2D simulation programs
- `TestGFemEnrichment`: Enrichment function testing
- `test_LinearElasticityConstitutive`: Unit tests with Catch2

### Material ID Convention
Material IDs are hardcoded constants (see `2D/GFem2D.cpp` lines 103-109):
```cpp
int volmat = 1;      // Volume domain
int BCb = 2;         // Bottom boundary
int BCr = 3;         // Right boundary
int BCt = 4;         // Top boundary
int BCl = 5;         // Left boundary
int cutmat = 6;      // Fracture cut
int fracedge = 7;    // Fracture edge/tip
```

### Simulation Types
Use enum to select physics (line 114):
```cpp
enum simultype {Darcy, Elast, DarcyNofrac, ElastNoFrac, DarcyOrthogonal};
```

## Common Pitfalls

1. **NeoPZ Must Be Installed First**: The build will fail if NeoPZ is not found. Check `find_package(NeoPZ REQUIRED ...)` paths.

2. **Mesh File Paths**: Executables read `.msh` files from the **binary directory**. CMakeLists.txt copies files using `configure_file(...COPYONLY)`.

3. **Connect Restraints**: For orthogonal enrichment, `ComputeConnectRestraints()` must be called before solving to properly couple enriched DOFs.

4. **Element Coloring**: The `BuildBlueRedElements()` function assigns colors based on fracture orientation. This affects which enrichment functions are applied.

5. **Debug Macros**: Code uses `DebugStop()` for assertions (NeoPZ convention). Enable with `-DPZDEBUG` for detailed checks.

## Testing & Debugging

### Running Tests
```bash
cd build
ctest  # or ./2D/test_LinearElasticityConstitutive
```

### Enable Logging
Compile NeoPZ with log4cxx support:
```bash
cmake -DUSING_LOG4CXX=ON ..
```
Then use `TPZLogger::InitializePZLOG()` in code.

### Mesh Visualization
View VTK outputs in ParaView to inspect:
- Mesh refinement quality
- Element colors (blue/red/white)
- Solution fields and discontinuities

## Key Files to Understand

- `NeoPZ_extensions/TPZGFemCompMesh.{h,cpp}`: Main enriched mesh class
- `NeoPZ_extensions/TPZGFemEnrichment.{h,cpp}`: Enrichment function computation
- `NeoPZ_extensions/TPZFrac2d.{h,cpp}`: Fracture displacement models
- `2D/GFem2D.cpp`: Complete 2D workflow example
- `3D/GFem3D.cpp`: Complete 3D workflow example
