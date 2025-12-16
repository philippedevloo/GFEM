# Quick Reference: Implementing TPZElasticity2D Tests

## Build & Run Commands

```bash
# Build
cd build && make test_TPZElasticity2D

# Run all tests
./2D/test_TPZElasticity2D

# Run with verbose output
./2D/test_TPZElasticity2D -s

# Run specific test case
./2D/test_TPZElasticity2D "Basic Material Properties"

# Run by tag
./2D/test_TPZElasticity2D "[constitutive]"

# List all tests
./2D/test_TPZElasticity2D --list-tests

# CTest integration
ctest -R TPZElasticity2DTests -V
```

## Catch2 Assertion Cheatsheet

```cpp
// Exact comparisons
REQUIRE(value == expected);
REQUIRE(mat.Dimension() == 2);

// Floating point comparisons
REQUIRE(result == Approx(expected).epsilon(1e-6));
REQUIRE(D(0,0) == Approx(210000.0).margin(0.01));

// Boolean checks
REQUIRE(mat.Id() == matid);
REQUIRE_FALSE(error_occurred);

// Exception handling
REQUIRE_THROWS(mat.SomeMethod());
REQUIRE_NOTHROW(mat.SafeMethod());
```

## TPZElasticity2D Common Methods

```cpp
// Constructor
TPZElasticity2D mat(matid, E, nu, fx, fy, planestress);

// Basic properties
mat.Id();              // Material ID
mat.Dimension();       // Should be 2
mat.NStateVariables(); // Number of state variables
mat.E();               // Young's modulus (if available)
mat.Nu();              // Poisson's ratio (if available)

// Check NeoPZ source for actual method names:
// - Contribute methods for element assembly
// - Solution methods for post-processing
// - Variable index methods for output
```

## Test Implementation Pattern

```cpp
SECTION("Your Test Description") {
    // 1. ARRANGE - Set up test data
    TPZElasticity2D mat(matid, E, nu, fx, fy, 1);
    TPZFNMatrix<9, STATE> matrix(3, 3, 0.);
    
    // 2. ACT - Execute the code under test
    // mat.YourMethodToTest(matrix);
    
    // 3. ASSERT - Verify results
    REQUIRE(matrix(0,0) == Approx(expected).epsilon(1e-6));
}
```

## Key Testing Areas

### ✅ Already Implemented
- Material construction
- Basic property checks (ID, dimension, state variables)

### ⏳ TODO - Implement These
1. **Constitutive Relations**
   - D matrix for plane stress
   - D matrix for plane strain
   - Compare with analytical values

2. **Boundary Conditions**
   - Create BC objects
   - Test BC application
   - Verify BC types (Dirichlet/Neumann)

3. **Element Contributions**
   - Stiffness matrix assembly
   - Load vector with body forces
   - Integration with TPZMaterialData

4. **Solution Variables**
   - Variable index mapping
   - Solution extraction
   - Post-processing output

5. **Edge Cases**
   - Invalid material parameters
   - Numerical limits
   - Error handling

## Where to Find Information

1. **NeoPZ Material API**: Check `neopz_GFEM/Material/Elasticity/TPZElasticity2D.h`
2. **Example Tests**: See `test_LinearElasticityConstitutive.cpp`
3. **Catch2 Docs**: https://github.com/catchorg/Catch2/tree/devel/docs
4. **Implementation Guide**: `TEST_IMPLEMENTATION_GUIDE.md`

## Expected D Matrix Values

### Plane Stress (planestress=1):
```cpp
REAL factor = E / (1.0 - nu * nu);
D11 = D22 = factor
D12 = D21 = factor * nu
D33 = E / (2.0 * (1.0 + nu))  // Shear modulus
```

### Plane Strain (planestress=0):
```cpp
REAL factor = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
D11 = D22 = factor
D12 = D21 = factor * nu / (1.0 - nu)
D33 = E / (2.0 * (1.0 + nu))  // Shear modulus
```

## Debugging Tips

```bash
# Run with debugger
lldb ./2D/test_TPZElasticity2D
(lldb) run

# Add breakpoints in your test
# Set PZDEBUG for NeoPZ debug checks
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Print during tests
std::cout << "Debug value: " << value << std::endl;
INFO("Debug value: " << value);  // Catch2 INFO macro
```

## Status Tracking

Mark completed sections:
- [x] Basic Properties
- [ ] Constitutive Matrix - Plane Stress
- [ ] Constitutive Matrix - Plane Strain
- [ ] Dirichlet BC
- [ ] Neumann BC
- [ ] Stiffness Contribution
- [ ] Load Vector Contribution
- [ ] Variable Names
- [ ] Post-Processing
- [ ] Edge Cases

---

**Start here:** Open `test_TPZElasticity2D.cpp` and replace the first `REQUIRE(true)` with actual test code!
