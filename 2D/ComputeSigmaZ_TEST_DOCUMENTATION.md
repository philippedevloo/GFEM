# ComputeSigmaZ Test Documentation

## Test Added: TPZLinearElasticityConstitutive - ComputeSigmaZ Verification

### Purpose
This test verifies the correctness of the `ComputeSigmaZ` method by comparing 2D plane strain analysis with full 3D elasticity analysis under plane strain conditions (epsilon_z = 0).

### Test Strategy

The test creates:
1. **3D orthotropic material** with randomly oriented principal axes (using `SetPrincipalAxes2D`)
2. **Equivalent 2D plane strain material** with the same properties and orientation
3. **Random in-plane strain** with epsilon_z = 0 (plane strain condition)

Then it:
- Applies the strain to the 3D material and computes the full 3D stress tensor
- Extracts sigma_z from the 3D analysis
- Applies the in-plane strain to the 2D material and computes in-plane stresses
- Uses `ComputeSigmaZ` to calculate the out-of-plane stress
- Compares all results for consistency

### Test Sections

#### 1. Orthotropic Material with Random 2D Rotation
- Tests the primary use case with full orthotropic properties
- Random rotation angle for principal axes orientation
- Verifies both in-plane stress match and sigma_z calculation

#### 2. Isotropic Material
- Simplified case with isotropic properties
- Verifies analytical formula: sigma_z = nu * (sigma_xx + sigma_yy)
- No rotation needed due to isotropy

#### 3. Plane Stress Check
- Verifies that `ComputeSigmaZ` returns zero for plane stress conditions
- Edge case verification

#### 4. Multiple Random Rotation Angles
- Robustness test with 10 different random configurations
- Tests various orientations and strain states
- Ensures consistency across different material orientations

#### 5. Transversely Isotropic Material
- Material with xy-plane isotropy and different z properties
- Common engineering material model
- Additional validation of orthotropic handling

### Key Verification Points

✅ **In-plane stress consistency**: 2D plane strain stresses match 3D analysis  
✅ **ComputeSigmaZ accuracy**: Method output matches sigma_z from 3D analysis  
✅ **Rotation handling**: Correct behavior with arbitrary 2D rotations  
✅ **Material types**: Works for isotropic, orthotropic, and transversely isotropic  
✅ **Plane stress**: Returns zero as expected  

### Mathematical Background

For plane strain (epsilon_z = 0), the out-of-plane stress is:

**Isotropic:**
```
sigma_z = nu * (sigma_xx + sigma_yy)
```

**Orthotropic (in principal axes):**
```
sigma_z = Ez * (epsilon_zz + nu_yz * epsilon_yy + nu_xz * epsilon_xx)
where epsilon_zz = 0 for plane strain
sigma_z = Ez * (nu_yz * epsilon_yy + nu_xz * epsilon_xx)
```

### Test Tolerance

- In-plane stress comparison: ±1e-4
- Sigma_z comparison: ±1e-4 (orthotropic), ±1e-6 (isotropic)
- Plane stress: ±1e-12

### Build and Run

```bash
cd build
make test_TPZElasticity2D
./2D/test_TPZElasticity2D "[sigmaz]"
```

Or run all ComputeSigmaZ tests:
```bash
./2D/test_TPZElasticity2D "[constitutive][sigmaz]"
```

### Related Methods Tested

- `SetOrthotropicProperties()` - Setting material properties
- `SetPrincipalAxes2D()` - Setting 2D rotation of principal axes
- `SetPlaneStress()` - Toggling plane stress/strain
- `ComputeStiffnessMatrix()` - Building 3D and 2D stiffness matrices
- `ComputeSigmaZ()` - Main method under test

### Complementary Test

This test complements the `ComputeStrainZ` test which verifies:
- Plane stress (sigma_z = 0) → calculates epsilon_z
- Uses 3D compliance matrix

While this test verifies:
- Plane strain (epsilon_z = 0) → calculates sigma_z
- Uses 3D stiffness matrix

Together they verify the complete 2D/3D consistency for both plane stress and plane strain conditions.
