# Horizontal-resolution sensitivity of Antarctic landfast ice in standalone CICE

## Purpose

This project tests whether the Antarctic landfast-ice solution obtained with the existing **0.25-degree free-slip/static-lateral-drag configuration** is robust to a large change in horizontal model resolution. Two additional standalone CICE configurations are proposed at the following degree resolutions:

```text
coarse resolution:       1.0
existing reference:      0.25
high resolution:         0.10
```

The first experiments should deliberately preserve the physical parameter set of the existing `Cs-high` experiment:

```fortran
boundary_condition = 'free_slip'
form_func           = 'static'
Cs                  = 1.0e-3
Ktens               = 0.2
e_yieldcurve        = 1.5
e_plasticpot        = 1.5
```

along with the same ice-strength, thermodynamic and forcing choices used by the 0.25° reference experiment. The primary question is therefore **not** initially:

> What parameters produce the best fast-ice area at each resolution?

It is:

> **Does the physical parameter set that produces a realistic Antarctic fast-ice solution at 0.25° remain effective when the model grid is made substantially coarser or finer?**

Only after this parameter-transfer experiment should resolution-specific adjustment of `Cs`, `Ktens`, `e_yieldcurve`, or `e_plasticpot` be considered.

---

## Scientific rationale

Antarctic landfast ice is an unusually strong test of horizontal-resolution dependence because its existence is controlled by interactions between processes occurring at very different spatial scales:

- coastline and embayment geometry;
- grounded icebergs and other pinning points;
- sea-ice stress transmission;
- tensile failure;
- shear and convergence;
- atmospheric stress;
- ocean stress;
- ice thickness and strength;
- lateral resistance from unresolved coastal geometry.

Changing horizontal resolution therefore changes more than the graphical representation of the ice edge. It changes the geometry on which the momentum equation is solved, the spatial gradients used to calculate deformation, the representation of coastal obstacles, the area represented by an individual fast-ice grid cell, and potentially the numerical convergence of the EVP solution. The proposed experiment is consequently a test of the **scale robustness of the fast-ice parameterisation system**.

### Why 1- and 0.1-degree?

These resolutions deliberately bracket the present 0.25° configuration over a large range. Relative to 0.25°:

```text
1.0°    = 4 × coarser linear grid spacing
         ≈ 16 × larger nominal grid-cell area

0.10°   = 2.5 × finer linear grid spacing
         ≈ 6.25 × smaller nominal grid-cell area
```

The full range from 1° to 0.1° is a factor of ten in nominal linear resolution. This is large enough that genuine resolution dependence should become distinguishable from small numerical differences between otherwise similar model grids. The two experiments also address different scientific questions.

### 1° experiment

The 1° configuration asks whether the lateral-drag/rheology system can represent Antarctic fast ice when much of the relevant coastal geometry is necessarily subgrid. At this resolution:

- narrow coastal embayments may occupy only part of a grid cell;
- island chains and grounded-iceberg complexes are increasingly subgrid;
- the fast-ice belt may be only one or a few model cells wide;
- atmospheric and ocean forcing are strongly spatially averaged;
- individual binary fast-ice classifications represent large areas;
- deformation and velocity gradients are substantially smoothed.

This is an important limit because it tests whether the form-factor parameterisation behaves as intended when unresolved coastal geometry becomes dominant.

### 0.1° experiment

The 0.1° configuration asks the complementary question. More coastal geometry becomes explicitly resolved, the fast-ice zone spans more cells, spatial gradients become sharper, and the model can represent substantially more structure in ice velocity, deformation and ocean forcing. The parameterisation must therefore operate in an environment in which some processes previously represented through effective subgrid resistance begin to become resolved. The high-resolution experiment tests whether the 0.25° solution remains valid as CICE approaches a much more explicitly resolved coastal mechanical system.

---

## Primary scientific questions

The first experiment set should address five questions.

### Q1. Is the 0.25° `Cs-high` solution resolution transferable?

With the same physical parameter values, do the 1°, 0.25° and 0.1° configurations produce broadly similar:

```text
circumpolar FIA
regional FIA
seasonal FIA cycle
fast-ice probability
fast-ice onset
seasonal maximum
breakout timing
```

relative to AF2020?

A positive result would be scientifically important because it would indicate that the parameterisation represents a reasonably scale-robust physical mechanism rather than compensating for a particular model grid.

### Q2. Are the rheological parameters resolution transferable?

The existing configuration uses:

```fortran
Ktens        = 0.2
e_yieldcurve = 1.5
e_plasticpot = 1.5
```

`Ktens` introduces tensile strength relative to compressive ice strength, while the elliptical yield-curve parameters control the relationship between the components of the internal stress state.

These are constitutive parameters.

There is therefore no strong reason to retune them *a priori* simply because grid spacing changes.

However, the discretised rheology is not resolution independent. Changing grid spacing changes resolved velocity gradients and strain rates and may alter EVP convergence. The same constitutive parameters can consequently produce a different effective mechanical response.

Testing them unchanged is therefore scientifically more useful than immediately retuning them.

### Q3. Is `Cs = 1.0e-3` resolution transferable?

The static lateral-drag coefficient should be treated separately from the form factor.

In the current implementation:

```text
Cs  -> strength scale of the static lateral-resistance closure
F2  -> spatial/geometrical weighting of that resistance
```

The CICE implementation describes the static branch using approximately:

```text
phi_static = Cs / (|u| + u0)
```

while the lateral-drag stress factor includes the local form factor:

```text
Ku ~ ice_mass × F2
```

`Cs` therefore does not directly encode the target-grid coastline geometry.

The clean first hypothesis is consequently:

> **If `F2` is generated consistently for each grid, the same `Cs` should first be tested at all three resolutions.**

Whether `Cs = 1.0e-3` actually remains appropriate is then an empirical result rather than an assumption.

### Q4. Does increasing resolution change FIA for the correct reason?

Similar total FIA does not necessarily imply equivalent physics.

For example, excessive lateral drag and excessive internal tensile strength can both produce a similar circumpolar FIA while producing very different:

```text
fast-ice probability maps
strain-rate fields
ice-strength fields
velocity distributions
breakout locations
regional seasonal cycles
```

The experiments must therefore diagnose the mechanism producing FIA rather than using FIA alone as the calibration target.

### Q5. How much of the 0.25° skill arises from unresolved geometry?

This may ultimately be the most interesting question. If substantially different resolutions reproduce similar AF2020 fast-ice distributions using the same `Cs`, `Ktens` and elliptical yield curve, the result would support the physical portability of the parameterisation. If parameter changes are required, their direction and magnitude can instead reveal how much of the 0.25° solution represents compensation for unresolved coastline, forcing or rheological structure.

---

# Stage R0: build resolution-specific model inputs

Before running the resolution experiments, all grid-dependent inputs must be constructed consistently.

## R0.1 CICE grids and masks

Independent CICE grids are required for:

```text
1.0°
0.10°
```

including the corresponding:

```text
T/U/E/N coordinates
grid metrics
grid-cell areas
land/ocean masks
bathymetry
C-grid masks
grid angles
```

The target-grid land mask should be generated consistently from the same high-resolution source geometry wherever possible. This is particularly important because changing the resolved coastline while simultaneously changing the form-factor source would make attribution difficult.

---

## R0.2 ERA5 atmospheric forcing

ERA5 must be remapped independently to each CICE grid. The forcing preprocessing should retain exactly the same physical variables and temporal sampling used by the 0.25° experiment. For vector variables:

```text
eastward/northward winds
        ↓
spatial remapping
        ↓
rotation onto target CICE grid
```

For scalar and flux variables, the same interpolation or conservative-remapping conventions should be used across resolutions.

### Important interpretation at 0.1°

ERA5 does not acquire additional atmospheric information simply because CICE is run at 0.1°. The 0.1° model may resolve finer sea-ice dynamics and coastline geometry, but the atmospheric forcing remains limited by the effective spatial resolution of ERA5. That is not a defect in the experiment. It simply means that improvements between 0.25° and 0.1° should not automatically be interpreted as resulting from better-resolved atmospheric forcing.

---

## R0.3 ORAS ocean forcing

ORAS fields should likewise be remapped independently to each target CICE grid. Particular care is required for:

```text
uocn
vocn
SST
SSS
mixed-layer depth
other ocean fields used by the standalone forcing system
```

Vector currents should be consistently rotated onto the target CICE grid. The 0.1° experiment is particularly interesting because the ORAS product used here contains substantially finer native spatial information than is retained by the existing 0.25° CICE grid. Consequently, the high-resolution experiment can retain more of the resolved ocean-current and coastal ocean structure. This means that the proposed experiments test the **effective resolution of the complete forced CICE system**, rather than an abstract numerical grid in isolation. If the 0.1° experiment later shows a large response, an optional follow-on experiment could impose a common 0.25° spatial bandwidth on the ORAS forcing before remapping it to 0.1°. That would help separate:

```text
CICE resolution effect

from

additional ocean-forcing resolution
```

but this is not required for the first experiment.

---

## R0.4 Form-factor fields

The form-factor NetCDF requires special treatment. The production 1° and 0.1° form-factor products should **not simply be obtained by bilinear interpolation of the existing 0.25° `F2` field**. Instead, the same high-resolution coastline and obstacle geometry used for the 0.25° product should be processed independently onto each target CICE grid.

Conceptually:

```text
high-resolution coastline / grounded-obstacle geometry
                  |
                  +------> 1.0° F2
                  |
                  +------> 0.25° F2
                  |
                  +------> 0.10° F2
```

using the same form-factor algorithm in every case. This distinction matters because the form factor represents coastline geometry relative to the dimensions and orientation of the target model cell. Changing target-grid spacing therefore changes the meaning of the aggregated geometry. The resulting files should preserve:

```text
F2x
F2y
target-grid coordinates
target-grid mask
source geometry version
form-factor algorithm version
normalisation convention
```

and use the same CICE:

```fortran
F2_map_method
```

at every resolution.

### Required F2 diagnostics

Before any CICE integration, compare the three form-factor products using:

```text
F2x/F2y maps
F2 magnitude maps
non-zero F2 area
number of affected cells
area-weighted F2 distribution
regional F2 distributions
integrated coastline/form-factor measure
```

A large change in these statistics is not necessarily an error. It may be the physically expected consequence of changing which coastline structure is resolved versus subgrid. It must, however, be quantified before interpreting changes in FIA.

---

## R0.5 Initial conditions and restart state

Any initial or restart state must also be compatible with the new grids. Where practical, all three resolutions should represent the same physical initial state. Restart remapping should preserve, as far as possible:

```text
ice area
ice volume
snow volume
category distribution
thermodynamic state
```

The first comparison period should not be interpreted until differences introduced purely by initial-condition remapping have sufficiently decayed.

---

# Stage R1: fixed-physics resolution experiment

The first science experiment should contain **no resolution-specific physical tuning**.

| Experiment | Resolution | `Ktens` | `e_yieldcurve` | `e_plasticpot` | `Cs` | Form factor |
|---|---:|---:|---:|---:|---:|---|
| reference | 0.25° | 0.2 | 1.5 | 1.5 | `1.0e-3` | native 0.25° F2 |
| R1-coarse | 1.0° | 0.2 | 1.5 | 1.5 | `1.0e-3` | regenerated 1° F2 |
| R1-high | 0.10° | 0.2 | 1.5 | 1.5 | `1.0e-3` | regenerated 0.1° F2 |

All other physical namelist choices should remain those of the 0.25° `Cs-high` reference wherever technically possible. This includes the same:

```text
ice-strength formulation
Pstar
Cstar
thermodynamics
atmospheric forcing product
ocean forcing product
free-slip boundary treatment
static form function
fast-ice classification algorithm
analysis years
```

The experiment therefore asks whether the existing physical configuration transfers across resolution without calibration.

---

# Numerical parameters that require reconsideration

Holding the **physics** fixed does not imply holding every numerical parameter fixed. Some namelist choices may need adjustment to ensure that the three integrations solve approximately the same physical equations.

## Model timestep

The current external model timestep should be tested at both new resolutions. The 1° experiment will generally be less restrictive. The 0.1° experiment may require a shorter timestep because finer spatial scales permit:

```text
larger local velocity gradients
larger deformation gradients
stronger local forcing gradients
more restrictive transport CFL conditions
```

A timestep reduction should be considered a numerical adjustment, not physical tuning.

---

## EVP subcycling: `ndte`

This requires particular attention. For standard EVP:

```text
EVP subcycle timestep ~ dt / ndte
```

Changing grid spacing while retaining the same `dt` and `ndte` does not guarantee equivalent numerical convergence. The 0.1° simulation should therefore undergo short convergence tests using at least two `ndte` values before the production integration is accepted. The objective is not to find an `ndte` that gives desirable FIA. The objective is to demonstrate that the resulting velocity, deformation and fast-ice solution is not materially controlled by inadequate EVP subcycling. This distinction is essential because otherwise apparent rheological resolution dependence could actually be solver-resolution dependence.

---

## EVP relaxation and regularisation parameters

The following should initially remain unchanged but be explicitly checked:

```fortran
arlx
brlx
elasticDamp
deltaminEVP
deltaminVP
```

They should only be changed if numerical diagnostics demonstrate that the existing values no longer provide an adequately converged or stable EVP solution. They should **not** be altered to improve FIA agreement.

---

## Transport

The same transport scheme should be retained. However, CFL behaviour and any resolution-dependent numerical diffusion should be diagnosed, particularly at 0.1°.

---

## Parallel decomposition

Changes to:

```text
nx_global
ny_global
block_size_x
block_size_y
processor decomposition
```

are computational rather than scientific and should be optimised independently for each grid. They should nevertheless be documented so that performance differences are not confused with physics changes.

---

# Expected physical effects of changing resolution

No monotonic relationship between grid resolution and FIA should be assumed. Several competing processes operate simultaneously.

## Expected behaviour at 1°

Coarser resolution may favour fast ice because:

```text
winds and ocean currents are spatially smoother
local stress maxima are reduced
velocity gradients are weaker
small mobile regions may disappear within coarse cells
```

but may suppress fast ice because:

```text
coastal embayments are poorly represented
narrow fast-ice corridors disappear
pinning geometry becomes increasingly subgrid
coastal stress transmission is poorly resolved
```

There is also a substantial sampling issue. At 1°, changing the state of one grid cell from mobile to fast can add or remove a very large physical area from FIA. Consequently, circumpolar FIA may change abruptly even where the underlying velocity difference is modest.

---

## Expected behaviour at 0.1°

Higher resolution may improve fast-ice representation because:

```text
coastal geometry is better resolved
fast-ice belts span more model cells
regional differences are better represented
ORAS ocean structure is better retained
fast-ice edges can occupy more realistic positions
```

but it may also make fast ice more difficult to maintain because:

```text
stress gradients are sharper
deformation becomes more spatially heterogeneous
localised failure becomes better resolved
ocean-current gradients are stronger
small mobile corridors can exist within formerly fast 0.25° cells
```

A lower FIA at 0.1° would therefore not automatically imply that the high-resolution experiment is worse. It could indicate that the 0.25° model was artificially locking spatially heterogeneous coastal ice into single large cells. Conversely, increased FIA could indicate that better-resolved geometry permits stronger mechanical anchoring. The spatial diagnostics are therefore as important as total area.

---

# Required fast-ice diagnostics

## Primary observational metrics

All three resolutions should be evaluated against AF2020 using:

```text
circumpolar FIA seasonal cycle
regional FIA seasonal cycles
annual/seasonal mean FIA
seasonal maximum FIA
timing of growth
timing of maximum extent
timing of breakout
fast-ice probability
```

The existing 0.25° `Cs-high` experiment remains the model reference.

---

## Resolution-aware AF2020 comparison

AF2020 should not simply be converted into an independent binary observation at each coarse model cell without retaining information about subgrid fast-ice coverage. For each target model grid, it is preferable to calculate the **fraction of the model-cell area classified as fast ice by AF2020**. This is especially important at 1°. A model cell containing 20% observed fast ice should not be observationally equivalent to one containing 100% fast ice simply because both contain some observed landfast ice. Two complementary evaluations are recommended:

```text
1. observations conservatively aggregated to each native model grid;

2. all experiments remapped to a common evaluation grid for direct
   inter-resolution FIP comparison.
```

Circumpolar and regional FIA should always be calculated using physical grid-cell area.

---

# Mechanical attribution diagnostics

FIA alone cannot determine whether the appropriate mechanism is being represented. At minimum, compare:

```text
ice speed
ice strength
strain-rate invariant
divergence
shear
lateral-drag stress
F2
fast-ice thickness
```

between the three resolutions. Particularly useful quantities are:

```text
F2-conditioned ice speed
F2-conditioned lateral-drag stress
fast vs mobile strain-rate distributions
distance-from-coast dependence
regional lateral-drag stress
regional internal stress
```

These diagnostics provide a way to distinguish the roles of the three principal parameter families.

### `Cs`

Primarily changes the magnitude of the imposed lateral resistance where form factors are active.

### `Ktens`

Changes the ability of the ice cover to transmit stress while resisting tensile failure.

### `e_yieldcurve` and `e_plasticpot`

Change the relationship between deformation and allowable internal stress within the elliptical rheology.

Therefore:

> **Similar FIA produced by different combinations of `Cs`, `Ktens` and ellipse aspect ratio should not be considered mechanically equivalent.**

---

# Interpretation of the R1 experiment

| Result | Interpretation | Next step |
|---|---|---|
| 1°, 0.25° and 0.1° all reproduce similar FIA and FIP | strong evidence that the parameterisation is scale robust | retain common parameters |
| FIA similar but FIP substantially different | compensating regional errors | diagnose geometry and stress before tuning |
| 1° differs strongly but 0.25° and 0.1° agree | coarse-grid representation becomes limiting | investigate 1° F2/geometry; avoid changing rheology first |
| 0.1° differs strongly but 1° and 0.25° agree | 0.25° parameters may include unresolved-scale compensation, or EVP convergence may differ | test numerical convergence and F2 before retuning |
| both new resolutions differ systematically in the same direction | likely parameter or preprocessing dependence tied to the 0.25° configuration | inspect common assumptions and form-factor scaling |
| changing `ndte` materially changes FIA | numerical convergence problem | resolve before interpreting physical parameters |
| lateral-drag stress changes strongly while internal-stress diagnostics remain similar | likely `F2`/`Cs` scale dependence | test `Cs` |
| internal deformation changes strongly away from F2 regions | likely rheology/resolution interaction | investigate `Ktens` and ellipse parameters |

---

# Stage R2: resolution-specific parameter sensitivity

Only after Stage R1 has been completed should physical parameters be adjusted.

The order of investigation should be deliberate.

## R2.1 Numerical convergence first

Before changing physical parameters:

```text
verify timestep
verify EVP subcycling
verify forcing interpolation
verify F2 construction
verify restart behaviour
```

A numerically unconverged 0.1° experiment should not be corrected by changing rheology.

---

## R2.2 Test `Cs`

`Cs` should probably be the first physical parameter tested if the primary difference is concentrated around form-factor regions.

A simple logarithmic bracket around the reference value is preferable to immediately searching a broad parameter space:

```text
weaker lateral resistance
Cs = 0.5 × reference

reference
Cs = 1.0 × reference

stronger lateral resistance
Cs = 2.0 × reference
```

The purpose is to determine whether the effective lateral resistance associated with a grid-specific `F2` field remains approximately invariant.

If a different `Cs` is required at every resolution, the relationship between `Cs`, `F2` and grid spacing becomes a scientific result in its own right.

---

## R2.3 Test tensile strength

If FIA differences are associated with tensile failure and stress transmission rather than only the coastal form-factor region, `Ktens` should then be tested independently.

The central question becomes:

> Does a common ratio of tensile-to-compressive strength generate similar fast-ice stability at different horizontal resolutions?

---

## R2.4 Test elliptical yield-curve aspect ratio

Only after separating the effects of lateral resistance and tensile strength should the elliptical yield curve be varied.

The key question is whether the relationship between shear/compressive deformation represented by:

```fortran
e_yieldcurve
e_plasticpot
```

must change as the model resolves progressively finer deformation structures.

These tests should be interpreted using deformation and stress diagnostics, not FIA alone.

---

# Parameter-identifiability problem

The experiment should explicitly avoid simultaneous tuning of:

```text
Cs
Ktens
e_yieldcurve
e_plasticpot
```

because these parameters can partially compensate for one another.

For example, realistic FIA could potentially be generated by:

```text
too much coastal resistance + weak internal ice

or

too little coastal resistance + overly strong tensile ice
```

Both can produce the same integrated area while representing different mechanics.

The resolution experiments provide an opportunity to distinguish these mechanisms precisely because coastal geometry and resolved deformation change so substantially between 1° and 0.1°.

---

# Recommended first experiment hierarchy

```text
R0
|
|-- build 1° grid
|-- build 0.1° grid
|
|-- regrid ERA5 independently
|-- regrid ORAS independently
|-- regenerate F2 independently
|-- prepare equivalent initial states
|
v
R1
|
|-- 1°   : original Cs-high physics
|-- 0.25°: existing Cs-high reference
|-- 0.1° : original Cs-high physics
|
|-- verify numerical convergence
|-- compare AF2020 FIA/FIP
|-- compare mechanical diagnostics
|
v
DECISION
|
|-- geometry/lateral-drag difference -> test Cs
|
|-- internal failure difference      -> test Ktens
|
|-- deformation/yield difference     -> test ellipse
|
`-- forcing-resolution ambiguity     -> common-bandwidth forcing test
```

---

# Success criteria

The first resolution experiment is successful if it establishes, rather than assumes, whether the present fast-ice parameterisation is scale robust.

The experiment should answer:

1. Can the existing 0.25° `Cs-high` physical configuration generate realistic Antarctic fast ice at 1°?

2. Can the same configuration generate realistic Antarctic fast ice at 0.1°?

3. Does `Ktens = 0.2` remain an appropriate tensile-strength ratio across the three resolutions?

4. Do `e_yieldcurve = e_plasticpot = 1.5` produce comparable internal mechanical behaviour?

5. Does `Cs = 1.0e-3` remain appropriate when `F2` is recomputed for the target grid?

6. Are differences in FIA associated primarily with:
   - resolved geometry,
   - lateral resistance,
   - internal rheology,
   - forcing resolution,
   - or EVP numerical convergence?

7. Can similar total FIA be demonstrated to arise from similar spatial and mechanical behaviour rather than compensating errors?

The most scientifically useful outcome is not necessarily that all three simulations produce identical FIA.

The useful result is establishing **why** FIA changes with resolution and whether the current lateral-drag and rheological parameters represent physical quantities that remain transferable across model scales.

---

# Relevant model implementation

The present CICE branch already provides the necessary conceptual separation.

In:

```text
cicecore/cicedyn/dynamics/ice_dyn_shared.F90
```

the tensile-strength parameter is defined through:

```text
T = Ktens × P
```

and the lateral-drag implementation contains separate:

```text
Cs
F2
```

components.

In:

```text
cicecore/cicedyn/infrastructure/ice_grid.F90
```

the model reads target-grid:

```text
F2x
F2y
```

fields and maps them to:

```text
F2E
F2N
```

on C-grid faces.

This separation is exactly what makes the proposed resolution experiment scientifically clean: **the constitutive sea-ice parameters can initially be held fixed while the geometric operator is reconstructed for each grid.**

---

# References

König Beatty, C. & Holland, D. M. (2010). *Modeling Landfast Sea Ice by Adding Tensile Strength*. Journal of Physical Oceanography, 40, 185–198. doi:10.1175/2009JPO4105.1

Liu, Y., Losch, M., Hutter, N. & Mu, L. (2022). *A New Parameterization of Coastal Drag to Simulate Landfast Ice in Deep Marginal Seas in the Arctic*. Journal of Geophysical Research: Oceans, 127, e2022JC018413. doi:10.1029/2022JC018413

Koldunov, N. V., Danilov, S., Sidorenko, D., Hutter, N., Losch, M., Goessling, H. F., Rakowsky, N., Scholz, P. & Jung, T. (2019). *Fast EVP Solutions in a High-Resolution Sea Ice Model*. Journal of Advances in Modeling Earth Systems. doi:10.1029/2018MS001485

Xu, S., Ma, J., Zhou, L., Zhang, Y., Liu, J. & Wang, B. (2021). *Comparison of sea ice kinematics at different resolutions modeled with a grid hierarchy in the Community Earth System Model*. Geoscientific Model Development, 14, 603–628. doi:10.5194/gmd-14-603-2021
