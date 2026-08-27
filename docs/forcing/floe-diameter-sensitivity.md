# Fixed floe-diameter sensitivity and lateral-melt geometry

## Purpose

This experiment tests the sensitivity of Antarctic fast ice to the fixed effective floe diameter used by Icepack in the parameterisation of lateral sea-ice melt.

The primary perturbation is:

```text
control:    floediam = 300 m
experiment: floediam = 100 m
```

with an optional smaller-floe bounding experiment:

```text
bounding:   floediam = 30 m
```

The experiment is run around the same standalone CICE6 free-slip/lateral-drag configuration used for the Paper 3 Antarctic fast-ice sensitivity experiments.

Strictly, this is **not an external forcing sensitivity**. It modifies an internal geometrical parameter controlling how available ocean heat is partitioned into lateral ice-area loss. It is nevertheless useful within the Paper 3 forcing-sensitivity programme because it tests whether an assumed sub-grid floe scale materially affects the seasonal evolution and persistence of Antarctic fast ice under otherwise identical atmospheric and oceanic forcing.

The experiment is intended primarily as a **bridge to the spectral-wave experiments**, where wave-induced breakup may make floe size prognostic or spatially variable. It also provides useful context for a later horizontal-resolution sensitivity suite, but should not be interpreted as a substitute for changing model resolution.

## Scientific rationale

The lateral melt scheme follows the geometrical argument used by Steele (1992), in which floe area and perimeter are related to a characteristic mean caliper diameter \(L\):

$$
S = \alpha L^2,
$$

with perimeter

$$
P = \pi L.
$$

Icepack retains this geometry using:

```fortran
floeshape = 0.66
floediam  = 300.0
```

where `floediam` is the effective floe diameter used for lateral melting.

The thermodynamic lateral melt calculation first computes a radial lateral melt rate from the ocean thermal forcing:

$$
w_{\mathrm{lat}} = m_1 \Delta T^{m_2},
$$

where

$$
\Delta T = \max(T_{\mathrm{SST}} - T_{\mathrm{bot}},0).
$$

The fractional ice-area loss associated with lateral melting is then calculated as

$$
r_{\mathrm{side}} = \frac{\pi w_{\mathrm{lat}}\Delta t} {\alpha L},
$$

which appears in Icepack as:

```fortran
wlat_loc = m1 * deltaT**m2
rside = wlat_loc * dt * pi / (floeshape * floediam)
```

Therefore, before energetic limiting,

$$
r_{\mathrm{side}} \propto \frac{1}{L}.
$$

Reducing `floediam` from 300 m to 100 m therefore increases the geometrical potential for fractional lateral area loss by a factor of three for the same ocean thermal forcing and model time step. his does **not** increase the radial melt velocity itself. Rather, the smaller assumed floes provide more ice perimeter per unit ice area. Icepack subsequently limits the combined basal and lateral melt demand by the available thermodynamic melt energy. The realised response therefore need not be exactly three times the lateral melt of the control, especially where the melt is strongly energy-limited.

## Fast-ice hypothesis

The experiment tests whether the fixed sub-grid floe geometry contributes materially to Antarctic fast-ice seasonality. The expected direct pathway is:

```text
smaller floediam
    ↓
greater floe perimeter per unit ice area
    ↓
larger potential fractional lateral melt
    ↓
faster concentration loss during melt conditions
    ↓
changes to ice mass, thickness distribution and strength
    ↓
changes to mobility and realised lateral-drag response
    ↓
possible changes to fast-ice retreat, persistence and subsequent freeze-up
```

The key point is that `floediam` acts directly during **lateral ablation**, not during winter lateral freezing. Any effect on fast-ice formation or growth timing is therefore indirect, through seasonal memory in concentration, open-water fraction, ocean heat content, thickness and dynamics. This makes the experiment particularly useful for determining whether an apparently incorrect fast-ice growth/retreat cycle could partly reflect an unrealistic fixed lateral-melt geometry rather than atmospheric or oceanic forcing alone.

## Relationship to the lateral-drag parameterisation

The coastal lateral-drag parameterisation should remain unchanged. The experiment should retain the same:

```text
boundary condition
lateral-drag switch
form function
coastline/form-factor fields
Cs / Cq / C_L parameters
rheology and strength parameters
```

used by the Paper 3 control. The fixed floe diameter is not a coefficient in the lateral-drag formulation itself. Any change in realised lateral drag should occur indirectly because altered lateral melt changes sea-ice concentration, thickness, mass, strength or mobility. This separation is scientifically useful: the experiment asks whether **thermodynamic floe geometry modifies the state on which the already-calibrated lateral-drag mechanism acts**.

## Experimental design

### F0: control

Use the existing Paper 3 free-slip/lateral-drag control:

```fortran
floediam = 300.0
```

No other physics or forcing changes.

### F1: smaller-floe experiment

Primary sensitivity:

```fortran
floediam = 100.0
```

All other namelist settings, forcing files, initial conditions, restart state and lateral-drag parameters remain identical to F0. The 100 m experiment is the preferred first perturbation because it gives a substantial but moderate change in effective floe perimeter:

```text
300 m → 100 m
potential lateral area-melt factor ≈ 3
```

before thermodynamic energy limiting.

### F2: small-floe bounding experiment

Optional bounding case:

```fortran
floediam = 30.0
```

This should be treated as a sensitivity bound rather than as a claim that 30 m is a representative Antarctic coastal floe diameter. Relative to the 300 m control:

```text
300 m → 30 m
potential lateral area-melt factor ≈ 10
```

before thermodynamic energy limiting. A weak fast-ice response even in this case would provide strong evidence that the fixed-floe lateral-melt assumption is not a first-order control on the simulated fast-ice seasonal cycle.


## Paper 3 diagnostics

The experiment should be evaluated using the same Paper 3 diagnostics as the external-forcing experiments so that the response is directly comparable. Primary fast-ice diagnostics:

```text
circumpolar FIA seasonal cycle
regional FIA seasonal cycle
FIP and FIP difference maps
fast-ice onset / retreat timing
fast-ice duration / persistence
```

Supporting sea-ice diagnostics:

```text
SIA / SIE
ice thickness
snow thickness
ice strength
ice speed
lateral melt
basal melt
ocean surface thermal forcing
```

Where possible, diagnose seasonal differences separately for:

```text
autumn freeze-up
winter growth / consolidation
spring transition
summer retreat
```

The most important mechanistic test is whether the response appears first during the melt season, as expected from the `floediam` formulation, and whether that anomaly then persists into the following fast-ice growth season.

## Interpretation constraints

This experiment should **not** be interpreted as a direct test of wave forcing. A smaller fixed `floediam` can mimic one thermodynamic consequence of having more fragmented ice — increased perimeter per unit area — but it does not represent:

```text
wave propagation
spectral attenuation
flexural fracture
event-driven breakup
directional swell
wave-induced stress
a prognostic floe-size distribution
```

Similarly, it should not be interpreted as a horizontal-resolution experiment. `floediam` is a sub-grid thermodynamic scale and is not set by the CICE grid spacing. The experiment instead answers the narrower question:

> How sensitive is the simulated Antarctic fast-ice seasonal cycle to the unresolved floe perimeter assumed by the standard lateral-melt parameterisation?

## Relationship to spectral-wave experiments

This is most naturally treated as a precursor to the spectral-wave pathway. The sequence is:

```text
fixed floe geometry
    ↓
floediam sensitivity
    ↓
diagnose whether floe perimeter matters
    ↓
wave-forcing diagnostics
    ↓
wave-induced breakup / FSD response
    ↓
active spectral wave–ice coupling
```

There are two useful outcomes.

### If the 100 m / 30 m experiments matter strongly

A strong response would show that Antarctic fast-ice seasonality is sensitive to floe perimeter and lateral area loss. This strengthens the motivation for replacing a spatially and temporally fixed floe diameter with a physically varying floe-size distribution influenced by wave breakup.

### If the experiments matter weakly

A weak response would also be valuable. It would suggest that simply increasing floe perimeter is unlikely to explain the fast-ice seasonal bias, and that the important wave pathway — if present — is more likely to involve mechanical breakup, mobility, strength reduction or loss of coastal attachment rather than thermodynamic lateral melt alone. The experiment therefore helps constrain the mechanism before implementing much more expensive spectral wave coupling.

## Relationship to horizontal model resolution

A later horizontal-resolution suite is scientifically distinct but complementary.

Suggested resolutions:

```text
coarse:       ~1 degree
baseline:     1/4 degree
eddy-rich:    ~1/10 degree
```

The purpose would be to test the sensitivity of Antarctic fast ice to the **resolved spatial scale** of:

```text
coastline geometry
narrow embayments and channels
bathymetric structure
coastal form factors
sea-ice deformation
ocean-current gradients
atmospheric coastal gradients
ice-edge geometry
```

This differs from the `floediam` experiment, which changes the assumed **sub-grid floe perimeter** while keeping the resolved grid unchanged. For a clean resolution study, first keep the physical floe parameter fixed:

```text
1 degree:    floediam = 300 m
1/4 degree:  floediam = 300 m
1/10 degree: floediam = 300 m
```

and regenerate the required grid-dependent forcing, coastline/form-factor and bathymetric fields consistently at each resolution. Only after establishing the pure resolution response should a cross-sensitivity be considered, for example:

| Grid | 300 m | 100 m |
|---|---:|---:|
| ~1° | control | optional |
| 1/4° | control | primary floe test |
| ~1/10° | control | optional |

This would distinguish:

```text
resolved-scale sensitivity
from
sub-grid floe-geometry sensitivity
```

without conflating the two.

## Recommended Paper 3 placement

Within the Paper 3 experiment hierarchy:

```text
FORCING / PROCESS SENSITIVITIES
│
├── atmospheric forcing
│   ├── precipitation / snowfall
│   └── boundary-layer forcing
│
├── ocean forcing
│   ├── mixed-layer depth / ocean heat capacity
│   └── tides
│
├── wave forcing
│   └── spectral wave boundary
│
└── process-scale bridge
    └── fixed floe-diameter / lateral-melt geometry
```

This makes the conceptual break explicit rather than pretending that `floediam` is an external forcing field. For the manuscript, it can be described as a **process sensitivity used to test the dependence of the forcing response on unresolved floe geometry**. A later model-resolution suite would sit outside this forcing hierarchy as a separate **structural/model-scale sensitivity**.

## Success criteria

- The 100 m experiment differs from the control only through `floediam`.
- Lateral-melt diagnostics respond in the expected direction during melt conditions.
- The model energy limiter remains numerically well behaved.
- Any FIA/FIP response can be traced through concentration, thickness, strength and mobility.
- Changes appear first in physically plausible melt-season conditions rather than as an unexplained winter-only response.
- The experiment provides a clear decision point for whether floe-size physics deserves active treatment in the spectral-wave experiments.
- The result is kept distinct from model-resolution sensitivity.

## Source-code context

The current Icepack implementation defines:

```fortran
floeshape = 0.66
floediam  = 300.0
```

and calculates lateral area loss using:

```fortran
wlat_loc = m1 * deltaT**m2
rside = wlat_loc * dt * pi / (floeshape * floediam)
```

The relevant implementation is in Icepack `icepack_therm_vertical.F90`, with the parameter definition in `icepack_parameters.F90`.

The CICE forcing-sensitivity branch currently pins the companion Icepack repository as a submodule, so the exact Icepack commit used by each experiment should be recorded with the model configuration.

## Reference

Steele, M. (1992). Sea ice melting and floe geometry in a simple ice-ocean model. *Journal of Geophysical Research*, **97**(C11), 17,729–17,738.
