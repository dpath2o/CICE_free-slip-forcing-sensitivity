# Dynamic ORAS-derived mixed-layer depth in CICE

## Purpose

This project tests whether Antarctic fast-ice growth, persistence and seasonal timing are sensitive to realistic spatial and temporal variation in upper-ocean mixed-layer depth (MLD). The first implementation should remain deliberately simple: prescribe `hmix(i,j,t)` from ORAS-derived fields and pass it through the existing CICE → Icepack mixed-layer pathway without rewriting Icepack thermodynamics.

The central distinction is important:

> **The existing Icepack ocean mixed-layer scheme can accept a different `hmix` on successive calls, but it is a slab heat-capacity scheme, not a prognostic entraining mixed-layer model.**

Changing `hmix` is therefore possible and scientifically interpretable, provided the first experiment is described correctly: it changes the represented ocean thermal mass beneath the ice. It does **not** explicitly account for the heat and salt associated with entrainment or detrainment when the prescribed mixed-layer base moves.

This document first establishes exactly where the scheme resides in Icepack, where `hmix` resides in CICE, how the two are connected, and what occurs when `hmix` changes between timesteps. The ORAS experiment design follows from that conceptual analysis.

---

## 1. Where the mixed-layer scheme actually lives

The CICE branch used here is linked to the Icepack submodule at commit:

```text
Icepack commit: a5b5ebe63f986dbda86d5b2ef91426811619d018
```

The relevant implementation is split cleanly between CICE and Icepack:

- **CICE owns the gridded field** `hmix(i,j,iblk)` and supplies forcing.

- **CICE calls Icepack** once for each thermodynamic grid cell through `ocean_mixed_layer()`.

- **Icepack consumes the instantaneous `hmix` value** and uses it in the slab heat balance and freeze/melt-potential calculation.
- Icepack does not carry a separate previous mixed-layer depth and does not calculate `dhmix/dt`.

![CICE to Icepack mixed-layer call path](figures/cice-icepack-mixed-layer-path.svg)

### 1.1 Icepack: `icepack_ocn_mixed_layer`

The mixed-layer heat-balance routine is:

```text
icepack/columnphysics/icepack_ocean.F90
    module icepack_ocean
        subroutine icepack_ocn_mixed_layer(...)
```

In this routine, `hmix` is an input argument:

```fortran
real (kind=dbl_kind), intent(in) :: &
   ...
   hmix      , & ! mixed layer depth (m)
   ...
   dt            ! time step (s)
```

That is the first key point: `hmix` is **not** an Icepack prognostic state variable in this scheme. Icepack receives the value supplied by CICE for the current call. The surface/ice/deep-ocean heat fluxes change **SST** according to:

```fortran
sst = sst + dt * ( (fsens_ocn + flat_ocn + flwout_ocn + flw + swabs) * (c1-aice) + fhocn + fswthru) / (cprho*hmix)
sst = sst - qdp*dt/(cprho*hmix)
```

Schematically,

$$
\Delta T = \frac{F_{\mathrm{net}}\,\Delta t}{\rho c_p H},
$$

where

- $H = \mathrm{hmix}$,
- $\rho c_p H$ is the heat capacity per unit horizontal area of the represented ocean slab,
- $F_{\mathrm{net}}$ is the net combination of atmospheric, ice-ocean, penetrating-shortwave and deep-ocean heat fluxes using the CICE/Icepack sign conventions.

Ok, so I then think that for the same heat flux:

- a **shallower** mixed layer changes temperature more rapidly;
- a **deeper** mixed layer changes temperature more slowly.

Following that, Icepack then computes the energy available for freezing or melting:

```fortran
frzmlt = (Tf-sst)*cprho*hmix/dt
frzmlt = min(max(frzmlt,-frzmlt_max),frzmlt_max)
if (sst <= Tf) sst = Tf
```

The same current `hmix` is therefore used both to calculate the **SST* tendency and to convert any supercooling relative to `Tf` into a freezing potential.

### 1.2 CICE: `hmix` is a full T-grid field

CICE declares `hmix_0` as the initial/fixed standalone depth and `hmix` as a three-dimensional block field in:

```text
cicecore/cicedyn/general/ice_flux.F90
```

In that file:

```fortran
real(kind=dbl_kind), public :: &
   hmix_0         ! initial (or fixed if standalone) mixed layer depth
real(kind=dbl_kind), dimension(:,:,:), allocatable, public :: &
   ...
   hmix    , &     ! mixed layer depth (m)
   ...
```

My version of the standalone initialization sets:

```fortran
hmix(:,:,:) = hmix_0
```

### 1.3 CICE to Icepack (and then back to CICE): the instantaneous grid-cell value is passed on every call

The CICE wrapper is:

```text
cicecore/cicedyn/general/ice_step_mod.F90
    subroutine ocean_mixed_layer(dt, iblk)
```

For each thermodynamic grid cell it calls Icepack with the current value:

```fortran
call icepack_ocn_mixed_layer( &
     ...
     hmix  = hmix(i,j,iblk), &
     Tf    = Tf  (i,j,iblk), &
     qdp   = qdp (i,j,iblk), &
     frzmlt= frzmlt(i,j,iblk), &
     dt    = dt)
```

There is no requirement in that interface that the value supplied at timestep $n+1$ be equal to the value supplied at timestep $n$.

### 1.4 Time-varying MLD is already represented in the CICE forcing architecture

The existing *NCAR-style* ocean forcing pathway in:

```text
cicecore/cicedyn/general/ice_forcing.F90
```

includes:

```fortran
data vname / 'T', 'S', 'hblt', 'U', 'V', 'dhdx', 'dhdy', 'qdp' /
```

and maps the boundary/mixed-layer depth field (`hblt`) into CICE `hmix`. This is useful about the intended architecture: a spatially and temporally varying mixed-layer depth is not foreign to the CICE to Icepack interface. However, the legacy NCAR pathway imposes a lower bound equivalent to the default mixed-layer depth:

```fortran
hmix(i,j,iblk) = max(mixed_layer_depth_default, work1(i,j,iblk))
```

with `mixed_layer_depth_default = 60 m` in my experiments thus. **That 60-m floor should not be copied into the ORAS Antarctic experiment.** It would remove much of the shallow coastal and summer MLD variability that this experiment is intended to test. The ORAS pathway *might* instead use an explicit, scientifically chosen `hmix_min`, `hmix_max`, and bathymetry-aware limit.

---

## 2. What happens if `hmix` changes suddenly between timesteps?

The thought experiment is deliberately severe. Suppose CICE uses $H_1$ at timestep $n$, then the forcing system changes `hmix` to $H_2$ before timestep $n+1$.

![Conceptual step change in prescribed mixed-layer depth](figures/mixed-layer-step-change.svg)


### 2.1 Timestep \(n\): Icepack uses \(H_1\)

At timestep $n$, Icepack receives:

```text
sst  = T_n
hmix = H_1
```

and applies the heat fluxes using heat capacity $C_1 = \rho c_p H_1$. For a net heat input $F_{\mathrm{net}}$,

$$
\Delta T_1 = \frac{F_{net}\Delta t}{\rho c_p H_1}.
$$

If cooling lowers SST below the freezing point, the same $H_1$ is used in the calculation of `frzmlt` before SST is reset to `Tf`.

### 2.2 Between timesteps: CICE changes $H_1 \rightarrow H_2$

If external forcing changes:

```text
hmix: H_1 -> H_2
```

CICE does **not** automatically apply a compensating heat flux associated with that change in slab volume. The SST state is simply carried forward in the normal way. At the next Icepack call the model *sees*:

```text
sst  = T_n+1 carried from the previous model state
hmix = H_2 supplied by CICE forcing
```

There is no hidden Icepack operation of the form:

```text
new MLD -> entrain water -> mix temperature -> conserve mixed-layer heat content
```

and no corresponding *detrainment* operation if the layer becomes shallower. I think this distinction is pretty important as it does not inherrently translate to some numerical instability implied simply by changing `hmix`. Neither is there a hidden physical mechanism that accounts for the energetic consequence of moving the mixed-layer base.

### 2.3 Timestep $n+1$: Icepack uses $H_2$

At the next call the heat capacity is simply $C_2 = \rho c_p H_2$, so subsequent forcing produces

$$
\Delta T_2 = \frac{F_{\mathrm{net}}\Delta t}{\rho c_p H_2}.
$$

If $H_2 > H_1$, the model has instantaneously changed to a larger represented thermal reservoir and subsequent SST tendencies are smaller for a given flux. If $H_2 < H_1$, the represented thermal reservoir is smaller and subsequent SST tendencies are larger. The important physical question is what is omitted at the instant the represented slab depth changes.

---

## 3. Energy interpretation of a changing prescribed slab

### 3.1 The existing Icepack slab equation

Ignoring sign-convention detail and collecting the explicit fluxes into $F_{\mathrm{net}}$, the existing scheme is approximately:

$$
H\frac{dT}{dt} = \frac{F_{\mathrm{net}}}{\rho c_p}.
$$

The current value of $H$ changes the thermal inertia, but Icepack contains no explicit term involving $dH/dt$.

### 3.2 Apparent heat-content change when the prescribed depth changes

A useful diagnostic is the mixed-layer sensible-heat anomaly relative to the local freezing temperature:

$$
E'_{\mathrm{ml}} = \rho c_p H(T - T_f).
$$

If `hmix` changes instantaneously from $H_1$ to $H_2$ while temperature is carried forward, the represented heat anomaly changes by

$$
\Delta E'_{\mathrm{geom}} = \rho c_p (T - T_f)(H_2 - H_1).
$$

In the current slab scheme, this change is **not** balanced by an explicit model heat flux. I suppose the two qualifications that matter are:
1. If $T \approx T_f$, then $E'_{\mathrm{ml}}$ and this geometric heat-anomaly jump can be small even for a substantial change in depth.
2. A deepening real mixed layer can nevertheless entrain water whose temperature differs from the existing slab, so the missing entrainment heat flux can still matter even when the surface layer itself is near freezing.

### 3.3 What a physically variable-depth mixed layer would add

For a simple deepening slab, a more complete temperature equation contains an entrainment term of the form

$$
H\frac{dT}{dt} = \frac{F_{\mathrm{net}}}{\rho c_p} + (T_b - T)\frac{dH}{dt}, \qquad \frac{dH}{dt}>0,
$$

where $T_b$ is the temperature of water incorporated from immediately beneath the mixed layer. The corresponding entrainment heat flux scale is

$$ Q_{\mathrm{ent}} \approx \rho c_p (T_b - T_{\mathrm{ml}})\max\left(\frac{dH}{dt},0\right).
$$

The present Icepack slab does not calculate this term. A full prognostic mixed-layer treatment would also need buoyancy physics. Also shoaling requires its own interpretation rather than simply applying the deepening expression with the opposite sign. Therefore experiment M1 (below) probably should **not** be called anything close to a prognostic mixed-layer experiment. It is more precisely:

> **a prescribed variable-depth slab experiment that tests sensitivity to changing upper-ocean heat capacity, while omitting explicit entrainment/detrainment thermodynamics.**

### 3.4 Why M1 is sensible?

The omission above is a limitation, not an algebraic failure of the CICE/Icepack heat budget. Within each thermodynamic call:

- one value of `hmix` is used consistently;
- all explicit heat-flux SST tendencies use that value;
- `frzmlt` uses the same value;
- the resulting SST/freezing calculation is well defined.

The experiment becomes questionable only if the omitted entrainment/detrainment energy is comparable to, or larger than, the explicit fluxes that control the result. That is an empirical diagnostic question and can be quantified from ORAS temperature profiles before moving to a more complex model.

---

## 4. Interaction with SST restoring in this CICE branch

The *AFIM* ocean-forcing pathway restores SST by Newtonian temperature relaxation:

```fortran
sst = sst + (sst_interp - sst) * dt / trest
```

The **temperature restoring timescale** therefore does not itself depend on `hmix`:

$$
\frac{dT}{dt}\bigg|_{\mathrm{restore}} = \frac{T_{\mathrm{ORAS}} - T}{\tau}.
$$

However, the heat flux energetically equivalent to that temperature increment does depend on slab depth:

$$
Q_{\mathrm{restore}} = \rho c_p H\frac{T_{\mathrm{ORAS}} - T}{\tau}.
$$

This distinction should be kept explicit when interpreting experiments M1 to M4 (below):

- changing `hmix` does **not** change the specified restoring time $\tau$;
- changing `hmix` **does** change the amount of implied heat per unit area associated with a given restored temperature increment.

`Q_restore` should therefore be diagnosed whenever SST restoring is active.

---

## 5. Implications for fast-ice growth and seasonality

A variable `hmix` should not be assumed to translate monotonically into faster or slower mid-winter basal growth. If the slab is already at the local freezing temperature and a flux would supercool it, Icepack first produces a temperature tendency proportional to $1/H$, but then converts the resulting supercooling into `frzmlt` using a factor proportional to $H$. To first order those dependencies can cancel:

$$
\Delta T \propto \frac{1}{H},
\qquad
F_{\mathrm{frz}} \propto H\Delta T.
$$

The strongest MLD sensitivity may therefore appear in periods when the slab retains sensible heat above freezing, particularly:

- autumn cooling and freeze-up timing;
- shoulder-season thermal memory;
- spring warming and retreat timing;
- regions where subsurface heat is important;
- regions/times where SST restoring or `qdp` supplies substantial heat.

That is still directly relevant to Antarctic landfast-ice growth rate, persistence and seasonal timing, but it makes the expected mechanism more precise than simply assuming “shallower MLD = faster ice growth.”

---

## 6. Experiments

First I will determine whether the simple slab's response to realistic MLD variability is large enough to matter and whether the omitted entrainment heat term is first or second order.

| Stage | Name | Description | Interpretation |
|---|---|---|---|
| M0 | `fixed-control` | current fixed `hmix_0` and existing ocean pathway | baseline (completed with `hmix=20`) |
| M0S | `hmix-step-test` | short synthetic test with a deliberately abrupt local change, e.g. 20 to 80 m and reverse | verify exact numerical response and diagnostics; **not** a publishable experiment |
| M1 | `ORAS-hmix-14d` | ORAS-derived MLD low-pass filtered on ~14-day scales and continuously interpolated to model timesteps | variable slab heat capacity, no explicit entrainment |
| M2 | `ORAS-hmix-daily` | daily ORAS MLD, continuously interpolated | higher-frequency variable slab heat capacity |
| M3 | `ORAS-hmix-clim` | climatological seasonal-cycle MLD | isolate seasonal MLD structure from IAV/weather noise |
| M4 | `ORAS-hmix-sstsss` | dynamic MLD plus ORAS SST/SSS restoring | add thermodynamic ocean-state constraint |
| M5 | `ORAS-hmix-uocn` | add surface or mixed-layer-mean ORAS currents | add ocean-dynamic forcing pathway |
| M6 | `entraining/prognostic-ML` | add entrainment/detrainment or a prognostic mixed-layer scheme only if M1–M5 justify it | physically evolving mixed layer |

### Why M1 should be filtered/interpolated rather than updated as a staircase

The original concept of updating `hmix` every N days (seaonally or even as low as 14 days) is useful as a timescale, but there is little scientific benefit in imposing discrete 14-day jumps. A staircase creates artificial impulses in $dH/dt$ and maximizes the unrepresented geometric/entrainment energy change at each update.

The preferred M1 is therefore:

```text
ORAS MLD
  -> quality control / clipping
  -> ~14-day temporal low-pass or smoothing
  -> linear interpolation to every CICE thermodynamic timestep
  -> hmix(i,j,t)
```

M0S retains the deliberately sudden change as a transparent unit/physics test so that the response shown in the conceptual figure can be reproduced in the model.

---

## 7. ORAS preprocessing

Preprocess ORAS data offline in `shuga` and write CICE-ready files. Raw fields to amend to existing dataset. 

```text
mlotst_oras   mixed-layer thickness
thetao_oras   potential temperature with depth
so_oras       salinity with depth
uo_oras       eastward velocity with depth
vo_oras       northward velocity with depth
```

The temperature profile is particularly valuable even if M1 only forces `hmix`, because it allows an **offline estimate of the entrainment heat flux omitted by the simple slab scheme**. Do not download the full-depth global ocean unless a later diagnostic shows it is necessary. For the initial fast-ice-focused experiment, the upper 500 m should capture the MLD and most immediately relevant subsurface thermal structure in the target coastal regions.

---

## 8. CICE-ready implementation

Suggested file pattern:

```text
${ocn_data_dir}/oras_for_cice6_YYYY_MM.nc
```

Suggested variables:

| Variable | Meaning | Units | Initial use |
|---|---|---:|---|
| `hmix` | ORAS-derived mixed-layer depth | m | primary M1 forcing |
| `sst_oras` | ORAS near-surface temperature | °C or K, documented | diagnostics; later restoring |
| `sss_oras` | ORAS near-surface salinity | psu | diagnostics; later restoring / freezing temperature |
| `uocn_sfc` | ORAS surface eastward velocity | m s⁻¹ | later current forcing |
| `vocn_sfc` | ORAS surface northward velocity | m s⁻¹ | later current forcing |
| `uocn_mld_mean` | mixed-layer-mean eastward velocity | m s⁻¹ | later current forcing |
| `vocn_mld_mean` | mixed-layer-mean northward velocity | m s⁻¹ | later current forcing |
| `theta_mld_mean` | mixed-layer-mean potential temperature | °C or K | entrainment / OHC diagnostics |
| `theta_below_mld` | representative temperature immediately below MLD | °C or K | estimate `Q_ent` |
| `salt_mld_mean` | mixed-layer-mean salinity | psu | diagnostic |
| `ohc_mld` | upper-ocean heat content over MLD | J m⁻² | diagnostic |
| `profile_mld_sigma` | independently recomputed density MLD | m | quality check |

The initial CICE implementation only requires `hmix`. The additional fields are useful because they allow the limitations of M1 to be quantified rather than assumed.

---

## 9. Fortran implementation outline

### 9.1 Namelist controls

```fortran
logical :: use_dynamic_hmix
character(char_len) :: hmix_data_type     ! 'constant', 'oras', 'oras_clim'
character(char_len_long) :: hmix_data_dir
character(char_len_long) :: hmix_file_template
real(kind=dbl_kind) :: hmix_min
real(kind=dbl_kind) :: hmix_max
logical :: hmix_smooth_time
```

Defaults must preserve existing behaviour:

```fortran
use_dynamic_hmix = .false.
hmix_data_type   = 'constant'
hmix_0           = 60.0
```

### 9.2 Reader structure

Keep ORAS product knowledge in CICE forcing code or offline preprocessing, not in Icepack:

```text
ORAS preprocessing
    -> CICE forcing reader/interpolator
        -> hmix(i,j,iblk)
            -> CICE ocean_mixed_layer()
                -> Icepack icepack_ocn_mixed_layer()
```

A suitable reader pattern is:

```fortran
subroutine ORAS_files_monthly(yr, mon)
subroutine ORAS_hmix_data(dt)
```

with two records retained for temporal interpolation:

```fortran
hmix_data(nx_block,ny_block,2,max_blocks)
```

and mapped to the existing CICE field:

```fortran
hmix(:,:,:) = hmix_forcing(:,:,:)
```

### 9.3 Bounds and bathymetry

Use explicit controls:

```fortran
hmix = max(hmix_min, min(hmix, hmix_max))
hmix = min(hmix, local_bathymetry_safe_limit)
```

Do not silently inherit the legacy `max(60 m, hblt)` behaviour. The lower/upper bounds should be documented, tested, and accompanied by clipping counters.

---

## 10. Required diagnostics

The diagnostics should make it possible to distinguish a useful variable-heat-capacity sensitivity from artifacts caused by rapid prescribed depth changes.

### 10.1 MLD forcing diagnostics

```text
hmix
hmix_previous
DhmixDt = (hmix - hmix_previous) / dt
hmix min/max/mean by timestep
number of cells clipped at hmix_min
number of cells clipped at hmix_max
number of cells limited by bathymetry
regional hmix mean and percentiles by Antarctic sector
```

### 10.2 Mixed-layer thermodynamic diagnostics

```text
sst
Tf
sst_minus_Tf
frzmlt
qdp
fhocn
fswthru
net explicit mixed-layer heat flux used by Icepack
```
Calculate offline:

$$
E'_{\mathrm{ml}}=\rho c_p H(T - T_f)
$$

and the timestep change $\Delta E'_{\mathrm{ml}}$.

For M0S experiment, check explicitly that changing `hmix` alone does not create an Icepack flux impulse, while $E'_{\mathrm{ml}}$ changes geometrically when $T \ne T_f$.

### 10.3 Restoring-energy diagnostic

SST restoring diagnose

$$
Q_{\mathrm{restore}} = \rho c_p H\frac{T_{\mathrm{ORAS}} - T}{\tau}.
$$

This converts the temperature-restoring tendency into an energetically interpretable W m$^{-2}$ quantity.

### 10.4 Offline entrainment diagnostic

From ORAS MLD and temperature profiles estimate

$$
Q_{\mathrm{ent}} \approx \rho c_p (T_b - T_{\mathrm{ml}}) \max \left( \frac{dH}{dt},0 \right).
$$

Compare that magnitude and seasonal distribution with:

```text
|F_net|
|Q_restore|
|qdp|
```

particularly in cells with high fast-ice probability and during freeze-up/retreat periods. The purpose is not to force `Q_ent` in M1. It is to determine quantitatively whether omitting it is acceptable for the scientific question.

---

## 11. M0S: deliberately abrupt step test

Before a long M1 integration, perform a short controlled test in which one or a small set of ocean cells experiences a known MLD change while all other forcing is unchanged. Sequence:

```text
H = 20 m for several thermodynamic timesteps
H -> 80 m in one forcing update
hold H = 80 m
H -> 20 m
hold H = 20 m
```

Repeat under at least two SST states:

1. `sst ≈ Tf`, representing winter ice-covered conditions;
2. `sst > Tf`, representing a slab containing sensible heat.

Expected behaviour:
- `sst` itself should not jump solely because `hmix` changes;
- there should be no instantaneous compensating `fhocn`, `qdp`, or `frzmlt` pulse generated solely by the new depth;
- subsequent `sst` tendencies under the same explicit heat flux should scale approximately as `1/hmix`;
- `frzmlt` should use the new `hmix` consistently when SST reaches/surpasses the freezing constraint;
- the diagnosed `E'_ml` should expose the represented heat-content discontinuity when `sst != Tf`.

This test directly verifies the conceptual analysis above in the compiled CICE configuration.

---

## 12. Scientific interpretation of M1

The immediate hypothesis should be stated narrowly:

```text
Prescribing realistic ORAS-derived mixed-layer depth changes the effective
upper-ocean thermal inertia represented by the Icepack slab beneath Antarctic
coastal sea ice. This may alter freeze-up timing, persistence and retreat,
while explicit heat exchange associated with entrainment/detrainment across
the moving mixed-layer base remains omitted in M1.
```

If M1 is successful scientifically even if the fast-ice response is small. A weak response would indicate that the remaining seasonal timing/growth-rate bias is unlikely to be explained primarily by the slab heat-capacity assumption alone. A strong response should trigger the energy diagnostics above before it is interpreted as evidence for realistic mixed-layer physics.

---

## 14. Acceptance into paper criteria

### Software

- `use_dynamic_hmix = .false.` reproduces the fixed-control pathway bit-for-bit or practically equivalently.
- Dynamic `hmix` compiles and runs through the existing CICE → Icepack interface without Icepack source modification.
- `hmix` is finite, positive, bounded and bathymetry-aware in every wet cell.
- M0S reproduces the expected `1/hmix` thermal-inertia response.

### Energy diagnostics

- No unexplained instantaneous SST jump occurs when only `hmix` changes.
- `E'_ml`, `Q_restore`, and the offline `Q_ent` estimate are available for interpretation.
- The magnitude of omitted entrainment energy can be compared directly with the explicit mixed-layer heat fluxes.

### Science

- M1 changes only `hmix`; SST/SSS/u/v forcing remains unchanged relative to M0.
- M2–M5 add one forcing pathway at a time.
- A prognostic or entrainment-aware mixed layer is pursued only if the simpler experiments show that MLD physics materially affects the fast-ice result and that omitted entrainment is non-negligible.

---

## 15. Bottom line

There are two separate questions that should not be conflated:

1. **Can CICE/Icepack accept a mixed-layer depth that changes during a run?**  
   Yes. CICE stores `hmix` as a gridded field, existing forcing infrastructure already includes an `hblt -> hmix` pathway, and CICE passes the instantaneous grid-cell value into `icepack_ocn_mixed_layer()` on every call.

2. **Does the current Icepack slab fully conserve the heat content of a physically deepening/shoaling mixed layer?**
   No. The current scheme changes thermal capacity through the supplied `hmix`, but it does not contain an explicit entrainment/detrainment term associated with `dhmix/dt`.

That limitation does not make M1 nonsensical. It defines what M1 is: a controlled sensitivity to **prescribed upper-ocean thermal inertia**. The appropriate next step is to run that controlled experiment with explicit energy diagnostics, quantify the omitted entrainment term from ORAS profiles, and only then decide whether a more complete mixed-layer treatment is warranted.
