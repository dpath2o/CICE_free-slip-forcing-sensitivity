# CICE_free-slip-forcing-sensitivity

CICE free-slip lateral-drag forcing sensitivity experiments for Antarctic fast ice.

This repository tests whether the missing fast-ice-area (FIA) growth-rate response in the lateral-drag CICE experiments is primarily atmospheric, oceanic, or dynamical. The immediate focus is atmospheric forcing sensitivity using a revised monthly ERA5 forcing product produced by `shuga`.

The scientific control remains the lateral-drag/free-slip configuration. The forcing experiments should be interpreted as sensitivity tests around that configuration, not as a new tuning pathway that replaces the lateral-drag story.

## 1. Forcing-product contract

The new ERA5 forcing product is produced offline by `mawsons-chest/shuga` and written as monthly NetCDF files on the CICE grid.

Expected file pattern:

```text
${atm_data_dir}/era5_for_cice6_YYYY_MM.nc
```

Current example:

```text
/g/data/gv90/da1339/afim_input/ERA5/0p25/bilinear/monthly_cice6/era5_for_cice6_1994_10.nc
```

Each file contains hourly fields on the CICE T grid:

| Variable | Meaning | Units | CICE use |
|---|---|---:|---|
| `airtmp` | 2 m air temperature | K | direct atmospheric state |
| `spchmd` | 2 m specific humidity | kg kg-1 | direct atmospheric state |
| `pair` | surface pressure | Pa | future air-density / boundary-layer physics |
| `glbrad` | downward shortwave radiation | W m-2 | atmospheric flux |
| `dlwsfc` | downward longwave radiation | W m-2 | atmospheric flux |
| `ttlpcp` | total precipitation rate | kg m-2 s-1 | fallback total precipitation |
| `snowfall` | snowfall rate | kg m-2 s-1 | preferred snow input |
| `rainfall` | rainfall rate | kg m-2 s-1 | preferred rain input |
| `wndewd` | eastward 10 m wind | m s-1 | wind forcing |
| `wndnwd` | northward 10 m wind | m s-1 | wind forcing |
| `blh` | boundary-layer height | m | future optional boundary-layer sensitivity |
| `windgust` | 10 m wind gust | m s-1 | future optional wind-burst diagnostic/sensitivity |
| `wnd100ewd` | eastward 100 m wind | m s-1 | future optional boundary-layer / shear sensitivity |
| `wnd100nwd` | northward 100 m wind | m s-1 | future optional boundary-layer / shear sensitivity |

The first implementation step should only prove that CICE can read these monthly files and reproduce the existing ERA5-forced behaviour using the core fields. Perturbations should be layered on only after that backbone is compile- and run-tested.

## 2. Implementation principles

1. Preserve the existing lateral-drag and tide-capable code path.
2. Avoid Icepack changes for the atmospheric forcing experiments.
3. Keep the first Fortran change focused on file discovery and variable read-in only.
4. Keep perturbation controls namelist-driven through existing or minimally extended `forcing_nml` variables.
5. Keep diagnostic history additions optional and avoid them unless a field must be verified from CICE output.
6. Prefer monthly files to annual files to avoid very large forcing files and to make production/recovery easier.
7. Do not rely on zarr in Fortran; NetCDF remains the CICE interface.

## 3. Fortran backbone: monthly ERA5 read-in

### 3.1 Current baseline

The current `ice_forcing.F90` already contains:

- `F_ERA5` as the ERA5 forcing filename.
- `fsw_data`, `fsnow_data`, `Tair_data`, `uatm_data`, `vatm_data`, `Qa_data`, `rhoa_data`, `flw_data`, and `frain_data` as two-record atmospheric buffers.
- `era5_mod_var` and `era5_mod_fac` as existing ERA5 perturbation controls.
- `hmix_0` and ocean mixed-layer infrastructure for later ocean-forcing work.

The first modification should extend this structure rather than replacing it.

### 3.2 New filename logic

Add or revise the ERA5 filename helper so that CICE opens monthly files:

```text
${atm_data_dir}/era5_for_cice6_YYYY_MM.nc
```

where `YYYY` is the forcing-cycle-adjusted year and `MM` is the model month.

Recommended helper concept:

```fortran
subroutine ERA5_files_monthly(yr, mon)
```

Responsibilities:

- build `F_ERA5`;
- use `fyear` / `fyear_init` / `ycycle` consistently with existing forcing-cycle logic;
- print the resolved file on the master task when `debug_forcing` is enabled;
- keep the monthly filename logic isolated from the physics.

### 3.3 Monthly record logic

For hourly monthly files:

```text
maxrec = 24 * daymo(mmonth)
recnum = 24 * (mday - 1) + int(real(msec,kind=dbl_kind) / c3600) + 1
```

Initial implementation:

- read record `recnum`;
- for state variables, optionally read `recnum+1` within the same month and interpolate later;
- for flux variables, use the current hourly record without interpolation;
- handle month-end simply first by persistence at the final record.

A more complete second pass can implement a next-month look-ahead at month boundaries. The first goal is compile/runtime success for September-December 1994.

### 3.4 Variables to read in the backbone

The first backbone commit should read:

```text
airtmp
spchmd
glbrad
dlwsfc
ttlpcp
snowfall
rainfall
wndewd
wndnwd
```

It may also read but not yet use:

```text
pair
blh
windgust
wnd100ewd
wnd100nwd
```

The first run should still populate CICE fields as:

```fortran
Tair  <- airtmp
Qa    <- spchmd
fsw   <- glbrad
flw   <- dlwsfc
uatm  <- wndewd
vatm  <- wndnwd
fsnow <- ttlpcp
```

That preserves the existing total-precipitation behaviour before the rain/snow split is activated.

## 4. Stage 1 experiment: monthly ERA5 backbone

### Case name

```text
ERA5-monthly-core
```

### Purpose

Verify that CICE can read the new monthly ERA5 files and reproduce the current behaviour using the existing forcing semantics.

### Expected changes

- Monthly filename helper.
- Monthly record indexing.
- Read core variables.
- Keep total precipitation mapped to `fsnow` initially, matching the older ERA5 path.
- No atmospheric perturbation yet.
- No boundary-layer physics yet.

### Success criteria

- Compile succeeds.
- CICE runs through September-December 1994.
- Atmospheric forcing min/max diagnostics look physically sane.
- No NaNs/Inf in `Tair`, `Qa`, `fsw`, `flw`, `fsnow`, `uatm`, `vatm`.
- Fast-ice diagnostics are not yet interpreted scientifically; this is an IO/control test.

## 5. Stage 2 experiment: separate rainfall and snowfall

### Case name

```text
ERA5-monthly-phase
```

### Purpose

Use ERA5 precipitation phase information instead of relying on CICE to infer phase from near-surface air temperature alone.

### Expected changes

Read:

```text
snowfall
rainfall
```

Then map:

```fortran
fsnow <- snowfall
frain <- rainfall
```

rather than:

```fortran
fsnow <- ttlpcp
```

If CICE atmospheric forcing currently lacks a direct `frain` use in this branch, add the least invasive plumbing needed to pass rain separately into the existing flux pathway.

### Science rationale

Near coastal Antarctica, temperature can hover near freezing. Small temperature biases can misclassify precipitation phase. ERA5-provided snowfall and rainfall allow a more physically direct test of snow loading, insulation, albedo, snow-ice formation, and growth-rate impacts.

### Success criteria

- Total precipitation sanity check: `ttlpcp approximately equals snowfall + rainfall`.
- `fsnow` and `frain` remain non-negative.
- FIA response can be compared against the total-precipitation backbone.

## 6. Stage 3 experiment: Antarctic coastal snowfall scaling

### Case names

```text
ERA5-snow50coast
ERA5-snow75coast
ERA5-snow125coast
```

### Purpose

Test whether fast-ice growth-rate biases are sensitive to coastal snowfall amount.

### Expected control

Use existing ERA5 perturbation controls where possible:

```fortran
era5_mod_var = 'snow_ant_coast'
era5_mod_fac = 0.5
```

### Coastal mask

Define a deliberately blunt Antarctic coastal-ocean selector:

- T-grid ocean cell: `hm > 0.5`.
- Antarctic/Southern Ocean latitude limit, for example `TLAT <= -50 deg`.
- Adjacent to land in at least one of the four cardinal T-grid neighbours.

This is a sensitivity experiment, not a production-quality coastal downscaling product.

### Application point

Apply the multiplicative factor after reading/interpolating `snowfall`, before mapping it to `fsnow` and before the final land mask.

## 7. Stage 4 experiment: wind perturbations

### Case names

```text
ERA5-wind10coast
ERA5-windgust-diagnostic
ERA5-wind100blend
```

### Purpose

Assess whether fast-ice formation/growth is sensitive to wind forcing amplitude or boundary-layer wind structure.

### Candidate perturbations

#### 7.1 Coastal 10 m wind scaling

```fortran
era5_mod_var = 'wind_ant_coast'
era5_mod_fac = 1.1
```

Apply to both `uatm` and `vatm` over the same Antarctic coastal-ocean mask.

#### 7.2 Wind gust diagnostic

Read `windgust` but do not initially force CICE with it.

Possible diagnostic:

```text
gust_factor = windgust / sqrt(wndewd^2 + wndnwd^2)
```

This is useful for identifying wind burstiness but should not be used directly as vector wind forcing.

#### 7.3 100 m wind / shear sensitivity

Read `wnd100ewd` and `wnd100nwd`.

Possible future experiments:

- use 100 m winds directly;
- blend 10 m and 100 m winds;
- diagnose low-level shear.

Do not implement these before the basic monthly reader, phase separation, and snowfall-scaling cases are validated.

## 8. Stage 5 experiment: boundary-layer-aware forcing

### Case names

```text
ERA5-rhoa-pair
ERA5-blh-diagnostic
ERA5-blh-stability-proxy
```

### Purpose

Introduce optional boundary-layer information without immediately changing core CICE thermodynamics.

### 8.1 Surface pressure / air density

Read `pair` and compute a more variable air density:

```fortran
rhoa = pair / (Rdry * Tair)
```

or, later:

```fortran
rhoa = pair / (Rdry * T_virtual)
```

where virtual temperature depends on specific humidity.

This is a low-risk improvement because wind stress and turbulent fluxes are sensitive to air density.

### 8.2 Boundary-layer height

Read `blh` as a diagnostic field first.

Potential later uses:

- modulate wind coupling under shallow stable boundary layers;
- diagnose katabatic/coastal boundary-layer regimes;
- define regimes for future parameterised wind perturbations.

This should remain optional and off by default.

## 9. Suggested commit sequence

Use small commits so each stage compiles independently.

```text
1. forcing: document ERA5 monthly forcing experiment plan
2. forcing: add monthly ERA5 filename helper
3. forcing: read monthly ERA5 core variables
4. forcing: read optional ERA5 boundary-layer variables
5. forcing: separate ERA5 rainfall and snowfall
6. forcing: add Antarctic coastal mask helper
7. forcing: add coastal snowfall scaling
8. forcing: add coastal wind scaling
9. forcing: add optional surface-pressure air-density path
```

Each commit should be compile-tested before moving to the next stage.

## 10. Initial run plan

Available monthly files:

```text
1994-09
1994-10
1994-11
1994-12
```

Recommended first CICE smoke test:

```text
start date: 1994-09-01
end date:   1994-12-31
```

Recommended order:

1. `ERA5-monthly-core`
2. `ERA5-monthly-phase`
3. `ERA5-snow50coast`
4. `ERA5-wind10coast`
5. `ERA5-rhoa-pair`

Do not interpret the science until the `ERA5-monthly-core` and `ERA5-monthly-phase` cases are confirmed numerically sane.

## 11. Notes for future ocean-forcing experiments

The atmospheric experiments should be followed by an ocean mixed-layer experiment, not mixed into the first atmospheric branch.

Future ocean path:

- monthly climatological mixed-layer depth;
- bathymetry-aware mixed-layer-depth bounds;
- ORAS-derived daily mixed-layer depth;
- ORAS-derived mixed-layer-averaged temperature, salinity, and velocity;
- eventually fully coupled or ROMS-coupled experiments.

Those belong in a later branch once the atmospheric forcing interface is stable.
