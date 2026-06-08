# ERA5 precipitation phase and Antarctic coastal snowfall forcing

## Purpose

This project tests whether Antarctic fast-ice thermodynamics are sensitive to the treatment of precipitation phase and coastal snowfall amount in standalone CICE6. The motivation is narrow: CICE should receive snowfall and rainfall as physically distinct atmospheric inputs when those fields are available from ERA5, rather than forcing all precipitation through a single total-precipitation pathway.

This is not a new tuning pathway for the lateral-drag model. It is a controlled forcing sensitivity around the free-slip/lateral-drag configuration.

## Scientific rationale

Near the Antarctic coast, 2 m air temperature can sit close to the freezing point. Small atmospheric temperature biases can therefore misclassify precipitation phase if CICE is required to infer rainfall versus snowfall internally. Separating ERA5 precipitation into `snowfall` and `rainfall` allows direct tests of:

- snow loading on sea ice;
- insulating effects of snow on basal growth;
- snow-albedo feedbacks;
- snow-ice formation pathways;
- regional coastal sensitivity in cells geometrically favourable to landfast ice.

The branch should document the implementation and numerical sanity of the forcing path, but should avoid publishing paper-level judgements about FIA skill or growth-rate impacts.

## Forcing-product inputs

Monthly ERA5 files are produced by `shuga` and should follow:

```text
${atm_data_dir}/era5_for_cice6_YYYY_MM.nc
```

The precipitation-relevant variables are:

| Variable | Meaning | Units | Use |
|---|---|---:|---|
| `ttlpcp` | total precipitation rate | kg m-2 s-1 | fallback / conservation check |
| `snowfall` | snowfall rate | kg m-2 s-1 | preferred `fsnow` forcing |
| `rainfall` | rainfall rate | kg m-2 s-1 | preferred `frain` forcing |

Core atmospheric fields needed by the same reader:

```text
airtmp, spchmd, glbrad, dlwsfc, wndewd, wndnwd
```

## Stage P0: monthly ERA5 backbone

### Goal

Verify that CICE can read the monthly ERA5 files and reproduce the previous total-precipitation semantics.

### Initial mapping

```fortran
Tair  <- airtmp
Qa    <- spchmd
fsw   <- glbrad
flw   <- dlwsfc
uatm  <- wndewd
vatm  <- wndnwd
fsnow <- ttlpcp
frain <- 0
```

### Filename logic

```fortran
subroutine ERA5_files_monthly(yr, mon)
```

Responsibilities:

- build `${atm_data_dir}/era5_for_cice6_YYYY_MM.nc`;
- respect `fyear`, `fyear_init`, `fyear_final`, and `ycycle`;
- keep file naming separate from physics;
- print the resolved file when `debug_forcing` is enabled.

### Record logic

For hourly monthly files:

```fortran
maxrec = 24 * daymo(mmonth)
recnum = 24 * (mday - 1) + int(real(msec,kind=dbl_kind) / c3600) + 1
```

Fluxes can initially use the current hourly record. State variables can either persist within the hour or interpolate between `recnum` and `recnum+1`, with month-end persistence used as the first simple boundary treatment.

## Stage P1: phase-separated precipitation

### Goal

Route ERA5 precipitation phase directly into CICE fields.

### Mapping

```fortran
fsnow <- snowfall
frain <- rainfall
```

### Sanity checks

- `snowfall >= 0`
- `rainfall >= 0`
- `ttlpcp >= 0`
- `abs(ttlpcp - (snowfall + rainfall))` is small enough for ERA5 accumulation/disaggregation conventions
- no NaN/Inf in `fsnow` or `frain`

## Stage P2: Antarctic coastal/form-factor snowfall scaling

### Goal

Allow controlled perturbation of snowfall in cells that are both Antarctic coastal-ocean cells and geometrically favourable to landfast ice.

### Control variables

Use existing perturbation controls where possible:

```fortran
era5_mod_var = 'snow_ant_lfi'
era5_mod_fac = 0.5
```

Recommended accepted aliases:

```text
snow_ant_coast
snow_ant_lfi
snow_f2
```

### Mask definition

The first crude coastal mask can use CICE grid information:

```text
ocean cell:      hm > 0.5
Antarctic band:  TLAT <= -50 deg
coastal cell:    at least one cardinal T-grid neighbour is land/boundary
```

The more fast-ice-centric mask should combine coastline proximity with form-factor magnitude:

```text
F2mag = sqrt(F2E^2 + F2N^2)
LFI-favourable if F2mag >= F2_threshold
```

A conservative first threshold should be diagnostic, not tuned. Suggested implementation:

```fortran
if (is_ant_coastal_ocean_cell(i,j,iblk) .and. F2mag(i,j,iblk) >= F2_min) then
   snowfall(i,j,iblk) = era5_mod_fac * snowfall(i,j,iblk)
endif
```

Apply the factor after the field is read and before the final mapping to `fsnow`.

## Required diagnostics

Do not add many history fields by default. Prefer master-task/global min/max and optional debug output:

```text
min/max/mean fsnow before scaling
min/max/mean fsnow after scaling
total snow mass before/after scaling
number/area of scaled cells
minimum/maximum F2mag in scaled cells
```

Optional history variables, off by default:

```text
fsnow_raw
fsnow_scaled
snow_scale_mask
F2mag_forcing_mask
```

## Success criteria

- Monthly ERA5 reader compiles and runs.
- Phase-separated precipitation path produces non-negative `fsnow` and `frain`.
- Total precipitation sanity checks are documented.
- Snowfall scaling affects only intended Antarctic coastal/form-factor cells.
- The repository documents that the branch has produced numerically reasonable output, without interpreting paper-level fast-ice skill.
