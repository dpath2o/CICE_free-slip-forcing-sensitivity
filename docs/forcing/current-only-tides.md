# Current-only tidal forcing in standalone CICE

## Purpose

This project introduces barotropic tidal-current forcing into standalone CICE while preserving the free-slip/lateral-drag fast-ice configuration. The first implementation deliberately uses **current-only** tides. Tidal SSH-gradient forcing is kept disabled so that the ocean-current stress pathway can be isolated before adding a second momentum source.

## Scientific rationale

Antarctic landfast ice can sit close to a mechanical mobility threshold. A daily-mean ocean-current product can miss sub-daily tidal mobilisation, especially near coastlines, grounded icebergs, ice tongues, and shallow continental-shelf regions. Current-only tides provide a way to test whether repeated sub-daily ocean-stress perturbations can reduce over-persistence, alter breakout timing, or modify seasonal FIA maxima without changing the atmospheric or thermodynamic forcing.

## Existing implementation state

The current implementation is based on CATS2008_v2023 harmonic tidal currents remapped to the CICE grid. The test files and notes in `CICE_free-slip-tides` document a current-only pathway with:

```fortran
tide_data_type             = 'harmonic'
tide_data_format           = 'CICE_TMD3'
tide_use_currents          = .true.
tide_use_ssh               = .false.
tide_use_bathymetry_limit  = .true.
tide_curr_fac              = 1.0
tide_speed_cap             = 2.0
tide_wct_min               = 20.0
tide_wct_full              = 80.0
tide_h_eff_min             = 50.0
tide_depth_mismatch        = 4.0
tide_ramp_days             = 20.0
```

The current-only branch has already produced stable CICE history output in short and longer diagnostic integrations. The repository-ready documentation should focus on the implementation pathway and diagnostics, not on final paper-level skill.

## Tide forcing file contract

The CATS prototype file should contain harmonic constituents and bathymetry/validity fields on the CICE grid. The expected conceptual content is:

```text
hRe, hIm      optional SSH harmonic components
URe, UIm      eastward tidal-current harmonic components
VRe, VIm      northward tidal-current harmonic components
wct           CATS water-column thickness
cats_mask     CATS ocean/validity mask
omega         constituent angular frequency
phase         constituent phase metadata
alpha         nodal/constituent metadata where needed
amplitude     constituent metadata where needed
```

The current-only path requires `URe/UIm/VRe/VIm`, `omega`, and sufficient mask/depth information to avoid shallow-water blow-ups.

## Fortran pathway

### Initialisation

1. Read tide namelist controls.
2. Open the harmonic tide file.
3. Read constituent metadata and current components.
4. Read `wct` and `cats_mask`.
5. Build a bathymetry-aware limiter using CICE bathymetry and CATS water-column thickness.
6. Initialise a half-cosine ramp to avoid shocking the ice state at branch start.

### Runtime update

At each forcing timestep:

```fortran
utide = sum_constituents( URe*cos(theta) - UIm*sin(theta) )
vtide = sum_constituents( VRe*cos(theta) - VIm*sin(theta) )
```

Then apply:

```fortran
utide = tide_ramp * tide_curr_fac * limiter * utide
vtide = tide_ramp * tide_curr_fac * limiter * vtide
```

and add the tide perturbation to the background ocean current used by CICE:

```fortran
uocn = uocn_background + utide
vocn = vocn_background + vtide
```

The SSH-gradient contribution should remain disabled in this project:

```fortran
tide_use_ssh = .false.
```

SSH forcing should become a separate experiment only after the current-only pathway has been diagnosed in isolation.

## Bathymetry-aware limiter

The limiter should prevent conversion of transport-like tidal information into unrealistically large velocity in poorly resolved shallow cells.

Inputs:

```text
CATS wct
CATS mask
CICE bathymetry
CICE ocean mask
```

Controls:

```text
tide_wct_min
tide_wct_full
tide_h_eff_min
tide_depth_mismatch
tide_speed_cap
```

Recommended behaviour:

- zero tides where CATS or CICE says the cell is invalid;
- taper tides smoothly between `tide_wct_min` and `tide_wct_full`;
- use an effective minimum depth to avoid large velocity from small water-column thickness;
- detect large CATS/CICE bathymetry mismatch;
- retain a high emergency speed cap as a safety limit, not as the primary physics.

Diagnostic counters should be globally reduced across MPI tasks before being interpreted:

```text
n_valid_global
n_tapered_global
n_zeroed_global
n_mismatch_global
```

## Required tide diagnostics

The next model-side requirement is to make the sub-daily tide signal visible without requiring every downstream analysis to reconstruct it indirectly.

Recommended optional history fields:

```text
utide, vtide
tide_speed
tide_speed_rms_daily
tide_speed_max_daily
tide_exceed_count_daily
uocn_background, vocn_background
uocn_total, vocn_total
```

Recommended stress diagnostics:

```text
abs_tau_ocn_control
abs_delta_tau_ocn
preferred ratio: area_mean(abs_delta_tau_ocn) / area_mean(abs_tau_ocn_control)
near-threshold abs_delta_tau_ocn
```

The ratio of area means is preferred for headline force-balance interpretation. Area means of local pointwise ratios should be kept as a vulnerability diagnostic only, because they become large when the control denominator is small.

## Next experiment: year-long hourly-output tide analysis

The next major experiment should be a year-long current-only integration with hourly history output sufficient to diagnose whether sub-daily mobility survives the persistence classifier.

Minimum hourly fields:

```text
aice, hi, uvel, vvel, uocn, vocn,
strocnx, strocny, strairx, strairy,
strintx, strinty, strength, divu, shear,
utide, vtide, tide_speed
```

If full global hourly output is too expensive, use a dual-stream compromise:

```text
hourly: Antarctic diagnostic subset or minimal tide/velocity fields
daily: standard averaged ice history fields
```

## `shuga` analysis requirement: binary-days from hourly sea ice

The model should not need to classify fast ice internally at the first stage. Instead, `shuga` should gain an hourly-aware classifier.

Recommended workflow:

1. Load hourly `aice`, `uvel`, `vvel`, and optionally `tide_speed`.
2. Compute hourly ice speed:

```python
speed = sqrt(uvel**2 + vvel**2)
```

3. Compute an hourly raw fast-ice mask:

```python
FI_hourly = (aice >= aice_threshold) & (speed <= speed_threshold)
```

4. Collapse hourly masks to daily binary states using several documented options:

```text
majority-fast:    day is fast if >= 12 hours are fast
strict-fast:      day is fast if all valid hours are fast
tide-sensitive:   day is mobile if >= 1 or >= 3 hours exceed the speed threshold
```

5. Apply the existing binary-days persistence rule to the derived daily mask:

```text
window = 11 days
minimum fast days = 9
```

6. Compare the resulting `FI_mask` with the daily-output classifier to identify how much hourly tidal mobilisation is filtered by daily averaging.

## Decision point before SSH tides

Only add SSH-gradient forcing after the current-only experiment answers these questions:

- Does current-only tidal forcing change hourly or daily raw mobility in near-threshold cells?
- Does that mobility survive the binary-days persistence rule?
- Are changes regionally concentrated in tidally active shelf/coastal sectors?
- Does the effect appear as timing/persistence rather than only as instantaneous speed noise?

If yes, SSH forcing should be added as a separate, traceable experiment rather than blended into the current-only pathway.
