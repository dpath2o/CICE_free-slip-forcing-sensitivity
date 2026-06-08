# ERA5 atmospheric boundary-layer forcing

## Purpose

This project explores whether Antarctic fast-ice seasonality is sensitive to near-surface wind structure and atmospheric boundary-layer state, beyond the standard 10 m wind and fixed/diagnostic air-density pathway normally used by CICE.

This is scientifically distinct from the precipitation-phase project. Precipitation is primarily a thermodynamic/snow-loading test. Boundary-layer forcing is primarily a momentum, burstiness, and seasonality test.

## ERA5 variables already available in the monthly forcing product

| Variable | Meaning | Units | Initial use |
|---|---|---:|---|
| `pair` | surface pressure | Pa | air-density diagnostic / optional forcing |
| `blh` | boundary-layer height | m | diagnostic regime variable |
| `windgust` | 10 m wind gust | m s-1 | burstiness diagnostic |
| `wnd100ewd` | eastward 100 m wind | m s-1 | shear / alternative wind diagnostic |
| `wnd100nwd` | northward 100 m wind | m s-1 | shear / alternative wind diagnostic |
| `wndewd` | eastward 10 m wind | m s-1 | current wind forcing |
| `wndnwd` | northward 10 m wind | m s-1 | current wind forcing |

## Stage B0: read but do not perturb

The first stage is diagnostic only. CICE should read the additional fields and verify that they are physically sane, but the model should continue to use the existing 10 m winds and existing turbulent-flux path.

Recommended diagnostics:

```text
wind10_speed = sqrt(wndewd^2 + wndnwd^2)
wind100_speed = sqrt(wnd100ewd^2 + wnd100nwd^2)
gust_factor = windgust / max(wind10_speed, eps)
low_level_shear = wind100_speed - wind10_speed
pair_min_max_mean
blh_min_max_mean
```

## Stage B1: surface-pressure air density

### Rationale

Wind stress and turbulent fluxes depend on air density. A fixed or overly simple air-density field may under-represent synoptic and coastal pressure effects.

### Implementation

Start with dry air density:

```fortran
rhoa = pair / (Rdry * Tair)
```

Then optionally use virtual temperature:

```fortran
Tv   = Tair * (1.0 + 0.61 * Qa)
rhoa = pair / (Rdry * Tv)
```

Keep this off by default:

```fortran
use_era5_pair_rhoa = .false.
```

or via an existing perturbation string:

```fortran
era5_mod_var = 'rhoa_pair'
```

## Stage B2: coastal 10 m wind scaling

### Rationale

If fast-ice seasonality is controlled by marginal mobility and episodic breakout, modest coastal wind-stress perturbations may be more relevant than seasonal-mean thermodynamic perturbations.

### Implementation

Use the same Antarctic coastal/form-factor mask as the snowfall experiment, but apply a vector-preserving scale:

```fortran
uatm = era5_mod_fac * uatm
vatm = era5_mod_fac * vatm
```

only inside the selected coastal/LFI-favourable cells.

Candidate control:

```fortran
era5_mod_var = 'wind_ant_lfi'
era5_mod_fac = 1.1
```

## Stage B3: gust and wind-burst diagnostics

`windgust` should not be used directly as a vector wind forcing because it lacks direction. It is more useful as a burstiness diagnostic or a scalar regime indicator.

Recommended first diagnostic:

```text
gust_factor = windgust / max(wind10_speed, eps)
```

Potential future perturbation, if justified:

```fortran
wind_speed_eff = wind10_speed * (1 + alpha_gust * max(gust_factor - 1, 0))
```

with the original 10 m wind direction retained. This should remain conceptual until the diagnostic distribution is understood.

## Stage B4: 100 m wind / shear sensitivity

The 100 m winds can be used in three escalating ways:

1. diagnostic only: compare 100 m and 10 m speed/direction;
2. blended forcing: `u_eff = (1 - w) * u10 + w * u100`;
3. regime dependent: increase `w` under shallow/stable boundary-layer conditions inferred from `blh`.

Candidate blend:

```fortran
u_eff = (1.0 - alpha_100) * wndewd  + alpha_100 * wnd100ewd
v_eff = (1.0 - alpha_100) * wndnwd  + alpha_100 * wnd100nwd
```

Keep `alpha_100 = 0` by default.

## Stage B5: boundary-layer height regime diagnostics

`blh` should first be used to identify coastal regimes, not to directly scale winds.

Potential diagnostics:

```text
mean wind speed where blh < 100 m
mean gust_factor where blh < 100 m
frequency of shallow-boundary-layer states by Antarctic sector
relationship between shallow BLH and raw fast-ice mobility
```

Only after these diagnostics are understood should `blh` be used in forcing.

## Success criteria

- Additional ERA5 boundary-layer fields are read with no NaNs/Inf.
- Diagnostics are produced without changing the control solution.
- Optional air-density forcing can be switched on and off cleanly.
- Any wind perturbation is vector-consistent and mask-limited.
- Boundary-layer tests remain separate from precipitation-phase tests.
