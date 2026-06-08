# Dynamic ORAS-derived mixed-layer depth in CICE

## Purpose

This project tests whether Antarctic fast-ice growth and persistence are sensitive to realistic spatial and temporal variation in upper-ocean mixed-layer depth. The first implementation should be diagnostic-forcing-light: prescribe `hmix(i,j,t)` from ORAS-derived fields and keep the rest of the ocean pathway as close as possible to the current standalone CICE configuration.

This should precede any full prognostic mixed-layer rewrite.

## Scientific rationale

A fixed or weakly varying mixed-layer depth imposes an idealised ocean heat capacity under the ice. In coastal Antarctica, the depth of the mixed layer can vary strongly by region and season, especially across shelf seas, polynyas, stratified meltwater-influenced sectors, and deep winter convection regions. If the model's remaining fast-ice growth-rate or seasonal timing behaviour is controlled by ocean heat capacity or restoring strength, then spatially and temporally varying `hmix` is a more direct next experiment than further atmospheric tuning.

The immediate hypothesis is:

```text
A realistic ORAS-derived hmix changes the effective ocean heat capacity and restoring response beneath coastal Antarctic sea ice, thereby altering fast-ice growth, persistence, or retreat timing.
```

## Recommended staging

Do not begin with a full Petty-style prognostic mixed-layer implementation. Start with prescribed ORAS-derived `hmix`.

| Stage | Name | Description |
|---|---|---|
| M0 | fixed-control | current fixed `hmix_0` / existing restoring path |
| M1 | ORAS-hmix-14d | ORAS `hmix` updated on the same temporal order as `trestore` |
| M2 | ORAS-hmix-daily | daily ORAS `hmix` update |
| M3 | ORAS-hmix-clim | climatological seasonal-cycle `hmix` |
| M4 | ORAS-hmix-sstsss | dynamic `hmix` plus ORAS SST/SSS restoring |
| M5 | ORAS-hmix-uocn | add surface or mixed-layer-mean ORAS currents |
| M6 | prognostic-ML | later prognostic mixed-layer scheme, only if M1-M5 justify it |

## ORAS preprocessing contract

Preprocess ORAS/GREP data offline in `shuga` and write CICE-ready monthly files.

Suggested raw fields:

```text
mlotst_oras   mixed-layer thickness
thetao_oras   potential temperature with depth
so_oras       salinity with depth
uo_oras       eastward velocity with depth
vo_oras       northward velocity with depth
```

Recommended download subset for first tests:

```text
longitude:  global, -180 to 180
latitude:   -90 to -45, or -90 to -50 for smaller tests
depth:      0 to 500 m initially
time:       monthly files, daily records if available
```

Do not download the full-depth global ocean unless a later diagnostic proves it is necessary. For the first fast-ice-focused forcing experiment, 0-500 m should capture most upper-ocean heat-capacity and mixed-layer variability relevant to the CICE boundary condition.

## CICE-ready file contract

Suggested file pattern:

```text
${ocn_data_dir}/oras_for_cice6_YYYY_MM.nc
```

Suggested variables:

| Variable | Meaning | Units | Initial use |
|---|---|---:|---|
| `hmix` | ORAS-derived mixed-layer depth | m | primary forcing |
| `sst_oras` | ORAS near-surface temperature | degC or K, documented | optional restoring |
| `sss_oras` | ORAS near-surface salinity | psu | optional restoring / freezing temperature |
| `uocn_sfc` | ORAS surface eastward velocity | m s-1 | optional current forcing |
| `vocn_sfc` | ORAS surface northward velocity | m s-1 | optional current forcing |
| `uocn_mld_mean` | mixed-layer-mean eastward velocity | m s-1 | later current forcing |
| `vocn_mld_mean` | mixed-layer-mean northward velocity | m s-1 | later current forcing |
| `theta_mld_mean` | mixed-layer-mean temperature | degC or K, documented | diagnostic / later restoring |
| `salt_mld_mean` | mixed-layer-mean salinity | psu | diagnostic / later restoring |
| `ohc_mld` | upper-ocean heat content over MLD | J m-2 | diagnostic |
| `profile_mld_sigma` | independently recomputed density MLD | m | quality check |

The initial CICE implementation only needs `hmix`. The other fields should be written early because they are valuable for diagnostics and later stages.

## Fortran implementation outline

### Namelist controls

```fortran
logical :: use_dynamic_hmix
character(char_len) :: hmix_data_type     ! 'constant', 'oras', 'oras_clim'
character(char_len_long) :: hmix_data_dir
character(char_len_long) :: hmix_file_template
real(kind=dbl_kind) :: hmix_min
real(kind=dbl_kind) :: hmix_max
logical :: hmix_smooth_time
```

Defaults should preserve current behaviour:

```fortran
use_dynamic_hmix = .false.
hmix_data_type   = 'constant'
hmix_0           = 60.0
```

### Reader structure

Use the existing ocean-forcing infrastructure as much as possible:

```fortran
subroutine ORAS_files_monthly(yr, mon)
subroutine ORAS_hmix_data(dt)
```

Read two time records when interpolation is required:

```fortran
hmix_data(nx_block,ny_block,2,max_blocks)
```

then map to the existing CICE field:

```fortran
hmix(:,:,:) = hmix_forcing(:,:,:)
```

with bounds:

```fortran
hmix = max(hmix_min, min(hmix, hmix_max))
hmix = min(hmix, local_bathymetry_safe_limit)
```

### Update cadence

The first controlled experiment should update `hmix` at the same temporal order as the current ocean restoring timescale, e.g. 14 days. Daily update should be a second experiment, not the default.

## Interaction with Icepack

Avoid Icepack changes unless the existing interface assumes `hmix` is scalar or namelist-only. The preferred structure is:

```text
ORAS preprocessing -> CICE forcing reader -> hmix(i,j,t) -> existing Icepack/oceanmixed path
```

Icepack should receive a physically meaningful `hmix` field. It should not know about ORAS variable names, vertical coordinates, regridding, or Copernicus product conventions.

## Diagnostics

Required:

```text
hmix min/max/mean by timestep
number of clipped cells at hmix_min
number of clipped cells at hmix_max
number of cells limited by bathymetry
regional hmix mean by Antarctic sector
```

Optional history fields:

```text
hmix
sst_oras
sss_oras
uocn_sfc, vocn_sfc
uocn_mld_mean, vocn_mld_mean
ohcmld
```

## 1/12-degree product policy

Start with 0.25-degree ORAS/GREP. A 1/12-degree product should be considered only as an offline diagnostic sensitivity. The useful test is not whether 1/12-degree fields look better natively, but whether they still differ meaningfully after remapping to the 1/4-degree CICE grid.

Recommended offline comparison:

```text
0.25-degree ORAS -> CICE grid
1/12-degree product -> CICE grid
compare hmix, OHC, SST/SSS, surface currents in Antarctic coastal/LFI-favourable cells
```

If the remapped fields differ weakly, stay with 0.25 degree. If they differ strongly in fast-ice-relevant coastal sectors, then a targeted 1/12-degree forcing experiment may be justified.

## Success criteria

- Dynamic `hmix` reader compiles and runs with `use_dynamic_hmix = .true.`.
- Setting `use_dynamic_hmix = .false.` is bit-for-bit or practically equivalent to the control path.
- ORAS `hmix` is bounded, bathymetry-aware, and physically sane.
- The first experiment changes only `hmix`, not SST/SSS/u/v, so attribution is clean.
- Later SST/SSS/u/v stages are introduced one at a time.
