# CICE_free-slip-forcing-sensitivity

`CICE_free-slip-forcing-sensitivity` is a development and experiment branch for testing whether the remaining Antarctic landfast-sea-ice response in the free-slip/lateral-drag CICE6 configuration is controlled by missing or over-simplified external forcing. The branch is not intended to replace the lateral-drag formulation. Instead, it treats the lateral-drag/free-slip model as the scientific control and layers targeted forcing experiments around it in `ice_forcing.F90`.

The forcing work is organised into five major projects:

1. **ERA5 precipitation phase and coastal snowfall forcing**: a thermodynamic fast-ice sensitivity pathway. This project modernises the ERA5 reader, separates rainfall and snowfall into CICE-aware fields, and provides a controlled way to perturb snowfall over Antarctic coastal/form-factor-favourable cells.
2. **ERA5 atmospheric boundary-layer forcing**: a seasonality and wind-coupling pathway. This project uses the expanded ERA5 monthly forcing product to expose surface pressure, boundary-layer height, wind gusts, and 100 m winds for diagnostics and later optional forcing perturbations.
3. **Current-only tidal forcing**: a sub-daily ocean-current pathway. This project introduces CATS2008_v2023 harmonic tidal currents into standalone CICE while keeping tidal SSH-gradient forcing disabled until the current-only mechanism has been isolated.
4. **Dynamic mixed-layer depth**: an ocean heat-capacity and restoring pathway. This project will replace a fixed or weakly varying `hmix` with ORAS-derived, time- and space-varying mixed-layer depth, with later options for mixed-layer-mean temperature, salinity, and velocity.
5. **Spectral wave forcing at the ice edge**: a marginal-ice-zone breakup and persistence pathway. This project will use CAWCR wave spectral data at the sea-ice boundary, following Noah Day's published wave-boundary concept more closely than the earlier crude `CICE_free-slip-waves` scaffold.

Two of these projects already have working code backbones: the monthly ERA5 precipitation/forcing reader and the current-only tide implementation. The remaining projects are staged as design documents so they can be implemented without entangling distinct physical hypotheses.

## Repository documentation map

| Project | Status | Main document |
|---|---|---|
| ERA5 precipitation phase and snowfall scaling | Implemented backbone; initial tests complete | [`docs/forcing/era5-precipitation.md`](docs/forcing/era5-precipitation.md) |
| ERA5 boundary-layer variables | Implemented in forcing product; CICE perturbations staged | [`docs/forcing/era5-boundary-layer.md`](docs/forcing/era5-boundary-layer.md) |
| Current-only tides | Implemented and tested through daily-history experiments | [`docs/forcing/current-only-tides.md`](docs/forcing/current-only-tides.md) |
| Dynamic ORAS mixed-layer depth | Conceptual design | [`docs/forcing/dynamic-mixed-layer-depth.md`](docs/forcing/dynamic-mixed-layer-depth.md) |
| Spectral wave forcing | Conceptual design | [`docs/forcing/spectral-wave-boundary.md`](docs/forcing/spectral-wave-boundary.md) |
| Cross-project implementation order | Living plan | [`docs/forcing/implementation-roadmap.md`](docs/forcing/implementation-roadmap.md) |

## Common implementation principles

1. Preserve the free-slip/lateral-drag control path and keep each forcing mechanism switchable from namelist settings.
2. Keep scientifically distinct pathways separable: precipitation, boundary-layer winds, tides, mixed-layer depth, and waves should each be independently attributable.
3. Prefer offline preprocessing in `mawsons-chest/shuga` and simple NetCDF ingestion in CICE. Do not rely on Zarr or Python-native data structures inside the Fortran model.
4. Keep Icepack changes out of the first implementation stage unless the physics explicitly requires an Icepack interface change.
5. Add diagnostic history only when needed to verify the forcing pathway or support downstream `shuga` attribution.
6. Avoid interpreting experiment skill in this repository. This repository should document the scientific rationale, implementation pathway, and numerical sanity checks; paper-level results and judgement should remain outside the public code branch until ready.

## Forcing-product contract

The current ERA5 forcing product is produced offline by [`mawsons-chest/shuga`](https://github.com/dpath2o/mawsons-chest/tree/main/shuga) and written as monthly NetCDF files on the CICE grid.

Expected file pattern:

```text
${atm_data_dir}/era5_for_cice6_YYYY_MM.nc
```

Example:

```text
/g/data/gv90/da1339/afim_input/ERA5/0p25/bilinear/monthly_cice6/era5_for_cice6_1994_10.nc
```

Core hourly fields:

| Variable | Meaning | Units | Initial CICE use |
|---|---|---:|---|
| `airtmp` | 2 m air temperature | K | atmospheric state |
| `spchmd` | 2 m specific humidity | kg kg-1 | atmospheric state |
| `pair` | surface pressure | Pa | optional air-density path |
| `glbrad` | downward shortwave radiation | W m-2 | atmospheric flux |
| `dlwsfc` | downward longwave radiation | W m-2 | atmospheric flux |
| `ttlpcp` | total precipitation rate | kg m-2 s-1 | fallback total precipitation |
| `snowfall` | snowfall rate | kg m-2 s-1 | preferred snow input |
| `rainfall` | rainfall rate | kg m-2 s-1 | preferred rain input |
| `wndewd` | eastward 10 m wind | m s-1 | wind forcing |
| `wndnwd` | northward 10 m wind | m s-1 | wind forcing |
| `blh` | boundary-layer height | m | diagnostic / future coupling |
| `windgust` | 10 m wind gust | m s-1 | diagnostic / future burstiness path |
| `wnd100ewd` | eastward 100 m wind | m s-1 | diagnostic / future shear path |
| `wnd100nwd` | northward 100 m wind | m s-1 | diagnostic / future shear path |

The first monthly-reader stage should reproduce the previous ERA5-forced behaviour using the core fields. Physical perturbations should be layered on only after the file discovery, record indexing, and field sanity checks are compile- and run-tested.

## Recommended branch philosophy

The branch should remain an implementation laboratory with small, reversible commits. Each commit should compile independently and should move only one forcing pathway forward. A useful sequence is:

```text
forcing: document integrated forcing-sensitivity roadmap
forcing: add/revise monthly ERA5 filename and record logic
forcing: read ERA5 phase-separated precipitation fields
forcing: add Antarctic coastal/form-factor snowfall scaling
forcing: read ERA5 boundary-layer diagnostic fields
forcing: document and preserve current-only tide implementation
forcing: add tide-current diagnostics needed for hourly persistence analysis
forcing: add ORAS-derived hmix reader and preprocessing contract
forcing: add wave-spectral boundary preprocessing contract
```
