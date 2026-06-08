# Source notes

This documentation bundle integrates material from the following project sources.

## User-provided README draft

The uploaded draft README for `CICE_free-slip-forcing-sensitivity` provided the original ERA5 monthly forcing contract, precipitation-stage plan, boundary-layer variable list, and the high-level framing of atmospheric, oceanic, tide, and wave forcing work.

Key items carried forward:

- monthly ERA5 files written by `mawsons-chest/shuga`;
- variables `airtmp`, `spchmd`, `pair`, `glbrad`, `dlwsfc`, `ttlpcp`, `snowfall`, `rainfall`, `wndewd`, `wndnwd`, `blh`, `windgust`, `wnd100ewd`, `wnd100nwd`;
- monthly file pattern `era5_for_cice6_YYYY_MM.nc`;
- ERA5 precipitation phase separation;
- Antarctic coastal/form-factor snowfall scaling;
- future mixed-layer-depth and wave-forcing pathways.

## CICE_free-slip-tides documentation

Integrated from:

```text
https://github.com/dpath2o/CICE_free-slip-tides/tree/docs/current-only-tide-experiments/docs/experiments/current-only-tides
```

Important implementation details retained:

- current-only CATS2008_v2023 harmonic tide pathway;
- `tide_use_currents = .true.`;
- `tide_use_ssh = .false.`;
- bathymetry-aware limiter;
- half-cosine ramp;
- diagnostic emphasis on ocean-stress perturbation and near-threshold fast-ice mobility;
- need to test persistence explicitly under the binary-days classifier.

The current README list includes experiment notes for 10-day, one-month, one-year, and later higher-frequency current-only tide experiments. Large history files, restarts, and Zarr stores are intentionally excluded from those notes.

## Project chat context

Integrated project context retained:

- tide implementation should remain current-only before adding SSH forcing;
- bathymetry-aware tide forcing is preferred over a hard cap as the primary limiter;
- long tide/no-tide comparisons require output fields sufficient for fast-ice classification and dynamical attribution;
- `shuga` should be extended to classify binary-days fast ice from higher-frequency sea-ice output.

## CICE_free-slip-waves scaffold

Integrated from:

```text
https://github.com/dpath2o/CICE_free-slip-waves
```

Important scaffold items retained:

- `F_WAVE`;
- `wave_spec_dir`;
- `wave_spec_file`;
- `wave_file_template`;
- `wave_data(nx_block,ny_block,nfreq,2,max_blocks)`;
- public `get_wave_spec` hook;
- `icepack_init_wave` import path.

These are treated as useful code-path hints, not as the final scientific design.

## CAWCR / Noah Day wave data

The requested design should be based primarily on the CAWCR spectral wave data and Noah Day's published boundary-edge implementation concept:

```text
https://zenodo.org/records/11081611
```

Before coding, verify the exact Zenodo metadata, variable names, units, direction convention, time coordinate, and spectral dimensions. The documentation here intentionally defines a stable CICE-facing contract rather than assuming the raw Zenodo variable names.

## Companion repository links

```text
mawsons-chest / shuga:
https://github.com/dpath2o/mawsons-chest/tree/main/shuga

CICE_free-slip-tides current-only experiment notes:
https://github.com/dpath2o/CICE_free-slip-tides/tree/docs/current-only-tide-experiments/docs/experiments/current-only-tides

CICE_free-slip-waves:
https://github.com/dpath2o/CICE_free-slip-waves

CAWCR/Noah Day wave data record:
https://zenodo.org/records/11081611
```
