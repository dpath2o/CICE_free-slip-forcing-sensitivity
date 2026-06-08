# Spectral wave forcing at the sea-ice boundary

## Purpose

This project will introduce wave forcing at the Antarctic sea-ice boundary using CAWCR wave spectral data, with design priority given to Noah Day's published boundary-edge spectral workflow. The aim is to test whether ocean swell and wave-induced breakup influence Antarctic fast-ice persistence, retreat timing, and seasonal maximum FIA.

This is a conceptual design document. The existing `CICE_free-slip-waves` repository should be treated as a useful scaffold, not as the final implementation design.

## Scientific rationale

Fast ice can be mechanically stable under local wind and ocean-current forcing but still vulnerable to remotely generated swell. A daily or synoptic wave field acting at the sea-ice boundary may influence:

- marginal ice-zone breakup;
- floe-size distribution and effective ice strength;
- loss of coastal attachment during high-swell events;
- transition from persistent fast ice to mobile pack ice;
- regional retreat timing, especially in swell-exposed but geometrically sheltered sectors.

This pathway is distinct from tides. Tides are a sub-daily current/stress perturbation. Waves are a boundary-propagating spectral energy and breakup pathway.

## Existing crude scaffold in `CICE_free-slip-waves`

The previous wave branch already contains useful hooks:

```fortran
F_WAVE
wave_spec_dir
wave_spec_file
wave_file_template
wave_data(nx_block,ny_block,nfreq,2,max_blocks)
get_wave_spec
icepack_init_wave
```

This indicates a reasonable first attempt: read wave spectral information in `ice_forcing.F90` and pass it toward Icepack wave machinery.

However, the next version should not simply force wave spectra everywhere on the grid. It should use a sea-ice-boundary or marginal-ice-zone forcing concept, with spectral energy imposed at the ice edge and attenuated or diagnosed into the ice-covered region.

## Design priority

Follow the published CAWCR/Noah Day boundary-edge spectral approach as the primary scientific design. The previous crude CICE branch should provide:

- variable names and Fortran allocation patterns;
- proof that a wave reader can be inserted into `ice_forcing.F90`;
- a starting point for `wave_data` buffers;
- a warning about not over-simplifying spectral forcing into a single bulk wave scalar.

Before code implementation, verify the exact Zenodo record metadata and variable names for `https://zenodo.org/records/11081611`.

## Offline preprocessing contract

Preprocess CAWCR wave spectra in `shuga` into CICE-ready monthly files.

Suggested file pattern:

```text
${wave_spec_dir}/cawcr_wave_for_cice6_YYYY_MM.nc
```

Suggested dimensions:

```text
time
frequency
direction
nj
ni
```

Suggested variables:

| Variable | Meaning | Units | Use |
|---|---|---:|---|
| `efth` | directional wave energy spectrum | m2 s rad-1 or source documented | primary spectral input |
| `frequency` | wave frequency bins | Hz | spectral coordinate |
| `direction` | wave direction bins | deg or rad, convention documented | spectral coordinate |
| `hs` | significant wave height | m | diagnostic / bulk fallback |
| `tp` | peak period | s | diagnostic / bulk fallback |
| `tm` | mean period | s | diagnostic |
| `mwd` | mean wave direction | deg | diagnostic |
| `wave_edge_mask` | active ice-edge forcing cells | 0/1 | boundary forcing mask |
| `aice_edge_source` | ice-edge concentration used to define boundary | 0-1 | reproducibility |
| `dist_to_edge` | distance into ice from edge | m | attenuation/breakup diagnostic |

The exact field names should be adjusted to the CAWCR/Zenodo source product, but the CICE-facing file should be stable and documented.

## Ice-edge definition

The wave-boundary implementation needs a reproducible sea-ice edge. Candidate definitions:

```text
open ocean:       aice < 0.15
ice covered:      aice >= 0.15
active edge cell: ice-covered cell adjacent to open ocean
MIZ band:         ice-covered cells within N grid cells or D km of active edge
```

The first implementation should compute the edge offline from daily CICE or observational ice concentration, then write `wave_edge_mask` into the forcing file. A later implementation may compute the active edge online in CICE.

## Stage W0: diagnostic wave reader

### Goal

Read wave spectral files and verify the data path without changing CICE dynamics or thermodynamics.

Fortran requirements:

```fortran
logical :: use_wave_spectral_forcing
character(char_len_long) :: wave_spec_dir
character(char_len_long) :: wave_file_template
real(kind=dbl_kind), allocatable :: wave_data(:,:,:,:,:)
```

Diagnostic-only output:

```text
Hs, Tp, MWD
spectral energy integrated over frequency/direction
wave_edge_mask
```

## Stage W1: boundary-edge bulk diagnostic

Before activating a spectral breakup model, compute bulk diagnostics at the ice edge:

```text
Hs_edge
Tp_edge
wave_power_edge
wave_energy_flux_normal_to_edge
```

These diagnostics allow comparison of wave events with fast-ice mobility, retreat, and breakup without yet modifying the model state.

## Stage W2: spectral attenuation / breakup tendency

The first active implementation should be conservative. Candidate outputs:

```text
wave_breakup_tendency
wave_attenuation_length
wave_energy_into_ice
wave_fracture_mask
```

Possible model action, off by default:

```text
reduce effective floe size
increase mobility / reduce local ice strength proxy
alter a wave-fracture diagnostic used by Icepack if available
```

Do not directly delete fast ice or alter `FI_mask` inside CICE. Fast-ice classification should remain diagnostic/offline in `shuga` unless a later paper explicitly requires online coupling.

## Stage W3: active wave-ice coupling

Only after W0-W2 are verified should the branch modify model physics. Possible coupling choices:

1. pass spectra into existing Icepack wave-fracture machinery;
2. modify floe-size distribution if the active CICE/Icepack build supports it;
3. apply a wave-induced weakening or mobility tendency in marginal ice only;
4. add a breakup diagnostic that is later used by the classifier rather than by CICE dynamics.

## Required diagnostics

```text
wave_edge_mask area
Hs/Tp/MWD at active edge
spectral energy by sector
wave energy entering ice-covered cells
number/area of cells exceeding breakup threshold
co-occurrence with high ice speed or tide-current speed
co-occurrence with loss of binary-days FI_mask
```

## Relationship to tides and boundary-layer forcing

Waves should be implemented after the current-only tide diagnostics are stable. The first wave analysis should explicitly separate:

```text
tide-driven sub-daily mobility
wind-driven mobility / gust events
wave-driven edge breakup or retreat
```

This is important because all three mechanisms can produce increased ice speed or reduced persistence, but they represent different physics and different implementation pathways.

## Success criteria

- The CAWCR/Zenodo spectral data can be preprocessed into monthly CICE-grid files.
- CICE can read wave spectra or derived bulk fields without changing the control solution.
- The wave-edge mask is reproducible and sector-aware.
- Active coupling remains off until diagnostics show the forcing path is physically sane.
- The implementation follows the published boundary-edge spectral concept rather than the earlier crude whole-grid forcing scaffold.
