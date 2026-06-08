# Cross-project implementation roadmap

## Philosophy

This branch should test external forcing hypotheses one at a time around the same free-slip/lateral-drag CICE6 control. The code can share infrastructure in `ice_forcing.F90`, but the science should remain separable.

## Current state summary

| Project | State | Next action |
|---|---|---|
| ERA5 monthly reader | Working backbone | preserve and document monthly file contract |
| ERA5 precipitation phase | Working branch pathway | keep results out of repo; document implementation and sanity checks |
| ERA5 boundary-layer fields | forcing product available | read/diagnose before perturbing |
| Current-only tides | working current-only implementation | document implementation; run year-long hourly-output diagnostic |
| Dynamic mixed-layer depth | conceptual | build ORAS preprocessing and CICE `hmix` reader |
| Spectral waves | conceptual with crude prior scaffold | design around CAWCR/Noah Day boundary-edge spectra |

## Recommended implementation order

### Phase 1: preserve and document completed pathways

```text
1. Update README.md with project map and implementation principles.
2. Add ERA5 precipitation documentation.
3. Add current-only tide documentation.
4. Add source-notes document linking to upstream experiment notes.
```

### Phase 2: make the completed pathways reproducible

```text
5. Ensure monthly ERA5 reader is cleanly namelist-controlled.
6. Ensure rain/snow separation can be switched on/off.
7. Ensure coastal/form-factor snowfall scaling can be switched on/off.
8. Ensure tide-current forcing can be switched on/off independently of SSH forcing.
9. Add global MPI reductions for tide limiter counts.
```

### Phase 3: add diagnostics for seasonality mechanisms

```text
10. Add optional ERA5 boundary-layer diagnostics: pair, blh, gust factor, 100 m winds.
11. Add optional tide diagnostics: utide, vtide, tide_speed, daily RMS/max/exceedance counters.
12. Extend shuga to classify binary-days fast ice from hourly sea-ice output.
```

### Phase 4: ocean heat-capacity pathway

```text
13. Build ORAS preprocessing script in shuga.
14. Write monthly ORAS CICE-grid files with hmix and supporting diagnostics.
15. Add CICE dynamic-hmix reader.
16. Run fixed-hmix control and ORAS-hmix-14d experiment.
17. Only then test daily hmix, SST/SSS, and u/v additions.
```

### Phase 5: wave-boundary pathway

```text
18. Verify CAWCR/Zenodo spectral data layout.
19. Build shuga wave preprocessing for CICE-grid monthly files.
20. Add diagnostic wave reader in CICE.
21. Add wave-edge and bulk wave diagnostics.
22. Add active wave-ice coupling only after diagnostic runs are sane.
```

## Suggested run strategy

Use known-good initial conditions from the best lateral-drag experiment rather than restarting from zero global sea ice. Suggested short science window:

```text
primary branch year: 2004
optional extension:  2005
```

This keeps the forcing experiments short enough to be affordable while avoiding a fresh spin-up.

## Output strategy

Hourly global output is expensive but required for the first rigorous tide diagnostic. Use targeted streams where possible.

### For precipitation and boundary-layer tests

Daily output is usually sufficient:

```text
aice, hi, hs, uvel, vvel, fsnow, frain, Tair, uatm, vatm
```

### For current-only tide tests

Year-long hourly diagnostic output should include:

```text
aice, hi, uvel, vvel, uocn, vocn,
utide, vtide, tide_speed,
strocnx, strocny, strairx, strairy,
strintx, strinty, strength, divu, shear
```

### For dynamic mixed-layer tests

Daily output should include:

```text
hmix, sst, sss, Tf, qdp, frzmlt, uocn, vocn
```

### For wave tests

At least daily diagnostics:

```text
Hs, Tp, MWD, wave_edge_mask, wave_energy_edge,
wave_breakup_tendency, wave_fracture_mask
```

## `shuga` analysis additions

Required additions to `mawsons-chest/shuga`:

1. ERA5 monthly forcing builder already exists; keep extending it for boundary-layer variables.
2. ORAS forcing builder for `hmix`, SST/SSS, mixed-layer means, and OHC.
3. CAWCR wave spectral preprocessor for CICE-grid edge forcing files.
4. Hourly fast-ice classifier:

```text
hourly speed -> hourly raw FI -> daily binary state -> 11-day / 9-day binary-days FI_mask
```

5. Cross-mechanism attribution plots:

```text
tide_speed vs mobility
windgust / BLH vs mobility
wave energy at edge vs FI_mask loss
hmix / OHC vs growth rate
```

## Public repository boundary

This repository should include:

- scientific rationale;
- forcing file contracts;
- namelist controls;
- implementation stages;
- sanity-check expectations;
- links to source data and companion repos.

This repository should not include:

- unpublished paper-level result interpretation;
- raw history files, restart files, or large Zarr stores;
- claims that a forcing pathway improves or worsens FIA skill before manuscript figures are ready.
