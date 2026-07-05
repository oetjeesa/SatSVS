# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Satellite Service Volume Simulator (SatSVS) — a Python tool for satellite mission
analysis. It propagates satellite orbits, rotates ground stations and users in
ECI/ECF, computes visibility links between them, and runs one analysis per run
(coverage, Earth observation, communication link budget, or navigation DOP/accuracy).

See `readme.md` for the full catalogue of analyses and the complete `Config.xml`
schema — it is the authoritative user-facing reference for every analysis type and
its XML parameters.

## Running

The tool is run as a script, not a package. There is no build, lint, or test suite.
The working directory matters: nearly every path is hardcoded relative to `src/`
(e.g. `../input/Config.xml`, `../output/`), and `sys.path.append('../src')` appears
at the top of each module. **Run from inside `src/`:**

```
cd src
python main.py
```

There are no requirements/environment files. Key third-party dependencies: `numpy`,
`astropy`, `sgp4` (v2.x — the accelerated `sgp4.api.Satrec` API), `numba`,
`matplotlib` + `cartopy`, `geopandas`/`shapely`, and `itur` (ITU-R propagation
models, used only by `com_*` analyses). The HPOP propagator additionally needs
`orekit_jpype` + `jdk4py` and the Orekit data archive at `input/orekit-data.zip`.

## Configuration drives everything

A single `input/Config.xml` (capital C) defines the whole simulation: space segment
(constellations/satellites), ground segment (stations), user segment, simulation
window, and **exactly one** `<Analysis>` block. The `projects/` and
`input/example_conf_blocks/` directories hold reusable config fragments and full
example configs to copy from. Logs are written to `output/main.log`.

## Architecture

The flow is a classic load → time-loop → plot pipeline, all orchestrated in
`main.py`. The `AppConfig` object (`config.py`), referred to everywhere as `sm`
("state machine"), is the single mutable state container passed to every function and
analysis hook — it holds the satellite/station/user lists, the link matrices, and the
current loop time.

**Startup (`main.load_configuration`)** parses `Config.xml` four times — once each for
`load_satellites`, `load_stations`, `load_users`, `load_simulation` — then
`setup_links`. `load_simulation` is also where the analysis type string is mapped to a
concrete analysis class instance (a long `if` chain); **adding a new analysis requires
adding a branch here** plus instantiating the class.

**Segments (`segments.py`)** define the domain objects: `Satellite`, `Station`,
`User`, `Constellation`, and three link classes `Ground2SpaceLink`, `User2SpaceLink`,
`Space2SpaceLink`. Position is held redundantly in ECI and ECF; each object has
`det_posvel_eci*` / `det_posvel_ecf` methods called every epoch. Links have
`compute_link` (geometry) and `check_masking` (elevation mask + Earth-blockage test).

**Link matrices** are 2D lists (`sm.gr2sp[station][sat]`, `sm.usr2sp[user][sat]`,
`sm.sp2sp[sat][sat]`). `setup_links` precomputes a `link_in_use` flag per pair from the
`ReceiverConstellation` bitstring (a per-constellation '1'/'0' mask), so the time loop
only computes geometry for enabled pairs.

**Time loop (`main.run_time_loop`)** iterates `num_epoch` steps. Each epoch: convert
MJD→GMST, propagate every satellite (Keplerian, SGP4 or HPOP, selected by
`<OrbitPropagator>`; HPOP is the Orekit-based numerical propagator in
`propagation_hpop.py`, configured by the `<HPOP>` block and sampling per-satellite
dense ephemerides generated at load time), update stations/users, recompute in-view
lists (`idx_sat_in_view`, etc.), then call `sm.analysis.in_loop(sm)`. Setting
`OrbitsFromPreviousRun` reads cached ECI orbits from `output/orbits_internal.txt`
instead of re-propagating.

**Analyses** subclass `AnalysisBase` (`analysis.py`) and implement four hooks:
`read_config(node)`, `before_loop(sm)`, `in_loop(sm)`, `after_loop(sm)` — the loop in
`main.py` calls them at the matching phases, and `after_loop` is where plots
(`../output/<type>.png`) are produced. Concrete analyses live grouped by domain in
`analysis_cov.py`, `analysis_obs.py`, `analysis_com.py`, `analysis_nav.py`;
`analysis.py` also holds the cartopy map helpers (`make_map_cyl`, `make_map_polar`)
and shared plotting mixins (e.g. `AnalysisObs` swath/revisit plotting on a map).
`misc_fn.py` holds geometry/time helpers (MJD/GMST, coordinate
conversions); `constants.py` holds physical constants (`R_EARTH`, `GM_EARTH`, etc.).

### Adding a new analysis (per readme.md)

1. Add an analysis class at the end of the relevant `analysis_*.py`, using an existing
   class or `AnalysisBase` as a template.
2. Add the type-string → class branch in `config.py` `load_simulation`.

## Conventions

- Angles in config XML are **degrees**; everything internal is **radians** (conversion
  happens at load time in `config.py`). Distances are SI (metres).
- `ConstellationID` is 1-based; it indexes into the `ReceiverConstellation` bitstring
  as `id - 1`.
- Output PNGs and data dumps go to `output/`; the reference images checked into `docs/`
  are the expected results for each analysis.
