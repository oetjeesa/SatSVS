# SatSVS analysis test suite

One folder per analysis type from `readme.md` (25 analyses incl. orb_semi_major_axis),
plus the HPOP propagator benchmark. Each folder contains:

- `Config.xml` — the test scenario (self-contained; TLE files referenced by the
  config are copied into the folder as well)
- the plots/data files the run produced (copied from `output/`)
- `main.log` — the simulator log of the run

## How to run

```
py tests/run_test.py <test_name> [more names...]
```

The runner copies the test `Config.xml` to `input/Config.xml`, runs
`python main.py` from `src/` headless (MPLBACKEND=Agg), and copies every output
file the run touched back into the test folder. `tests/make_configs.py`
regenerates all test configs from scratch.

## Scenarios

| Test | Scenario |
|---|---|
| cov_* (most) | 24-satellite GPS Walker constellation (Keplerian), global 10x10 deg user grid or static users (Delft, Singapore), Kourou station |
| cov_ground_track, cov_depth_of_coverage | TerraSAR-X TLE + SGP4, polar stations (Svalbard/Kiruna/Inuvik) |
| obs_swath_conical | MetOp-A TLE + SGP4, 650 km conical scanner, 2x2 deg global grid, revisit + Numpy export |
| obs_swath_push_broom, obs_sza_push_broom | Sentinel-1 TLE + SGP4, 250–650 km right-looking push broom, revisit + NetCDF export |
| obs_sza_subsat | Sentinel-1 TLE + SGP4, 2 days, SZA vs latitude + Numpy export |
| com_gr2sp_budget, com_doppler | TerraSAR-X TLE + SGP4, X-band (8.025 GHz) downlink to Svalbard/Kiruna, ITU-R attenuation, QPSK required C/N0 |
| com_gr2sp_budget_interference | 2 co-planar LEO SSO satellites (Keplerian, 1 deg mean-anomaly separation), Ka-band (26.25 GHz) to Svalbard |
| com_sp2sp_budget | GPS SV1→SV2 inter-satellite link at 23 GHz (space-to-space links enabled) |
| nav_dilution_of_precision, nav_accuracy | GPS constellation, 10x10 deg global grid (elevation-dependent UERE profile for nav_accuracy) |
| pow_battery_depth_discharge | 700 km SSO satellite defined via LTAN (exercises the SSO/LTAN orbit path), payload latitude limit 60 deg |
| pow_eclipse_duration | 700 km SSO satellite, RAAN 140 deg → beta ~0 for max eclipses, 2 days at 30 s steps |
| dat_storage, dat_latency | 700 km SSO satellite, Svalbard + Inuvik downlink stations, 60 Mbps instrument / 1070 Mbps downlink |
| orb_semi_major_axis | 250 km LEO, HPOP propagator with full force model (NRLMSISE00 drag), 3 days — semi-major axis decay |
| hpop_benchmark | TerraSAR-X TLE propagated with HPOP through cov_satellite_pvt; benchmarked by benchmark_hpop.py against a two-body analytic orbit and the SGP4 reference trajectory (reference_orbit_sgp4.txt, regenerated with Config_sgp4_reference.xml) |

## Results — all 24 analyses PASS

| Test | Outcome (visually verified) |
|---|---|
| cov_ground_track | SSO ground track over one day, correct pattern; also the 3D globe render (Plot3D with pyvista): textured Earth, track, ECF orbit and satellite model in cov_ground_track_3d.png |
| cov_satellite_pvt | ECI pos/vel sinusoids, 12 h GPS period, orbits.txt written |
| cov_satellite_visible | 6–10 GPS satellites in view per user |
| cov_satellite_visible_grid | Min 5–7 satellites in view worldwide |
| cov_satellite_visible_id | Pass arcs of SV IDs 1–24 |
| cov_satellite_contour | MEO visibility contour at 10 deg mask |
| cov_satellite_sky_angles | Az/el passes for SV1 from Delft |
| cov_satellite_highest | Mean max-elevation map, symmetric bands |
| cov_depth_of_coverage | 0–2 stations in view along TerraSAR-X track |
| cov_pass_time | Mean pass 14000–27000 s, equatorial maximum |
| obs_swath_conical (+revisit) | Classic conical swath coverage/revisit diamonds; 3D swath ribbon render (Plot3D) in obs_swath_conical_3d.png |
| obs_swath_push_broom (+revisit) | 400 km strips; NetCDF export works; 3D swath ribbon render (Plot3D) in obs_swath_push_broom_3d.png |
| obs_sza_push_broom | Mean SZA in swath, N–S gradient |
| obs_sza_subsat (+lat, +lat_year) | Dawn-dusk daylight SZA, plots labelled |
| com_gr2sp_budget | C/N0 95–110 dBHz vs required, ITU-R attenuations |
| com_gr2sp_budget_interference | C/N0 with/without interferer, antenna patterns |
| com_sp2sp_budget | Constant co-planar ISL geometry, FSL 202.5 dB |
| com_doppler | ±180 kHz Doppler S-curves at 8.025 GHz |
| nav_dilution_of_precision | Max VDOP 1.8–5.2 worldwide |
| nav_accuracy | Mean horizontal accuracy 1.7–2.2 m (95%) |
| pow_battery_depth_discharge | Per-orbit eclipses, DoD swings 0–60% |
| pow_eclipse_duration | ~35 min eclipse per orbit (beta~0, 700 km) |
| dat_storage | SSR sawtooth 0–250 Gbit, downlink windows shaded |
| dat_latency | Mean 0.95 h, 95% 1.57 h, 100% < 2 h |
| orb_semi_major_axis | HPOP + drag: 18.1 km secular SMA decay over 3 days at ~250 km |
| hpop_benchmark | Two-body HPOP vs analytic Kepler: max 0.024 m/day (PASS); full-force HPOP vs SGP4 reference: RMS 6.7 km, max 10.8 km/day (PASS) |

## 3D plots

Every world-map analysis test enables `<Plot3D>True</Plot3D>` and stores the 3D
globe render as `<type>_3d.png` in its folder. The GPS grid tests also
demonstrate the flags: `cov_pass_time` uses `<ShowSatellite>False</ShowSatellite>`
(orbit shells only), `cov_satellite_visible_grid` uses
`<ShowOrbit>False</ShowOrbit>` (satellite models only).

## HPOP benchmark

```
py tests/run_test.py hpop_benchmark          # produces orbits.txt (HPOP trajectory)
py tests/hpop_benchmark/benchmark_hpop.py    # writes benchmark_report.txt + diff plot
```

## Bugs found and fixed during testing

1. **`src/misc_fn.py` — `calc_az_el_dist_batch`**: the batched space-to-space
   link kernel included the satellite-to-itself pair; the zero line-of-sight
   vector made `calc_az_el` divide 0/0 and numba raised `ZeroDivisionError`.
   Any run with `IncludeSpace2SpaceLinks=True` crashed (regression from the
   time-loop batching speedup). Zero-distance targets now get az/el = 0 and are
   skipped (such self-links are never in use).
2. **`src/misc_fn.py` — `kep2xyz` (and `newton_raphson`)**: the Newton-Raphson
   convergence test `|E[k+1]/E[k] - 1|` divides by zero when the mean anomaly
   is exactly 0 (e.g. requested epoch equals the element epoch with
   `MeanAnomaly=0`). Replaced with an absolute-increment test
   `|E[k+1] - E[k]| > 1e-13`, which is equally accurate and safe at M=0.
3. **`src/analysis_pow.py` — `AnalysisPowEclipseDuration.after_loop`**: when a
   simulation contains no eclipses (high-beta orbit), the code logged a warning
   but still indexed the empty array → `IndexError`. Now returns after the
   warning.
4. **`src/analysis_obs.py` — `plot_sza_latitude` / `plot_sza_latitude_year`**:
   `plt.savefig()` was called before the axis labels and grid were set, so the
   saved PNGs had unlabelled axes. Labels/grid now set before saving.
5. **`src/analysis_obs.py` — `AnalysisObsSzaSubSat.__init__`**: `self.save_output`
   was never initialised; a config without `<SaveOutput>` crashed with
   `AttributeError` in `plot_sza_latitude`. Now defaults to `None`.
6. **`src/analysis_com.py` — `AnalysisComSp2SpBudget.__init__`**: `self.sat_id1`
   was initialised twice (typo); a config without `<SatelliteID2>` crashed with
   `AttributeError` instead of failing cleanly. Second line now initialises
   `self.sat_id2`.

Scenario note (not a code bug): an SSO orbit plane near the terminator
(RAAN ~50 deg in February) is genuinely eclipse-free — the eclipse test uses
RAAN 140 deg so eclipses actually occur; the empty-eclipse crash it originally
triggered is fix #3.
