# SatSVS analysis test suite

One folder per analysis type from `readme.md` (35 analyses incl. the orb_ and
sat_ series), plus the HPOP propagator benchmark and a multi-analysis run
(several `<Analysis>` blocks evaluated in one simulation). Each folder contains:

- `Config.xml` — the test scenario (self-contained; TLE files referenced by the
  config are copied into the folder as well)
- the plots and CSV data dumps the run produced — the CSVs are the golden
  reference data of the automated regression tests
- `main.log` — the simulator log of the run

## Automated regression tests (pytest)

```
py -m pytest tests                  # full suite (takes a while: MP4/HPOP tests)
py -m pytest tests -m "not hpop"    # skip the Orekit-based tests
py -m pytest tests -k sat_thermal   # a single test
py -m pytest tests --update-golden  # refresh the golden CSVs from this run
```

Each test runs its `Config.xml` through the simulator (into a scratch output
directory, via the `main.py <config> --output-dir <dir>` command line) and
compares every produced CSV data dump against the golden copy in the test
folder, number by number with tolerances (rtol 1e-6). Comparing the numeric
dumps instead of the PNGs keeps the tests robust against matplotlib/cartopy
version changes. A test folder without golden CSVs is skipped with a hint.

## Refreshing the reference outputs

```
py tests/run_test.py <test_name> [more names...]
py tests/run_test.py --all
```

The runner executes `python main.py <config> --output-dir <scratch>` from
`src/` headless (MPLBACKEND=Agg) and copies every produced output file (plots,
CSV dumps, movies, main.log) into the test folder — i.e. it refreshes both the
golden CSVs and the reference plots. `tests/make_configs.py` regenerates all
test configs from scratch.

## Scenarios

| Test | Scenario |
|---|---|
| cov_* (most) | 24-satellite GPS Walker constellation (Keplerian), global 10x10 deg user grid or static users (Delft, Singapore), Kourou station |
| cov_ground_track, cov_depth_of_coverage | TerraSAR-X TLE + SGP4, polar stations (Svalbard/Kiruna/Inuvik) |
| obs_swath_conical | MetOp-A TLE + SGP4, 650 km conical scanner, 2x2 deg global grid, revisit + Numpy export |
| obs_swath_push_broom, obs_sza_push_broom | Sentinel-1 TLE + SGP4, 250–650 km right-looking push broom, revisit + NetCDF export |
| obs_sza_subsat | Sentinel-1 TLE + SGP4, 2 days, SZA vs latitude + Numpy export |
| obs_aoi_revisit | 700 km SSO with 250–650 km right-looking push broom, Polygon user grid over Iberia/France at 1 deg, 3 days at 60 s — AOI revisit map + coverage fill-up curve |
| obs_target_imaging | 700 km SSO (RAAN 140 deg → noon/midnight LTAN), 45 deg agility cone, 10 deg minimum sun elevation, 5 targets incl. polar-night Longyearbyen, 2 days at 30 s |
| com_gr2sp_budget, com_doppler | TerraSAR-X TLE + SGP4, X-band (8.025 GHz) downlink to Svalbard/Kiruna, ITU-R attenuation, QPSK required C/N0 |
| com_gr2sp_budget_interference | 2 co-planar LEO SSO satellites (Keplerian, 1 deg mean-anomaly separation), Ka-band (26.25 GHz) to Svalbard |
| com_sp2sp_budget | GPS SV1→SV2 inter-satellite link at 23 GHz (space-to-space links enabled) |
| com_contact_plan | 2 co-planar SSO satellites 10 deg apart in mean anomaly, Svalbard + Inuvik, 1 day at 30 s — pass table with overlapping passes at every station (conflict detection) |
| com_pfd | TerraSAR-X TLE + SGP4, X-band isoflux pattern to Svalbard (mask 0 deg), 20 W over 300 MHz — PFD vs Article 21 style mask |
| nav_dilution_of_precision, nav_accuracy | GPS constellation, 10x10 deg global grid (elevation-dependent UERE profile for nav_accuracy) |
| sat_battery_depth_discharge | 700 km SSO satellite defined via LTAN (exercises the SSO/LTAN orbit path), payload latitude limit 60 deg |
| sat_eclipse_duration | 700 km SSO satellite, RAAN 140 deg → beta ~0 for max eclipses, 2 days at 30 s steps |
| sat_data_storage, sat_data_latency | 700 km SSO satellite, Svalbard + Inuvik downlink stations, 60 Mbps instrument / 1070 Mbps downlink |
| orb_kepler_elements | 250 km LEO, HPOP propagator with full force model (NRLMSISE00 drag), 3 days — evolution of all osculating Kepler elements |
| orb_air_density | Same 250 km LEO/HPOP scenario, 1 day — NRLMSISE00 density at the satellite altitude |
| orb_disturbance_forces | Same 250 km LEO/HPOP scenario, 1 day — per-epoch magnitude of every enabled perturbation |
| orb_pole_wobble | 700 km SSO, light HPOP force model, 60 days at 1800 s — IERS polar motion |
| orb_deltav_element | 250 km LEO/HPOP with drag, 2 days — altitude kept at 240 km mean +/- 1 km deadband |
| orb_beta_angle | ISS-like 420 km / 51.6 deg orbit, SGP4 (J2 nodal regression), 4 months at 1 h — beta beat cycle and eclipse fraction |
| orb_lifetime | 500 km SSO, Keplerian, 1 day at 60 s — drag decay projection to re-entry, 25-year rule and deorbit delta-v |
| orb_environment | 800 km SSO, Keplerian, 1 day at 60 s — eccentric-dipole trapped radiation (SAA + polar horns), dose-depth curve, atomic oxygen, Gruen micrometeoroids |
| sat_thermal | 700 km SSO satellite, RAAN 140 deg (beta ~0 for max eclipses), 1 day at 60 s — single-node thermal balance with per-orbit eclipse saw-tooth |
| sat_aocs | Same 700 km SSO scenario — worst-case disturbance torques (magnetic dominant at 700 km) and momentum buildup |
| multi_analysis | GPS constellation, static users (Delft, Singapore): cov_satellite_visible plus cov_satellite_sky_angles for SV1 and SV7, all three in a single run (the repeated type writes cov_satellite_sky_angles_2.png) |
| hpop_benchmark | TerraSAR-X TLE propagated with HPOP through cov_satellite_pvt; benchmarked by benchmark_hpop.py against a two-body analytic orbit and the SGP4 reference trajectory (reference_orbit_sgp4.txt, regenerated with Config_sgp4_reference.xml) |

## Results — all 27 analyses PASS

| Test | Outcome (visually verified) |
|---|---|
| cov_ground_track | SSO ground track over one day, correct pattern; also the 3D globe render (Plot3D with pyvista): textured Earth with cloud layer (EarthClouds), starry Milky Way background, track, orbit and satellite model in cov_ground_track_3d.png; MP4 movies (2D map filling up + 3D fly-along with the camera circling the satellite) |
| cov_satellite_pvt | ECI pos/vel sinusoids, 12 h GPS period, orbits.txt written |
| cov_satellite_visible | 6–10 GPS satellites in view per user |
| cov_satellite_visible_grid | Min 5–7 satellites in view worldwide; MP4 movies (2D instantaneous field per epoch + 3D fly-along of the GPS constellation with the statistic field on the globe) |
| cov_satellite_visible_id | Pass arcs of SV IDs 1–24 |
| cov_satellite_contour | MEO visibility contour at 10 deg mask |
| cov_satellite_sky_angles | Az/el passes for SV1 from Delft |
| cov_satellite_highest | Mean max-elevation map, symmetric bands |
| cov_depth_of_coverage | 0–2 stations in view along TerraSAR-X track |
| cov_pass_time | Mean pass 14000–27000 s, equatorial maximum |
| obs_swath_conical (+revisit) | Smooth semi-transparent conical swath ribbons (overlaps darker) + revisit diamonds + longitude-averaged max/mean revisit vs latitude profile (_revisit_lat.png); 3D swath ribbon render (Plot3D) in obs_swath_conical_3d.png; MP4 movies (2D ribbons filling up + 3D fly-along with the growing ribbon) |
| obs_swath_push_broom (+revisit) | Smooth 400 km semi-transparent ribbon strips + revisit latitude profile (max=mean after 1 day: single gap per point); NetCDF export works; 3D swath ribbon render (Plot3D) in obs_swath_push_broom_3d.png |
| obs_sza_push_broom | Mean SZA in swath, N–S gradient |
| obs_sza_subsat (+lat, +lat_year) | Dawn-dusk daylight SZA, plots labelled |
| obs_aoi_revisit | 84.5% of the Iberia/France AOI covered in 3 days (50% after 18 h), max revisit 36–48 h with the east–west push-broom asymmetry visible on the map |
| obs_target_imaging | 1–2 daylight opportunities per mid-latitude target (2–5 min each, ~23 h revisit), 0 for Longyearbyen in polar night — daylight constraint verified |
| com_gr2sp_budget | C/N0 vs required with the GRASP antenna patterns (isoflux Tx .cut, 3 m dish Rx .grd): realistic elevation-dependent margin crossings; pattern verification plot in _antenna.png |
| com_gr2sp_budget_interference | C/N0 with/without interferer; nominal gains at the GRASP pattern peaks (34.5/63.2 dBi), interferer discriminated through the actual pattern sidelobes |
| com_sp2sp_budget | Constant co-planar ISL geometry, FSL 202.5 dB |
| com_contact_plan | 52 passes/day over the two polar stations; every pass flagged as a conflict (the 10-deg separated satellites are always in view together), 280/220 min/day contact at Svalbard/Inuvik |
| com_pfd | Isoflux-flat PFD ~-160 dB(W/m2/4kHz) above 25 deg elevation, rolling off below the pattern edge; worst margin +19.1 dB vs the -150/-140 mask (compliant X-band EESS downlink) |
| com_doppler | ±180 kHz Doppler S-curves at 8.025 GHz |
| nav_dilution_of_precision | Max VDOP 1.8–5.2 worldwide |
| nav_accuracy | Mean horizontal accuracy 1.7–2.2 m (95%) |
| sat_battery_depth_discharge | Per-orbit eclipses, DoD swings 0–60% |
| sat_eclipse_duration | ~35 min eclipse per orbit (beta~0, 700 km) |
| sat_data_storage | SSR sawtooth 0–250 Gbit, downlink windows shaded |
| sat_data_latency | Mean 0.95 h, 95% 1.57 h, 100% < 2 h |
| orb_kepler_elements | HPOP + drag: 6-panel Kepler element evolution — J2 SMA oscillation over a 15.4 km secular drag decay in 3 days at ~250 km, +0.97 deg/day RAAN drift |
| orb_air_density | Density 4–8e-11 kg/m3 oscillating anti-phase with the 240–265 km altitude oscillation (perigee peaks) |
| orb_disturbance_forces | Textbook hierarchy at 250 km: J2 ~2e-2, drag ~2e-5, Moon > Sun > solid tides ~1e-6..5e-7, SRP square wave dropping to zero in eclipse, central gravity 9 m/s2 reference |
| orb_pole_wobble | xp 0.09→0.13, yp 0.36→0.41 arcsec over 60 days — arc of the annual/Chandler circle with sub-daily EOP loops |
| orb_deltav_element | Controlled altitude saw-tooths inside the deadband (8 boost maneuvers of 0.59 m/s = v·da/2a) while the uncontrolled orbit decays to 231 km; 4.7 m/s in 2 days = 887 m/s/year, consistent with the drag makeup rate |
| orb_beta_angle | Classic beta beat of a 51.6-deg LEO under J2 regression: ±60..72 deg swings over ~60-day cycles, eclipse fraction up to 38% at low beta and eclipse-free peaks above the 70-deg shadow half-angle |
| orb_lifetime | 500 km, Cd·A/m 0.014 m2/kg: 3.64-year decay with the accelerating knee below 400 km (25-year rule compliant); deorbit delta-v 109 m/s matches the half-Hohmann hand calculation |
| orb_environment | SAA localized over the South Atlantic (18.5% of orbit time above 100 protons/cm2/s), electron horns at high latitude; TID 879 rad/yr behind 4 mm Al, proton fluence 3.7e9 /cm2/yr, kapton erosion 0.1 um per 5 yr at 800 km, 179 impacts >1 ug on 20 m2 in 5 yr — SPENVIS orders of magnitude |
| sat_thermal | Per-orbit temperature saw-tooth −63..−21 degC between the eclipse cooling and sunlit heating, converging to a limit cycle from the first-epoch equilibrium; hot/cold case equilibria −11/−63 degC match the (Q/εσA)^0.25 hand calculation |
| sat_aocs | Textbook torque hierarchy at 700 km: magnetic ~4.6e-5 N m dominant (two peaks per orbit at the magnetic poles), SRP 3.7e-6 dropping to zero in eclipse, gravity gradient 2.4e-6 constant, aero 1.2e-6; momentum ramp 0.24 N m s/orbit |
| multi_analysis | Three analyses in one run: cov_satellite_visible plus sky angles for SV1 and SV7 with independent metric memory (different pass patterns per satellite); repeated type numbered as cov_satellite_sky_angles_2.png |
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
