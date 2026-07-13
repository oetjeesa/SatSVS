# Project configurations

One directory per mission study, each with a self-contained `Config.xml` to
copy to `input/Config.xml` (or adapt). The ESA mission configs use
representative public orbit and instrument values (LTAN-defined SSO orbits
where applicable, epoch/simulation start 2026-02-01) with a Svalbard/Kiruna
ground segment; reduce the user grid or simulation window for quick looks.

## Copernicus Sentinel missions
- `Sentinel-1` — C-band SAR, 2 satellites, dawn-dusk SSO 693 km, IW swath
  (26–41 deg look angle), 12-day swath/revisit analysis
- `Sentinel-2` — MSI optical, 2 satellites, 10:30 LTDN SSO 786 km, 290 km
  nadir swath, 10-day revisit analysis
- `Sentinel-3` — OLCI/SLSTR, 2 satellites 140 deg apart, 10:00 LTDN SSO
  814.5 km, 1270 km swath
- `Sentinel-5P` — TROPOMI, 13:30 LTAN SSO 824 km, 2600 km swath (near-daily
  global coverage)
- `Sentinel-6` — radar altimeter in the 1336 km / 66 deg Jason reference
  orbit, 10-day repeat ground track

## Copernicus expansion (former High Priority Candidate) missions
- `CHIME` — hyperspectral imager, 2 satellites, 130 km swath
- `CIMR` — conically scanning microwave radiometer, ~1900 km swath
- `CO2M` — CO2 imaging spectrometer, 2 satellites, 250 km swath
- `CRISTAL` — polar ice altimeter, 88 deg non-SSO orbit, ground track
- `LSTM` — thermal infrared land surface temperature, 2 satellites, 384 km swath
- `ROSE-L` — L-band SAR on the Sentinel-1 orbit, ~260 km swath

## Earth Explorer missions
- `FLEX` — EE8, FLORIS spectrometer in tandem with Sentinel-3, 150 km swath
- `FORUM` — EE9, far-infrared sounder with MetOp-SG A, nadir ground track
- `Harmony` — EE10, two passive SAR companions +/-350 km along-track around a
  Sentinel-1 (3-satellite constellation)
- `Wivern` — EE11, conically scanning W-band Doppler radar, ~800 km swath

## Other studies
- `Ka_Interference` — Ka-band downlink interference between two co-planar LEOs
- `Polygon` — polygon-defined user segment example
- `S1-NG` — Sentinel-1 Next Generation swath/revisit study
- `S1-RCM` — Sentinel-1 with the Radarsat Constellation Mission (TLEs)
- `S1_ROSEL` / `S1_ROSEL_NISAR` — combined L-band SAR constellation studies
- `TANGO` — TANGO formation (TLE-based)
