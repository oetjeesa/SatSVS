# Project configurations

One directory per mission study, each with a self-contained `Config.xml` to
copy to `input/Config.xml` (or adapt). The mission configs use representative
public orbit and instrument values (LTAN-defined SSO orbits where applicable,
epoch/simulation start 2026-02-01) with a Svalbard/Kiruna ground segment;
reduce the user grid or simulation window for quick looks.

## Copernicus Sentinel missions
- `Sentinel-1` — C-band SAR, 2 satellites, dawn-dusk SSO 693 km, IW swath
  (26–41 deg look angle), 12-day swath/revisit analysis
- `Sentinel-2` — MSI optical, 2 satellites, 10:30 LTDN SSO 786 km, 290 km
  nadir swath, 10-day revisit analysis
- `Sentinel-3` — OLCI/SLSTR, 2 satellites 140 deg apart, 10:00 LTDN SSO
  814.5 km, 1270 km swath

Only these reference missions are kept in the repository; further mission
study folders (other Copernicus/Earth Explorer missions, interference and
formation studies) live outside version control.
