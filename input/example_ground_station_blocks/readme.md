# Example ground station blocks

Typical ground station networks, ready to copy into the `<Scenario>` section of
`input/Config.xml` (replacing the `<GroundSegment>` block):

- `GroundSegmentESTRACK.xml` — the ESA ESTRACK core network: the three 35 m
  deep-space antennas (New Norcia, Cebreros, Malargüe) plus the Kourou, Kiruna,
  Redu and Santa Maria stations.
- `GroundSegmentNASA_DSN.xml` — the NASA Deep Space Network: the Goldstone,
  Madrid (Robledo) and Canberra (Tidbinbilla) complexes, roughly 120 deg apart
  in longitude for continuous deep-space coverage.
- `GroundSegmentNASA_NEN.xml` — a representative subset of the NASA Near Earth
  Network (Near Space Network direct-to-Earth): Wallops, Fairbanks, White Sands,
  Svalbard, McMurdo and South Point, strong on polar coverage for LEO missions.

Notes:
- Coordinates and heights are approximate site values (deg, m).
- `<ReceiverConstellation>` is set for a single constellation (`1`); extend the
  flag string with one `1`/`0` per ConstellationID when the scenario has more
  constellations (e.g. `1111`).
- The elevation masks are typical defaults (5 deg for near-Earth support, 10 deg
  for the deep-space antennas); adapt them to the mission.
