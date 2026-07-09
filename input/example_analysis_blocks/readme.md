# Example analysis blocks

One file per analysis type, ready to copy into the `<SimulationManager>` section of
`input/Config.xml`. One or more `<Analysis>` blocks can be active in the same run:
all analyses are computed over the same propagated orbits, each writing its own
output files (named after the Type; repeated Types are numbered `_2`, `_3`, ...).

Values shown are sensible starting points; lines marked Optional can be removed.
See the main `readme.md` for the meaning and units of every parameter, and note
that some analyses need specific segments:

- `cov_pass_time`, `cov_satellite_highest`, `cov_satellite_visible_grid`,
  `nav_dilution_of_precision`, `nav_accuracy`: user segment of Type `Grid`
- `obs_*` swath analyses: instrument definition (`ObsSwath*`/`ObsIncidenceAngle*`)
  in the `<Constellation>` block; a user grid for the revisit statistics
- `com_sp2sp_budget`: `<IncludeSpace2SpaceLinks>True</IncludeSpace2SpaceLinks>`
- `nav_accuracy`: `<UERE>` list in the `<Constellation>` block
- World-map analyses additionally accept the `<Plot3D>` parameters (see the
  3D plots section of the main readme)
