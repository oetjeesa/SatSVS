# Satellite Service Volume Simulator
## Open Source satellite, ground station and user tool
### M. Tossaint - 2026 - v3

<img src="/docs/schema.png" alt="schema"/>

## Installation & first run
Download from github and install the dependencies from the included
__pyproject.toml__ (run from the repository root):

```
pip install -e .                # core dependencies
pip install -e .[hpop]          # + HPOP (Orekit) numerical propagator
pip install -e .[plot3d]        # + 3D globe plots (pyvista)
pip install -e .[movie]         # + MP4 movies (imageio/imageio-ffmpeg)
pip install -e .[test]          # + pytest for the regression test suite
```

The core dependencies are: numpy, pandas, scipy, numba, astropy, sgp4, matplotlib,
cartopy, xarray, netcdf4, geopandas, shapely and itur. Some features need extra
data files:
- For the HPOP orbit propagator only: orekit_jpype and jdk4py (bundled JVM), plus the
  Orekit physical data archive saved as __input/orekit-data.zip__
  (download from https://gitlab.orekit.org/orekit/orekit-data)
- For the 3D plots only: pyvista, plus an equirectangular Earth texture saved as
  __input/earth_texture.jpg__ (e.g. the public domain NASA Blue Marble image
  https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57752/land_shallow_topo_2048.jpg).
  The starry background uses __input/starmap.jpg__ (Milky Way panorama, ESO/S. Brunier,
  CC BY 4.0, https://www.eso.org/public/images/eso0932a/, shipped downscaled and
  dimmed for a subtle background — plain black background if absent) and the
  optional cloud layer (EarthClouds) uses __input/earth_clouds.jpg__
  (e.g. the NASA cloud composite
  https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57747/cloud_combined_2048.jpg)
- For the MP4 movies only (MP4 flag): imageio and imageio-ffmpeg

To run, edit the configuration file and run from inside `src/`:

```
cd src
python main.py                                        # uses ../input/Config.xml -> ../output/
python main.py path/to/Config.xml                     # any scenario file
python main.py path/to/Config.xml --output-dir out1   # per-run output directory
```

The optional command line arguments make batch/trade-study runs easy: every run
can read its own scenario file and write plots, data dumps and main.log into its
own output directory. The configuration is validated at startup — misspelled
tags produce a "did you mean" warning and missing/invalid required parameters
stop the run with a clear message instead of a mid-run crash.

Next to every plot the analysis writes the plotted numbers as a CSV file
(__output/<analysis_type>.csv__, with a header line naming the columns), so
results can be post-processed without rerunning the simulation.

## Introduction
Framework takes care of geometry computations, satellite propagation, ground station and user rotation in ECI/ECF.
It will also automatically compute links between stations and satellites, users and satellites, and between satellites.

### Main structure of the tool

<img src="/docs/satsvs_architecture.png" alt="schema"/>

### Configuration of the tool

Configuration of the tool can be done in the config.xml file where satellites, ground stations, users, simulation
parameters and analysis are defined. Reusable configuration fragments to copy from are available in
__input/example_conf_blocks__ (constellations and user segments), __input/example_analysis_blocks__
(one example per analysis type below) and __input/example_ground_station_blocks__ (typical ground station
networks: ESA ESTRACK, NASA DSN and NASA NEN).
Analysis can be added as wished, the baseline of analysis available are below 
(and explained further below):

### Orbit
- __orb_kepler_elements__: Evolution of the osculating Kepler elements over time
- __orb_air_density__: Atmospheric density at the satellite altitude over time (HPOP)
- __orb_disturbance_forces__: Magnitude of the perturbation accelerations over time (HPOP)
- __orb_pole_wobble__: Polar motion (wobble of the Earth rotation axis) over time (HPOP)
- __orb_deltav_element__: Station-keeping delta-v to hold one orbit element in a deadband
  (semi-major axis, eccentricity, inclination, RAAN, argument of perigee, mean
  anomaly), e.g. the orbital decay under atmospheric drag with the HPOP propagator
- __orb_beta_angle__: Solar beta angle and analytic eclipse fraction over time
- __orb_lifetime__: Orbital lifetime under drag, 25-year rule compliance and deorbit delta-v
- __orb_environment__: Space environment along the orbit (trapped radiation/SAA, dose vs. shielding, atomic oxygen, micrometeoroids)

### Coverage
- __cov_ground_track__: Ground Track of Satellite on map
- __cov_satellite_pvt__: Output satellite position and velocity in ECI (for later use)
- __cov_satellite_visible__: Number of satellites in view over time for user (also spacecraft user)
- __cov_satellite_visible_grid__: Number of satellites in view statistics for user grid
- __cov_satellite_visible_id__: Satellite IDs in view over time for selected user (also spacecraft user)
- __cov_satellite_contour__: Satellite ground visibility contour at StopDate
- __cov_satellite_sky_angles__: Satellite Azimuth and Elevation for selected user (also spacecraft user)
- __cov_satellite_highest__: Satellite elevation statistics for highest satellite in view for user grid
- __cov_depth_of_coverage__: Depth of station coverage (DOC) based on satellite ground track
- __cov_pass_time__: Satellite passes time statistics for user grid

### Earth observation
- __obs_swath_conical__: Swath coverage for satellite(s) with conical scanners
- __obs_swath_pushbroom__: Swath coverage for satellite(s)
- __obs_sza_push_broom__: Solar Zenith Angle (SZA) within the push broom swath of satellite(s)
- __obs_sza_subsat__: Solar Zenith Angle (SZA) for all satellite(s) sub satellite point
- __obs_aoi_revisit__: Revisit and coverage build-up statistics over a polygon area of interest
- __obs_target_imaging__: Imaging opportunities for point targets within a pointing agility cone (optional daylight constraint)

### Communication
- __com_gr2sp_budget__: For station-satellite received power, losses and C/N0
- __com_sp2sp_budget__: For satellite-satellite received power, losses and C/N0
- __com_doppler__: For satellite-station elevation and doppler
- __com_contact_plan__: Ground station pass table (AOS/LOS/duration/volume) with conflict detection
- __com_pfd__: Power flux density at the ground vs. ITU-R Article 21 style limit mask

### Navigation
- __nav_dillution_of_precision__: DOP values for user(s) (also spacecraft user)
- __nav_accuracy__: Navigation accuracy (UERE*DOP) values for user(s) (also spacecraft user)

### Satellite
- __sat_thermal__: Single-node thermal balance over orbit (solar/albedo/Earth IR/eclipses)
- __sat_aocs__: AOCS disturbance torques and momentum buildup over orbit
- __sat_battery_depth_discharge__: Battery state-of-charge / depth-of-discharge and power generation vs. draw over orbit
- __sat_eclipse_duration__: Eclipse duration per orbit over the simulation
- __sat_data_storage__: Solid State Recorder (SSR) fill level over orbit (recording vs. downlink)
- __sat_data_latency__: Data latency statistics from acquisition to ground reception
- __sat_drag_coefficient__: Drag coefficient from the satellite geometry (Sentman free-molecular panel method on the STL model)

## Configuration file
The configuration file is found in the input directory under the name __config.xml__. 
This file contains several main parts:
- The space segment in several constellations
- The ground segment with ground stations servicing the satellites
- The user segment receiving the satellite information at single locations or user grids
- The simulation parameters like start,stop,step time and analysis to be performed.

All units in the configuration file are in SI units, angles are in degrees.

### Space segment
The following xml is used to setup the space segment:
```
<?xml version="1.0" encoding="utf-8"?>
<!--Simulation Scenario-->
<Scenario>
    <SpaceSegment>
    <!--Defines the constellations and satellites that should be simulated-->
        <Constellation>
            <NumOfPlanes>3</NumOfPlanes>
            <NumOfSatellites>30</NumOfSatellites>
            <ConstellationID>1</ConstellationID>
            <ConstellationName>Galileo</ConstellationName>
            <Satellite>
                <SatelliteID>1</SatelliteID>
                <Plane>1</Plane>
                <EpochMJD>54465.5</EpochMJD>
                <SemiMajorAxis>29600318</SemiMajorAxis>
                <Eccentricity>0.00</Eccentricity>
                <Inclination>56.00</Inclination>
                <RAAN>0</RAAN>
                <ArgOfPerigee>0.00</ArgOfPerigee>
                <MeanAnomaly>0</MeanAnomaly>
            </Satellite>
        </Constellation>
    </SpaceSegment>

```
The satellite can also be defined as an Sun Synchronous Satellite (SSO) with:
- LTAN Local Time Ascending Node in hours
- Altitude in m
The inclination is computed based on the SSO assumption. Put the Epoch at midday, ArgOfPerigee at 0 and MeanAnomaly at 0,
so that the satellite passes the equator at ltan requested.
```
<Satellite>
    <SatelliteID>1</SatelliteID>
    <Plane>1</Plane>
    <EpochMJD>58945.15</EpochMJD>
    <Altitude>694000</Altitude>
    <Eccentricity>0.0001402</Eccentricity>
    <LTAN>22.25</LTAN>
    <ArgOfPerigee>0.0</ArgOfPerigee>
    <MeanAnomaly>0.0</MeanAnomaly>
</Satellite>
```
Orbit parameters are either Keplerian or can be defined as a list of satellites in a TLE file (in input directory):
```
<Constellation>
    <NumOfPlanes>3</NumOfPlanes>
    <NumOfSatellites>30</NumOfSatellites>
    <ConstellationID>1</ConstellationID>
    <ConstellationName>OrbComm</ConstellationName>
    <TLEFileName>input/tle.txt</TLEFileName>
</Constellation>
```

With the optional `<TLEFromCelestrak>` tag the latest TLE is force-downloaded from
CelesTrak before every run. The satellite is identified by its NORAD catalog number
(all digits) or by name (a name query may match several satellites — all are loaded):
```
<Constellation>
    <NumOfPlanes>1</NumOfPlanes>
    <NumOfSatellites>1</NumOfSatellites>
    <ConstellationID>1</ConstellationID>
    <ConstellationName>TerraSAR-X</ConstellationName>
    <TLEFromCelestrak>31698</TLEFromCelestrak>   <!-- or a name, e.g. TERRASAR-X -->
    <TLEFileName>../input/terrasarx.txt</TLEFileName>
</Constellation>
```
- When `<TLEFileName>` is also given, the downloaded TLE overwrites that file (which
  then doubles as the offline copy). Without it, the TLE is cached as
  `tle_celestrak_<identifier>.txt` next to the Config.xml.
- If the download fails (offline, unknown satellite), a warning is logged and the
  previously downloaded file is used instead; the run only stops when no file
  exists at all.

### Ground segment
The following xml is used to setup the ground segment:
```
  <GroundSegment>
    <Network>
      <NumStation>2</NumStation>
      <NetworkName>GCC</NetworkName>
      <GroundStation>
        <Type>GCC</Type>
        <ConstellationID>1</ConstellationID>
        <GroundStationID>1</GroundStationID>
        <GroundStationName>OBE</GroundStationName>
        <Latitude>48.0744</Latitude>
        <Longitude>11.262</Longitude>
        <Height>0</Height>
        <ReceiverConstellation>110</ReceiverConstellation>
        <ElevationMask>5</ElevationMask>
      </GroundStation>
    </Network>
  </GroundSegment>
```
Ground stations belong to a constellation <ConstellationID> and can be told to receive from multiple constellations 
<ReceiverConstellation> through a list of True/False separated by commas. The elevation mask is defined as a list of
values seperated by commas, dividing the azimuth circle in equal parts.

### User segment
The following xml is used to setup the user segment, here a quadrangle grid:
```
<UserSegment>
    <User>
        <Type>Grid</Type>
        <LatMin>-90</LatMin>
        <LatMax>90</LatMax>
        <LonMin>-180</LonMin>
        <LonMax>180</LonMax>
        <LatStep>10</LatStep>
        <LonStep>10</LonStep>
        <Height>0</Height>
        <ReceiverConstellation>111</ReceiverConstellation>
        <ElevationMask>5</ElevationMask>
    </User>
</UserSegment>
```
or as single static location:
```
<User>
    <Type>Static</Type>
    <Latitude>50</Latitude>
    <Longitude>5</Longitude>
    <Height>0</Height>
    <ReceiverConstellation>111</ReceiverConstellation>
    <ElevationMask>5</ElevationMask>
</User>
```
or as a polygon surface with manual points as tuples:
```
<User>
    <Type>Polygon</Type>
    <Name>Europe</Name>
    <LatStep>.5</LatStep>
    <LonStep>.5</LonStep>
    <PolygonList>(-24.11, 65.34),(-15.97, 52.45),(-14.89, 42.52),(-15.59, 37.36),(-20.89, 42.52),(-30.06, 43.62),(-36.40, 41.78),(-38.86, 35.15),(-27.21, 31.10),(-20.50, 23.74),(-7.47, 33.68),(6.98, 34.42),(19.68, 31.10),(30.26, 33.68),(39.77, 37.73),(44.35, 40.67),(38.70, 45.09),(37.63, 49.51),(27.04, 57.61),(33.37, 63.87),(39.71, 65.71),(36.88, 72.70),(31.58, 78.22),(13.60, 75.64),(1.28, 68.65),(-9.29, 62.02),(-14.59, 66.81)</PolygonList>
    <Height>0</Height>
    <ReceiverConstellation>1111</ReceiverConstellation>
    <ElevationMask>0</ElevationMask>
</User>
```
or as a polygon surface defined in a shapefile:
```
<User>
    <Type>Polygon</Type>
    <Name>Europe</Name>
    <LatStep>.5</LatStep>
    <LonStep>.5</LonStep>
    <PolygonFile>../input/polygon.shp</PolygonFile>
    <Height>0</Height>
    <ReceiverConstellation>1111</ReceiverConstellation>
    <ElevationMask>0</ElevationMask>
</User>
```
or as a spacecraft user through a TLE file:
```
<UserSegment>
    <Type>Spacecraft</Type>
    <TLEFileName>../input/example_tle_files/TLE_MetopA_2006_12_26.txt</TLEFileName>
    <ElevationMask>20</ElevationMask>
    <ReceiverConstellation>1000</ReceiverConstellation>
</UserSegment>
```
Users can be told to receive from multiple constellations <ReceiverConstellation> through a list of True/False 
separated by commas. The elevation mask is defined as a list of values (one or more) seperated by commas, 
dividing the azimuth circle in equal parts.

### Simulation parameters
The following xml is used to setup the simulation parameters:
```
<SimulationManager>
    <StartDate>2013-05-08 00:00:00</StartDate>
    <StopDate>2013-05-09 00:00:00</StopDate>
    <TimeStep>3600</TimeStep>
    <IncludeStation2SpaceLinks>True</IncludeStation2SpaceLinks>
    <IncludeUser2SpaceLinks>True</IncludeUser2SpaceLinks>
    <IncludeSpace2SpaceLinks>False</IncludeSpace2SpaceLinks>
    <OrbitsFromPreviousRun>False</OrbitsFromPreviousRun>
    <OrbitPropagator>SGP4</OrbitPropagator>
    <Report>True</Report>

    <Analysis>
          <Type>cov_ground_track</Type>
          <ConstellationID>1</ConstellationID>
    </Analysis>
    <Analysis>
          <Type>cov_satellite_visible</Type>
    </Analysis>
</SimulationManager>
```
The following explanations apply for the parameters:
- The Start/Stop time parameters are in UTC time and TimeStep in seconds. 
- The IncludeStation2SpaceLinks, etc. parameters determine whether links between different objects: sat, station and user are computed. Normally leave these to True so that all analysis works. Time could
be saved by disabling some. 
- The OrbitsFromPreviousRun flag (True/False) reuses the satellite ECI orbits cached in 'output/orbits_internal.txt' from a previous run instead of re-propagating, to save time when only the analysis changes.
- The OrbitPropagator determines which propagator to take: 'Keplerian', 'SGP4' or 'HPOP'.
- The optional Report flag (True/False, default False) writes a mission report
  __report.html__ into the output directory at the end of the run: one page with a
  scenario configuration summary (simulation parameters, force model,
  constellations and their satellite orbits, ground stations, users, analyses),
  the warnings, and per analysis the summary lines from the log, every plot and
  links to the CSV data dumps, pass tables and movies. The same page can be built
  later from any existing results directory (e.g. a projects/ mission folder)
  without re-running the simulation:
  `py src/report.py <results_dir> [--title "Mission name"] [--config Config.xml]`.
- One or more Analysis blocks can be given. All analyses run in the same simulation
  over the same propagated orbits, each keeping its own metric memory for the
  post-processing/plots in its own output files (named after the analysis Type).
  When the same Type appears more than once (e.g. the same analysis for two
  different satellites), the output files of the repeats are numbered
  __<type>_2.png__, __<type>_3.png__, ... Note that fixed-name data dumps
  (e.g. orbits.txt of cov_satellite_pvt, user_cov_swath.nc of the obs analyses)
  keep their name, so repeats of those analyses overwrite each other's dump.

### HPOP orbit propagator
'HPOP' selects the High Precision Orbit Propagation based on the Orekit library
(numerical integration with a configurable perturbation force model, see
src/propagation_hpop.py). It requires the python packages __orekit_jpype__ and
__jdk4py__ and the Orekit physical data archive at __input/orekit-data.zip__.
Satellites defined by a TLE file get their initial state from the TLE at the
simulation start; satellites defined by Kepler elements are integrated from their
EpochMJD (keep EpochMJD close to StartDate). Each satellite is integrated once over
the whole simulation window; the time loop samples the resulting dense ephemeris.

When (and only when) HPOP is selected, an additional <HPOP> block inside
SimulationManager configures the force model. Every entry is optional — the
defaults give the full force model shown here:
```
<HPOP>
    <IntegratorMinStep>0.001</IntegratorMinStep>
    <IntegratorMaxStep>300</IntegratorMaxStep>
    <IntegratorPositionTolerance>1.0</IntegratorPositionTolerance>
    <Mass>1000.0</Mass>

    <Geopotential>True</Geopotential>
    <GeopotentialDegree>21</GeopotentialDegree>
    <GeopotentialOrder>21</GeopotentialOrder>

    <EarthPoleRotation>True</EarthPoleRotation>

    <Drag>True</Drag>
    <DragArea>1.0</DragArea>
    <DragCd>2.2</DragCd>
    <DragModel>NRLMSISE00</DragModel>

    <SolarRadiationPressure>True</SolarRadiationPressure>
    <SRPArea>1.0</SRPArea>
    <SRPCr>1.5</SRPCr>

    <ThirdBodySun>True</ThirdBodySun>
    <ThirdBodyMoon>True</ThirdBodyMoon>
    <ThirdBodyPlanets>False</ThirdBodyPlanets>

    <SolidTides>True</SolidTides>
    <OceanTides>False</OceanTides>
    <OceanTidesDegree>4</OceanTidesDegree>
    <OceanTidesOrder>4</OceanTidesOrder>

    <Relativity>False</Relativity>
</HPOP>
```
The parameters are:
- IntegratorMinStep/IntegratorMaxStep: variable step size bounds in seconds of the
  Dormand-Prince 8(5,3) integrator.
- IntegratorPositionTolerance: position error tolerance in m controlling the step size.
- Mass: spacecraft mass in kg (used by drag and SRP accelerations).
- Geopotential: spherical harmonics gravity field (Holmes-Featherstone) with
  GeopotentialDegree x GeopotentialOrder resolution; the point-mass central
  attraction is always applied.
- EarthPoleRotation: True uses the full IERS 2010 Earth orientation (precession,
  nutation, UT1 and polar motion from the EOP files) for the Earth-fixed frame of
  the force models; False uses simplified equinox-based transforms without the
  pole corrections.
- Drag: atmospheric drag with DragArea (m^2), DragCd (drag coefficient) and
  DragModel one of 'NRLMSISE00', 'DTM2000' (both driven by the CSSI space weather
  file in orekit-data) or 'HarrisPriester' (static density model).
- SolarRadiationPressure: isotropic SRP with SRPArea (m^2) and reflection
  coefficient SRPCr, including Earth shadow.
- ThirdBodySun/ThirdBodyMoon/ThirdBodyPlanets: point mass third-body attraction
  (planets adds Venus, Mars and Jupiter).
- SolidTides/OceanTides: tidal corrections to the gravity field (ocean tides up to
  OceanTidesDegree x OceanTidesOrder).
- Relativity: Schwarzschild relativistic correction.

A benchmark of the HPOP propagator against an analytical two-body orbit (cm-level
agreement) and against an SGP4 reference trajectory (km-level agreement) is
available in tests/hpop_benchmark (run tests/hpop_benchmark/benchmark_hpop.py).

The analysis are described below:

## Analysis parameters
In order to run an analysis block it has to be uncommented and the parameters adapted. To add a new analysis to the code the following has to be performed:
- Add a new analysis class at the end of analysis.py
- Use as a template one of the other analysis classes above or the base analysis class definition
- Add the class instantiation at the end of config.py

### 2D map options
Every analysis that produces a 2D world map accepts these optional flags in its
`<Analysis>` block to decorate the map (all default off):
```
<EarthImage>True</EarthImage>
<ShowStations>True</ShowStations>
<ShowUsers>True</ShowUsers>
<ShowGroundTrack>True</ShowGroundTrack>
<Coastlines>False</Coastlines>
```
- EarthImage: draw the map background as an Earth photo (input/earth_texture.jpg,
  the same texture as the 3D globe, with the cartopy stock image as fallback)
  instead of plain coastlines.
- ShowStations: mark the ground station locations with red triangles (drawn
  underneath the ground tracks) and halo-outlined name labels.
- ShowUsers: mark the user locations with blue triangles (skipped with a warning
  for user grids with more than 200 points).
- ShowGroundTrack: draw the subsatellite ground track(s) of all satellites as thin
  grey lines, recorded during the run.
- Coastlines: set to False to turn the coastline outlines off (default on) — e.g.
  with EarthImage the photo already shows the coasts.

For example the cov_ground_track reference image below uses EarthImage,
ShowStations and ShowUsers; the obs_target_imaging image uses ShowGroundTrack.

### 3D plots
All analyses that produce a world map (cov_ground_track, cov_satellite_contour,
cov_satellite_visible_grid, cov_satellite_highest, cov_depth_of_coverage,
cov_pass_time, obs_swath_conical, obs_swath_push_broom, obs_sza_push_broom,
obs_sza_subsat, nav_dilution_of_precision, nav_accuracy) can additionally render
their result on a textured 3D Earth (needs pyvista and input/earth_texture.jpg,
see the installation section). The 3D plot is saved as
__output/<analysis_type>_3d.png__; when the tool runs with a display (not
headless) an interactive pyvista window opens first and the screenshot is saved
on closing it. The following parameters in the Analysis block control it:
```
<Plot3D>True</Plot3D>
<ShowSatellite>True</ShowSatellite>
<ShowOrbit>True</ShowOrbit>
<SatelliteModelFile>../input/my_satellite.stl</SatelliteModelFile>
<SatelliteModelScale>200000</SatelliteModelScale>
<EarthClouds>False</EarthClouds>
<MP4>False</MP4>
```
- Plot3D: enables the 3D plot (default False).
- ShowSatellite: draw a 3D satellite model at the last simulated epoch of every
  satellite, nadir pointing with the solar panel axis cross-track and a yellow
  line to the subsatellite point (default True). Without SatelliteModelFile a
  simple procedural bus + solar panels + antenna dish model is drawn.
- ShowOrbit: draw the satellite orbital track(s) as cyan lines (default True).
  The track is the inertial orbit path over the simulation, oriented as the
  Earth stands at the final epoch — i.e. the familiar orbit ellipse(s) (the
  path relative to the rotating Earth is what the ground track shows).
- SatelliteModelFile: optional satellite mesh file (STL/OBJ/PLY/VTK, +x flight
  direction, +y solar panel axis, +z towards nadir).
- SatelliteModelScale: size of the (hugely exaggerated) satellite model in
  metres, default 200000, so that it is visible at globe scale. For MEO
  constellations like GPS a larger value (e.g. 500000) is recommended.
- EarthClouds: adds a semi-transparent cloud layer (input/earth_clouds.jpg) on
  top of the Blue Marble globe (default False). The 3D scene always shows the
  starry Milky Way background when input/starmap.jpg is present.
- MP4: writes movies of the world maps (default False, needs imageio +
  imageio-ffmpeg). __output/<analysis_type>_2d.mp4__ shows the 2D map filling
  up over the simulation time (growing ground track, growing swath ribbons, or
  the instantaneous metric field per epoch on the user grid);
  __output/<analysis_type>_3d.mp4__ is a fly-along of the 3D scene: the camera
  circles once around the (first) satellite while it moves along its orbit,
  with the Earth in view and the analysis data growing on the globe. Supported
  by cov_ground_track, cov_depth_of_coverage, cov_pass_time,
  cov_satellite_visible_grid, cov_satellite_highest, nav_dilution_of_precision,
  nav_accuracy, obs_swath_conical and obs_swath_push_broom. Movies are capped
  at 360 frames (longer simulations are subsampled) at 20 frames per second.
Ground stations are always drawn as magenta markers. Near-Earth scenes use a
perspective camera above the (last) satellite; scenes whose orbits reach far
above the Earth (MEO/GEO, e.g. GNSS constellations) are rendered with a fitted
parallel projection instead, so the Earth and the orbits appear at their true
relative size.

The specific parameters for the existing analysis are given here below:

### cov_ground_track
Plots the ground track of one or more satellites over simulation time as a
continuous line (broken cleanly at the date line, matching the 3D view). The
following parameters are needed, to plot the ground track of satellites in a
constellation:
```
<Type>cov_ground_track</Type>
<ConstellationID>1</ConstellationID>
```
Optional is SatelliteID for selection of one satellite. For a TLE file it is the NORAD number. 
```
<SatelliteID>1</SatelliteID>
```
If this parameter is omitted all the satellites of the constellation are plotted.
<img src="/docs/cov_ground_track.png" alt="cov_ground_track"/>

Optionally the ground track can additionally be rendered in 3D with the red ground
track on the surface (see the 3D plots section above for all Plot3D parameters):
```
<Plot3D>True</Plot3D>
```
<img src="/docs/cov_ground_track_3d.png" alt="cov_ground_track_3d"/>

### cov_satellite_pvt
Plots the satellite position and velocity and outputs the position and velocity to the 'output/orbits.txt' file. 
The following parameters are needed:
```
<Type>cov_satellite_pvt</Type>
<ConstellationID>1</ConstellationID>
```
Optional is:
```
<SatelliteID>1</SatelliteID>
``` 
For a TLE file the SatelliteID is the NORAD number. If this parameter is omitted all the satellites of the 
constellation are output to file and the plot shows the first recorded satellite. Additionally an 
__/output/orbits.txt__ file is saved to disk.
<img src="/docs/cov_satellite_pvt.png" alt="cov_satellite_pvt"/>

### cov_satellite_visible
Plots the number of available satellites for the user(s). The following parameters are needed:
```
<Type>cov_satellite_visible</Type>
```
<img src="/docs/cov_satellite_visible.png" alt="cov_satellite_visible"/>

### cov_satellite_visible_grid
Plots the number of available satellites at a user grid. Statistics on plots can be: minimum, mean, maximum, std and
median. The following parameters are needed:
```
<Analysis>
      <Type>cov_satellite_visible_grid</Type>
      <Statistic>Min</Statistic>
</Analysis>
```
<img src="/docs/cov_satellite_visible_grid.png" alt="cov_satellite_visible_grid"/>

With <Plot3D>True</Plot3D> the statistic is also draped over the 3D globe (see the
3D plots section for all parameters; the example below uses <ShowOrbit>False</ShowOrbit>):
<img src="/docs/cov_satellite_visible_grid_3d.png" alt="cov_satellite_visible_grid_3d"/>

### cov_satellite_visible_id
Plots the satellite IDs in view over time for the first user. The following parameters are needed:
```
<Analysis>
    <Type>cov_satellite_visible_id</Type>
    <ConstellationID>1</ConstellationID>
</Analysis>
```
<img src="/docs/cov_satellite_visible_id.png" alt="cov_satellite_visible_id"/>

### cov_satellite_contour
Plots the satellite(s) ground contour on the world map. The following parameters are needed:
```
<Analysis>
    <Type>cov_satellite_contour</Type>
    <ConstellationID>1</ConstellationID>
    <ElevationMask>20</ElevationMask>
</Analysis>
```
Optional is:
```
    <SatelliteID>1</SatelliteID>
``` 
The elevation mask is for the user who has to receive the satellite. The satellite is selected by constellation ID and
satellite ID or multiple if satellite ID is omitted.
<img src="/docs/cov_satellite_contour.png" alt="cov_satellite_contour"/>

With <Plot3D>True</Plot3D> the visibility contour(s) at the end of the simulation are
drawn as coloured rings on the 3D globe (see the 3D plots section for all parameters):
<img src="/docs/cov_satellite_contour_3d.png" alt="cov_satellite_contour_3d"/>

### cov_satellite_sky_angles
Plots the satellite azimuth and elevation over time for the first user. The following parameters are needed:
```
<Analysis>
    <Type>cov_satellite_sky_angles</Type>
    <ConstellationID>3</ConstellationID>
    <SatelliteID>24307</SatelliteID>
</Analysis>
```
The satellite is selected by constellation ID and satellite ID.
<img src="/docs/cov_satellite_sky_angles.png" alt="cov_satellite_sky_angles"/>

### cov_depth_of_coverage
Plots the number of ground stations in view from the satellite over the orbit of the satellites. 
The following parameters are needed:
```
<Analysis>
  <Type>cov_depth_of_coverage</Type>
</Analysis>
```
The elevation mask is taken by the ground station setup.
<img src="/docs/cov_depth_of_coverage.png" alt="cov_depth_of_coverage"/>

With <Plot3D>True</Plot3D> the ground track is drawn on the 3D globe with the points
coloured by the number of stations in view (see the 3D plots section for all parameters):
<img src="/docs/cov_depth_of_coverage_3d.png" alt="cov_depth_of_coverage_3d"/>

### cov_pass_time
Plots the satellite constellation pass time statistics for a user grid. The following parameters are needed:
```
<Analysis>
    <Type>cov_pass_time</Type>
    <ConstellationID>1</ConstellationID>
    <Statistic>Mean</Statistic>
</Analysis>
```
<img src="/docs/cov_pass_time.png" alt="cov_pass_time"/>

With <Plot3D>True</Plot3D> the statistic is draped over the 3D globe (the example
below uses <ShowSatellite>False</ShowSatellite>, showing only the constellation
orbital tracks; see the 3D plots section for all parameters):
<img src="/docs/cov_pass_time_3d.png" alt="cov_pass_time_3d"/>

### cov_satellite_highest
Plots elevation for the highest satellite in view over a user grid. The following parameters are needed:
```
<Analysis>
    <Type>cov_satellite_highest</Type>
    <ConstellationID>1</ConstellationID>
    <Statistic>Mean</Statistic>
</Analysis>
```
<img src="/docs/cov_satellite_highest.png" alt="cov_satellite_highest"/>

With <Plot3D>True</Plot3D> the statistic is draped over the 3D globe together with
the constellation (see the 3D plots section for all parameters):
<img src="/docs/cov_satellite_highest_3d.png" alt="cov_satellite_highest_3d"/>

### obs_swath_conical
Plots the swath coverage for a conical scanner on one or more satellites defined in the space segment.
The swath is drawn as smooth semi-transparent ribbons built from the swath edge points of every pass
(overlapping passes show darker), matching the 3D globe render.
The user segment grid defines the granularity of the revisit statistics and the SaveOutput data
(not of the swath coverage map itself).
Typically a grid of 1x1 deg is sufficient otherwise for a complete globe the simulation will take lots of time.

The following parameters are needed:
```
<Analysis>
    <Type>obs_swath_conical</Type>
</Analysis>
```
In the constellation part of the space segment are defined the instrument characteristics:
```
<ObsSwathStop>650000.0</ObsSwathStop>
```
as above in meters, at the edge, or in degrees incidence angle:
```
<ObsIncidenceAngleStop>52.0</ObsIncidenceAngleStop>
```
The incidence angle is defined as the angle between the line-of-sight and the nadir vector from the satellite. 
This is not to be confused with the user observation zenith angle.

Optional in the analysis part are:
```
<PolarView>90</PolarView>
<Revisit>True</Revisit>
<Statistic>Mean</Statistic>
<SaveOutput>Numpy</SaveOutput>
```
- PolarView angle: This parameter can be given to see one part of the globe in an stereographic view, eg.  for the polar region.
- Revisit flag: This flag will enable revisit computation after the swath coverage. The statistic will determine
what kind of statistic is displayed per user location. Besides the revisit map, a longitude-averaged
max/mean revisit time versus latitude profile is saved as __output/<analysis_type>_revisit_lat.png__.
- SaveOutput: NetCDF or Numpy This flag will enable saving user swath coverage for every timestep.
- Plot3D flag: additionally renders the swath as a semi-transparent ribbon on a
  textured 3D Earth, saved as output/obs_swath_conical_3d.png (see the 3D plots
  section for all Plot3D parameters).

<img src="/docs/obs_swath_conical.png" alt="obs_swath_conical"/>
<img src="/docs/obs_swath_conical_revisit.png" alt="obs_swath_conical_revisit"/>
<img src="/docs/obs_swath_conical_revisit_lat.png" alt="obs_swath_conical_revisit_lat"/>
<img src="/docs/obs_swath_conical_3d.png" alt="obs_swath_conical_3d"/>

### obs_swath_push_broom
Plots the swath coverage for a push broom scanner on one or more satellites defined in the space segment.
The swath is drawn as smooth semi-transparent ribbons built from the swath edge points of every pass
(overlapping passes show darker), matching the 3D globe render.
The user segment grid defines the granularity of the revisit statistics and the SaveOutput data
(not of the swath coverage map itself).
Typically a grid of 1x1 deg is sufficient otherwise for a complete globe the simulation will take lots of time.

The following parameters are needed:
```
<Analysis>
    <Type>obs_swath_push_broom</Type>
</Analysis>
```
In the constellation part of the space segment are defined the instrument characteristics:
```
<ObsSwathStart>250000.0</ObsSwathStart>
<ObsSwathStop>650000.0</ObsSwathStop>
```
as above in meters, at the edge, or in degrees incidence angle:
```
<ObsIncidenceAngleStart>20.0</ObsIncidenceAngleStart>
<ObsIncidenceAngleStop>52.0</ObsIncidenceAngleStop>
```
The incidence angle is defined as the angle between the line-of-sight and the nadir vector from the satellite. 
This is not to be confused with the user observation zenith angle.
Positive values are for the right looking situation. 
Start is left most point, Stop is right most point, looking in direction of velocity vector.

Optional in the analysis part are:
```
<PolarView>60</PolarView>
<Revisit>True</Revisit>
<Statistic>Mean</Statistic>
<SaveOutput>NetCDF</SaveOutput>
```
- PolarView angle: This parameter can be given to see one part of the globe in an stereographic view, eg.  for the polar region.
  The value describes the minimum bounding latitude visible. When negative the area on the South Pole will be visible.
- Revisit flag: This flag will enable revisit computation after the swath coverage. The statistic will determine
what kind of statistic is displayed per user location. Besides the revisit map, a longitude-averaged
max/mean revisit time versus latitude profile is saved as __output/<analysis_type>_revisit_lat.png__.
- SaveOutput: NetCDF or Numpy This flag will enable saving user swath coverage for every timestep.
- Plot3D flag: additionally renders the swath as a semi-transparent ribbon on a
  textured 3D Earth, saved as output/obs_swath_push_broom_3d.png (see the 3D plots
  section for all Plot3D parameters).

<img src="/docs/obs_swath_push_broom.png" alt="cov_satellite_push_broom"/>
<img src="/docs/obs_swath_push_broom_revisit.png" alt="cov_satellite_push_broom_revisit"/>
<img src="/docs/obs_swath_push_broom_revisit_lat.png" alt="cov_satellite_push_broom_revisit_lat"/>
<img src="/docs/obs_swath_push_broom_3d.png" alt="obs_swath_push_broom_3d"/>


### obs_sza_push_broom
Plots the mean Solar Zenith Angle (SZA) for the user locations falling within the push broom 
swath of one or more satellites. It combines the push broom swath geometry (as in 
obs_swath_push_broom) with a per-location solar zenith angle computation, so it shows the 
illumination conditions under which each location is observed.
The user segment is used to define the grid of analysis and defines the granularity of the result.

_Note: this analysis is in an early stage and the solar angle computation makes it slow; use a 
coarse user grid._

The following parameters are needed:
```
<Analysis>
    <Type>obs_sza_push_broom</Type>
</Analysis>
```
In the constellation part of the space segment the instrument characteristics are defined the same
way as for obs_swath_push_broom, either in meters at the edge:
```
<ObsSwathStart>250000.0</ObsSwathStart>
<ObsSwathStop>650000.0</ObsSwathStop>
```
or in degrees incidence angle:
```
<ObsIncidenceAngleStart>20.0</ObsIncidenceAngleStart>
<ObsIncidenceAngleStop>52.0</ObsIncidenceAngleStop>
```

Optional in the analysis part are:
```
<PolarView>60</PolarView>
<SaveOutput>NetCDF</SaveOutput>
<Plot3D>True</Plot3D>
```
- PolarView angle: This parameter can be given to see one part of the globe in an stereographic view, eg. for the polar region.
  When negative the area on the South Pole will be visible.
- SaveOutput: NetCDF or Numpy This flag will enable saving the user swath coverage for every timestep.
- Plot3D flag: additionally renders the mean SZA of the observed locations as
  coloured points on the 3D globe (see the 3D plots section for all parameters):

<img src="/docs/obs_sza_push_broom_3d.png" alt="obs_sza_push_broom_3d"/>


### obs_sza_subsat
Plots the Solar Zenith Angle SZA for the subsatellite point. 
All the satellites defined in the config will be used.

The following parameters are needed:
```
<Analysis>
    <Type>obs_sza_subsat</Type>
</Analysis>
```
<img src="/docs/obs_sza_subsat.png" alt="obs_sza_subsat"/>

Optional in the analysis part are:
```
<PolarView>90</PolarView>
<RangeLatitude>-80,80,10</RangeLatitude>
<SaveOutput>Numpy</SaveOutput>
```
- PolarView angle: This parameter can be given to see one part of the globe in an stereographic view, eg.  for the polar region.
  When negative the area on the South Pole will be visible.
- RangeLatitude: Min_Lat, Max_Lat and Lat_Step in degrees, will enable plots vs. latitude.
- SaveOutput: Will write to file the SZA vs latitude averaged over the simulation time.
- Plot3D (True/False): additionally renders the subsatellite SZA points on the 3D
  globe (see the 3D plots section for all parameters):

<img src="/docs/obs_sza_subsat_3d.png" alt="obs_sza_subsat_3d"/>

<img src="/docs/obs_sza_subsat_lat.png" alt="obs_sza_subsat_lat"/>
<img src="/docs/obs_sza_subsat_lat_year.png" alt="obs_sza_subsat_lat_year"/>    


### obs_aoi_revisit
Revisit and coverage build-up statistics over an area of interest. The AOI is the
configured user segment — typically a Type Polygon user (a grid clipped to an inline
`<PolygonList>` or a `<PolygonFile>` shapefile, see the user segment section), but a
regional Grid works too. The instrument is the push-broom swath defined in the
`<Constellation>` block (ObsSwathStart/Stop or ObsIncidenceAngleStart/Stop), evaluated
with the same machinery as obs_swath_push_broom. Produces a map of the revisit
statistic per AOI grid point zoomed to the AOI, the fraction of the AOI covered at
least once versus time (fill-up curve), and aggregate numbers in the log (coverage
percentage, time to 50/90/99% coverage, mean and worst revisit gap).

The following parameters are needed:
```
<Analysis>
    <Type>obs_aoi_revisit</Type>
</Analysis>
```

Optional in the analysis part is:
```
    <Statistic>Max</Statistic>
```
- Statistic: min / mean / max / std / median revisit interval shown on the map (default max).

<img src="/docs/obs_aoi_revisit.png" alt="obs_aoi_revisit"/>
<img src="/docs/obs_aoi_revisit_coverage.png" alt="obs_aoi_revisit_coverage"/>


### obs_target_imaging
Imaging opportunities over a list of point targets for an agile satellite: a target
can be imaged when it lies within MaxOffNadir degrees of the satellite nadir direction
(the pointing agility cone) and, optionally, while the Sun is at least MinSunElevation
degrees above the target horizon (optical imaging daylight constraint). The log reports
the opportunity windows per target (start, duration, best off-nadir angle) and the
per-target counts and revisit gaps; the plot shows the targets coloured by opportunity
count. The CSV dumps contain the opportunity table and the per-target statistics.

The following parameters are needed:
```
<Analysis>
    <Type>obs_target_imaging</Type>
    <MaxOffNadir>45</MaxOffNadir>
    <Target>Rome, 41.9, 12.5</Target>
    <Target>Delft, 52.0, 4.36</Target>
</Analysis>
```
Parameters are:
- MaxOffNadir: pointing agility cone half-angle in degrees.
- Target: one per target, "Name, latitude_deg, longitude_deg" (repeat the tag per target).

Optional in the analysis part are:
```
    <MinSunElevation>10</MinSunElevation>
    <ConstellationID>1</ConstellationID>
    <TargetFile>../input/targets.csv</TargetFile>
```
- MinSunElevation: minimum Sun elevation at the target in degrees (default: no
  daylight constraint).
- ConstellationID: restrict the imaging satellites to one constellation (default all).
- TargetFile: CSV file with one "Name, lat_deg, lon_deg" line per target, as an
  alternative (or in addition) to inline `<Target>` tags.

<img src="/docs/obs_target_imaging.png" alt="obs_target_imaging"/>


### Antenna patterns for the com analyses (GRASP .cut/.grd)
The link-budget analyses (com_gr2sp_budget, com_sp2sp_budget,
com_gr2sp_budget_interference) can use real antenna patterns instead of the
scalar gains, read from the TICRA GRASP file formats:
```
<TransmitAntennaPatternFile>../input/example_antenna_patterns/isoflux_x_band.cut</TransmitAntennaPatternFile>
<ReceiveAntennaPatternFile>../input/example_antenna_patterns/dish_x_band_3m.grd</ReceiveAntennaPatternFile>
```
- Tabulated cut files (__.cut__, polar cuts ICUT=1) and ASCII grid files
  (__.grd__, uv grids IGRID=1) are supported; the far field must use the
  standard GRASP power normalisation (gain dBi = 10log10 of the summed squared
  field components). The 2D pattern is power-averaged over the phi cuts /
  grid directions into a gain versus off-boresight-angle curve.
- When a pattern file is given, the corresponding scalar gain
  (TransmitGaindB/ReceiveGaindB) is replaced. The ground station antenna is
  assumed to track the satellite (peak gain). In com_gr2sp_budget the
  satellite antenna is nadir-pointed by default — the gain follows the
  off-nadir angle to the station every epoch, which is the right model for
  isoflux and horn antennas; set
  `<SatelliteAntennaPointing>Tracking</SatelliteAntennaPointing>` for a
  steered antenna (peak gain). In com_sp2sp_budget both ISL antennas track
  each other (peak gains). In com_gr2sp_budget_interference the pattern files
  replace the analytic dish patterns: peak gain on the nominal link, and the
  interferer discriminated by the leader/interferer separation angle through
  the actual pattern sidelobes.
- The loaded pattern(s) are plotted to __output/<analysis_type>_antenna.png__
  for verification.

Ready-to-use example patterns (LEO X-band isoflux, X-band 3 m and Ka-band
0.25 m / 6.8 m dishes) are provided in __input/example_antenna_patterns__,
generated by make_patterns.py there from realistic aperture models — patterns
exported from GRASP or compatible antenna tools can be dropped in directly.

### com_gr2sp_budget
Plots the link budget parameters for a certain ground station to all satellites. 
The models used are coming from ITU-R python package itur, which implements:
- ITU-R P.453-13: The radio refractive index: its formula and refractivity data
- ITU-R P.618-13: Propagation data and prediction methods required for the design of Earth-space telecommunication systems
- ITU-R P.676-11: Attenuation by atmospheric gases 
- ITU-R P.835-6: Reference Standard Atmospheres 
- ITU-R P.836-6: Water vapour: surface density and total columnar content 
- ITU-R P.837-7: Characteristics of precipitation for propagation modelling 
- ITU-R P.838-3: Specific attenuation model for rain for use in prediction methods 
- ITU-R P.839-4: Rain height model for prediction methods. 
- ITU-R P.840-7: Attenuation due to clouds and fog
- ITU-R P.1144-7: Interpolation methods for the geophysical properties used to compute propagation effects 
- ITU-R P.1511-1: Topography for Earth-to-space propagation modelling 
- ITU-R P.1853-1: Tropospheric attenuation time series synthesis

The following parameters are needed:
```
<Analysis>
    <Type>com_gr2sp_budget</Type>
    <GroundStationID>1</GroundStationID>
    <TransmitterObject>Satellite</TransmitterObject>
    <CarrierFrequency>10e9</CarrierFrequency>
    <TransmitPowerW>10</TransmitPowerW>
    <TransmitLossesdB>2</TransmitLossesdB>
    <TransmitGaindB>20</TransmitGaindB>
    <ReceiveGaindB>64</ReceiveGaindB>
    <ReceiveLossesdB>3</ReceiveLossesdB>
    <ReceiveTempK>290</ReceiveTempK>
    <PExceedPerc>0.01</PExceedPerc>
    <IncludeGas>True</IncludeGas>
    <IncludeRain>True</IncludeRain>
    <IncludeScintillation>False</IncludeScintillation>
    <IncludeClouds>False</IncludeClouds>
</Analysis>
```
Parameters are:
- GroundStationID: Station to be used, refer to the ground segment part.
- TransmitterObject: Satellite or Ground Station, which one is transmitting
- CarrierFrequency: Carrier frequency of signal in Hz
- TransmitPowerW: Transmit power of transmitter in W
- TransmitLossesdB: All transmit losses in dB
- TransmitGaindB: Transmit gain of antenna in dB
- ReceiveGaindB: Receive gain of antenna in dB
- ReceiveLossesdB: All receive losses in dB
- ReceiveTempK: Receive Temperature in K
- PExceedPerc: Probability to exceed attenuation values in %.
- IncludeGas: Whether gas attenuation should be included True/False
- IncludeRain: Whether rain attenuation should be included True/False
- IncludeScintillation: Whether scintillation attenuation should be included True/False
- IncludeClouds: Whether cloud attenuation should be included True/False

Optional in the analysis part are:
```
    <ModulationType>BPSK</ModulationType>
    <BitErrorRate>1e-5</BitErrorRate>
    <DataRateBitPerSec>1e6</DataRateBitPerSec>
    <TransmitAntennaPatternFile>../input/example_antenna_patterns/isoflux_x_band.cut</TransmitAntennaPatternFile>
    <ReceiveAntennaPatternFile>../input/example_antenna_patterns/dish_x_band_3m.grd</ReceiveAntennaPatternFile>
    <SatelliteAntennaPointing>Nadir</SatelliteAntennaPointing>
```

Some parameters can be entered to get the required CN0:
- ModulationType: 'BPSK' or 'QPSK'
- BitErrorRate: bit error rate required
- DataRateBitPerSec: datarate required

The antenna pattern files replace the scalar gains (see the Antenna patterns
section above); SatelliteAntennaPointing is Nadir (default) or Tracking.

<img src="/docs/com_gr2sp_budget.png" alt="com_gr2sp_budget"/>

### com_gr2sp_budget_interference
Plots the link budget parameters for a certain ground station to all satellites as in com_gr2sp_budget,
but now includes interference from a second satellite following closely the nominal link satellite. 

The following parameters are needed:
```
<Analysis>
    <Type>com_gr2sp_budget</Type>
    <GroundStationID>1</GroundStationID>
    <TransmitterObject>Satellite</TransmitterObject>
    <CarrierFrequency>10e9</CarrierFrequency>
    <BandWidth>675e6</BandWidth>
    <TransmitPowerW>10</TransmitPowerW>
    <TransmitLossesdB>2</TransmitLossesdB>
    <TransmitGaindB>20</TransmitGaindB>
    <TransmitAntennaDiameter>0.25</TransmitAntennaDiameter>
    <ReceiveGaindB>64</ReceiveGaindB>
    <ReceiveAntennaDiameter>6</ReceiveAntennaDiameter>
    <ReceiveLossesdB>3</ReceiveLossesdB>
    <ReceiveTempK>290</ReceiveTempK>
    <PExceedPerc>0.5</PExceedPerc>
    <IncludeGas>True</IncludeGas>
    <IncludeRain>True</IncludeRain>
    <IncludeScintillation>False</IncludeScintillation>
    <IncludeClouds>False</IncludeClouds>
</Analysis>
```
Parameters are:
- GroundStationID: Station to be used, refer to the ground segment part.
- TransmitterObject: Satellite or Ground Station, which one is transmitting
- CarrierFrequency: Carrier frequency of signal in Hz
- BandWidth: Signal bandwith in Hz
- TransmitPowerW: Transmit power of transmitter in W
- TransmitLossesdB: All transmit losses in dB
- TransmitGaindB: Transmit gain of antenna in dB (theoretical pattern assumed), or,
  TransmitGainManualdB: list of tuples in string from 0-180 off boresight gain values
- TransmitAntennaDiameter: Antenna diameter in m (not used to compute Gain)
- ReceiveGaindB: Receive gain of antenna in dB
- ReceiveLossesdB: All receive losses in dB
- ReceiveTempK: Receive Temperature in K
- ReceiveAntennaDiameter: Antenna diameter in m (not used to compute Gain)
- PExceedPerc: Probability to exceed attenuation values in %.

Optional in the analysis part are:
```
    <ModulationType>BPSK</ModulationType>
    <BitErrorRate>1e-5</BitErrorRate>
    <DataRateBitPerSec>1e6</DataRateBitPerSec>
```

Some parameters can be entered to get the required CN0:
- ModulationType: 'BPSK' or 'QPSK'
- BitErrorRate: bit error rate required
- DataRateBitPerSec: datarate required

<img src="/docs/com_gr2sp_budget_interference.png" alt="com_gr2sp_budget_interference"/>


### com_sp2sp_budget
Plots the link budget parameters for a certain satellite to another satellite. 

The following parameters are needed:
```
<Analysis>
    <Type>com_sp2sp_budget</Type>
    <SatelliteID1>1</SatelliteID1>
    <SatelliteID2>2</SatelliteID2>
    <CarrierFrequency>10e9</CarrierFrequency>
    <TransmitPowerW>10</TransmitPowerW>
    <TransmitLossesdB>2</TransmitLossesdB>
    <TransmitGaindB>20</TransmitGaindB>
    <ReceiveGaindB>20</ReceiveGaindB>
    <ReceiveLossesdB>2</ReceiveLossesdB>
    <ReceiveTempK>290</ReceiveTempK>
</Analysis>
```
Parameters are:
- SatelliteID1: Satellite to be used as transmitter.
- SatelliteID2: Satellite to be used as receiver.
- CarrierFrequency: Carrier frequency of signal in Hz
- TransmitPowerW: Transmit power of transmitter in W
- TransmitLossesdB: All transmit losses in dB
- TransmitGaindB: Transmit gain of antenna in dB

Optional in the analysis part are:
```
    <ModulationType>BPSK</ModulationType>
    <BitErrorRate>1e-5</BitErrorRate>
    <DataRateBitPerSec>1e6</DataRateBitPerSec>
```

Some parameters can be entered to get the required CN0:
- ModulationType: 'BPSK' or 'QPSK'
- BitErrorRate: bit error rate required
- DataRateBitPerSec: datarate required

<img src="/docs/com_sp2sp_budget.png" alt="com_sp2sp_budget"/>


### com_doppler
Plots the doppler shift in Hz for the station to satellites. 

The following parameters are needed:
```
<Analysis>
    <Type>com_doppler</Type>
    <StationID>1</StationID>
    <CarrierFrequency>10e9</CarrierFrequency>
</Analysis>
```
Parameters are:
- StationID: station to be selected.
- CarrierFrequency: Carrier frequency of signal in Hz

<img src="/docs/com_doppler.png" alt="com_doppler"/>


### com_contact_plan
Ground station contact plan: every station-satellite pass over the simulation
window as a table with AOS, LOS, duration, maximum elevation, downlinkable data
volume and a station-conflict flag when two passes overlap at the same station
(one antenna cannot track two satellites). Three outputs are produced: the CSV
data dump, a human-readable pass table __output/com_contact_plan.txt__ with UTC
times, and a pass-timeline plot (one lane per station with a sub-lane per
satellite; overlapping passes are edged in red). The log summarises the number
of passes, contact minutes per day, mean/max pass duration, data volume per day
and the number of conflicts per station.

All parameters are optional:
```
<Analysis>
    <Type>com_contact_plan</Type>
    <GroundStationID>0</GroundStationID>
    <ConstellationID>0</ConstellationID>
    <MinDuration>60</MinDuration>
    <DownlinkRateMbps>1070</DownlinkRateMbps>
</Analysis>
```
- GroundStationID/ConstellationID: restrict the plan to one station and/or one
  constellation (0 or omitted: all).
- MinDuration: passes shorter than this many seconds are dropped (default 0).
- DownlinkRateMbps: when given, each pass carries its downlinkable volume
  (duration x rate) in Gbit and the log reports Gbit/day per station.

<img src="/docs/com_contact_plan.png" alt="com_contact_plan"/>

### com_pfd
Power flux density produced at the ground by the satellite emissions, versus
elevation, against an ITU-R Article 21 style limit mask — the classic
regulatory compliance check of a downlink. Per epoch the PFD at the ground
station is the EIRP spectral density (transmit power spread over BandWidth,
taken in the ReferenceBandwidth) plus the satellite antenna gain at the
epoch's off-nadir angle (nadir-pointed GRASP pattern file, or a fixed gain)
minus the spreading loss 10log10(4 pi d^2); with several satellites in view
the worst margin is kept. The limit mask follows the standard Article 21
shape: PfdLimit up to 5 deg elevation, +0.5 dB per degree between 5 and 25
deg, PfdLimit+10 above 25 deg. Atmospheric attenuation is not credited
(conservative, as in the Radio Regulations). The log reports the worst margin
and the elevation where it occurs.

The following parameters are needed:
```
<Analysis>
    <Type>com_pfd</Type>
    <CarrierFrequency>8025e6</CarrierFrequency>
    <TransmitPowerW>20</TransmitPowerW>
    <BandWidth>300e6</BandWidth>
    <TransmitAntennaPatternFile>../input/example_antenna_patterns/isoflux_x_band.cut</TransmitAntennaPatternFile>
</Analysis>
```
- CarrierFrequency: carrier frequency in Hz (only used for reporting).
- TransmitPowerW: transmit power in W.
- BandWidth: occupied bandwidth in Hz the transmit power is spread over.
- TransmitAntennaPatternFile or TransmitGaindB: satellite antenna gain, either
  a GRASP .cut/.grd pattern evaluated at the off-nadir angle or a fixed dBi value.

Optional are:
```
    <GroundStationID>1</GroundStationID>
    <TransmitLossesdB>1</TransmitLossesdB>
    <ReferenceBandwidth>4000</ReferenceBandwidth>
    <PfdLimit>-150</PfdLimit>
```
- GroundStationID: measurement point (default: the first ground station).
- TransmitLossesdB: losses between amplifier and antenna (default 0).
- ReferenceBandwidth: ITU reference bandwidth in Hz (default 4000, i.e. 4 kHz
  as applicable below 15 GHz; use 1e6 above).
- PfdLimit: the low-elevation limit in dB(W/m2) in the reference bandwidth
  (default -150, the Article 21 value for the 8025-8400 MHz EESS band).

<img src="/docs/com_pfd.png" alt="com_pfd"/>


### nav_dilution_of_precision
Plots navigation dilution of precision for a user grid. 

The following parameters are needed:
```
<Analysis>
    <Type>nav_dilution_of_precision</Type>
    <Direction>Ver</Direction>
    <Statistic>Max</Statistic>
</Analysis>
```
Parameters are:
- Direction: Hor/Ver/Pos direction of interest.
- Statistic: Min/Mean/Max/Median/Std statistic of interest

<img src="/docs/nav_dilution_of_precision.png" alt="nav_dilution_of_precision"/>

With <Plot3D>True</Plot3D> the DOP statistic is draped over the 3D globe (the example
below uses <ShowSatellite>False</ShowSatellite>; see the 3D plots section for all
parameters):
<img src="/docs/nav_dilution_of_precision_3d.png" alt="nav_dilution_of_precision_3d"/>


### nav_accuracy
Plots navigation accuracy for a user grid based on on sqrt(uere)*DOP computations. 

The following parameters are needed:
```
<Analysis>
    <Type>nav_dilution_of_precision</Type>
    <Direction>Ver</Direction>
    <Statistic>Max</Statistic>
</Analysis>
```
Parameters are:
- Direction: Hor/Ver/Pos direction of interest.
- Statistic: Min/Mean/Max/Median/Std statistic of interest

Additionally the constellation needs to be supplied with an elevation dependent list of uere values:
```
    <UERE>1.72,1.72,1.17,1.02,0.92,0.92,0.85,0.85,0.81,0.81,0.80,0.80,0.79,0.79,0.79,0.79,0.79,0.79</UERE>
```

<img src="/docs/nav_accuracy.png" alt="nav_accuracy"/>

With <Plot3D>True</Plot3D> the accuracy statistic is draped over the 3D globe (see
the 3D plots section for all parameters):
<img src="/docs/nav_accuracy_3d.png" alt="nav_accuracy_3d"/>


### sat_battery_depth_discharge
Models the satellite power subsystem over the simulation: it determines for each epoch whether the
satellite is in eclipse (geometric Earth-shadow check against the Sun direction), computes the solar
power generated, subtracts the bus and payload power draw, and integrates the battery state-of-charge.
The result plots the battery Depth-of-Discharge (DoD) in % together with the generated and drawn power.
The analysis uses the first satellite defined in the space segment.

The following parameters are needed:
```
<Analysis>
    <Type>sat_battery_depth_discharge</Type>
    <BatteryCapacityWh>500</BatteryCapacityWh>
    <InitialSoC>1.0</InitialSoC>
    <SolarPanelArea>5.0</SolarPanelArea>
    <PanelEfficiency>0.3</PanelEfficiency>
    <BasePowerDrawW>200</BasePowerDrawW>
    <InstrumentPowerDrawW>400</InstrumentPowerDrawW>
</Analysis>
```
Parameters are:
- BatteryCapacityWh: Battery capacity in Wh.
- InitialSoC: Initial state-of-charge as a fraction (1.0 = full).
- SolarPanelArea: Solar panel area in m^2.
- PanelEfficiency: Solar panel efficiency as a fraction (a solar constant of 1361 W/m^2 is assumed).
- BasePowerDrawW: Continuous bus power draw in W.
- InstrumentPowerDrawW: Additional payload power draw in W when the instrument is active.

Optional in the analysis part is:
```
    <PayloadLatitudeLimit>60</PayloadLatitudeLimit>
```
- PayloadLatitudeLimit: The payload only draws power when the satellite is below this absolute
  latitude in degrees (default 90, i.e. always on).

<img src="/docs/sat_battery_depth_discharge.png" alt="sat_battery_depth_discharge"/>


### sat_eclipse_duration
Detects eclipse entry/exit for the first satellite using the same geometric Earth-shadow check and
plots the duration in minutes of each eclipse over the simulation (day-of-year on the x-axis).
A small time step is recommended so eclipse transitions are captured accurately.

The following parameters are needed:
```
<Analysis>
    <Type>sat_eclipse_duration</Type>
</Analysis>
```

<img src="/docs/sat_eclipse_duration.png" alt="sat_eclipse_duration"/>


### sat_data_storage
Models the on-board Solid State Recorder (SSR). At each epoch data is recorded when the satellite is
below the payload latitude limit, and downlinked when a ground station is in view (using the
station-to-space links, so IncludeStation2SpaceLinks must be True). The SSR fill level is integrated
and clipped between zero and the capacity, and plotted over time with downlink-active periods shaded.
The analysis uses the first satellite defined in the space segment.

The following parameters are needed:
```
<Analysis>
    <Type>sat_data_storage</Type>
    <SSRCapacityGbits>2000</SSRCapacityGbits>
    <InitialFillGbits>0</InitialFillGbits>
    <InstrumentRateMbps>500</InstrumentRateMbps>
    <DownlinkRateMbps>1200</DownlinkRateMbps>
</Analysis>
```
Parameters are:
- SSRCapacityGbits: Recorder capacity in Gbits.
- InitialFillGbits: Initial amount of data stored in Gbits.
- InstrumentRateMbps: Instrument data generation rate in Mbps (while recording).
- DownlinkRateMbps: Downlink rate in Mbps (while a ground station is in view).

Optional in the analysis part is:
```
    <PayloadLatitudeLimit>60</PayloadLatitudeLimit>
```
- PayloadLatitudeLimit: The instrument only records when the satellite is below this absolute
  latitude in degrees (default 90, i.e. always recording).

<img src="/docs/sat_data_storage.png" alt="sat_data_storage"/>


### sat_data_latency
Extends the SSR model of sat_data_storage with a first-in-first-out data queue to compute the latency
between data acquisition and reception on the ground (orbit time until downlink plus a fixed ground
processing delay). It plots the latency time series and a histogram, and reports mean, 95th
percentile and the percentage of data received within 2 hours.
The analysis uses the first satellite defined in the space segment and requires station-to-space
links to be enabled.

The following parameters are needed:
```
<Analysis>
    <Type>sat_data_latency</Type>
    <SSRCapacityGbits>2000</SSRCapacityGbits>
    <InitialFillGbits>0</InitialFillGbits>
    <InstrumentRateMbps>500</InstrumentRateMbps>
    <DownlinkRateMbps>1200</DownlinkRateMbps>
</Analysis>
```
Parameters are the same as for sat_data_storage.

Optional in the analysis part are:
```
    <PayloadLatitudeLimit>60</PayloadLatitudeLimit>
    <GroundProcessingMin>15</GroundProcessingMin>
```
- PayloadLatitudeLimit: The instrument only records when the satellite is below this absolute
  latitude in degrees (default 90).
- GroundProcessingMin: Fixed ground processing delay in minutes added to each latency value (default 0).

<img src="/docs/sat_data_latency_timeseries.png" alt="sat_data_latency_timeseries"/>
<img src="/docs/sat_data_latency_histogram.png" alt="sat_data_latency_histogram"/>


### orb_kepler_elements
Plots the evolution of all osculating Kepler elements of the satellite(s) over the
simulation time, computed each epoch from the ECI state vector: semi-major axis,
eccentricity, inclination, RAAN, argument of perigee and mean anomaly (one plot
per element, orb_kepler_elements_semi_major_axis.png etc.). Run with the HPOP
propagator this shows the perturbation effects — e.g. the
drag decay of the semi-major axis (the log reports the secular change, averaging
out the J2 short-period oscillation) and the J2 RAAN drift. It also works with the
other propagators: constant elements for Keplerian, mean-element variations for
SGP4. For near-circular orbits the argument of perigee and mean anomaly are noisy
by nature (the perigee direction is poorly defined at low eccentricity).

The following parameters are needed:
```
<Analysis>
    <Type>orb_kepler_elements</Type>
</Analysis>
```
Optional are, to select one constellation or one satellite:
```
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
```
<img src="/docs/orb_kepler_elements_semi_major_axis.png" alt="orb_kepler_elements_semi_major_axis"/>
<img src="/docs/orb_kepler_elements_eccentricity.png" alt="orb_kepler_elements_eccentricity"/>
<img src="/docs/orb_kepler_elements_inclination.png" alt="orb_kepler_elements_inclination"/>
<img src="/docs/orb_kepler_elements_raan.png" alt="orb_kepler_elements_raan"/>
<img src="/docs/orb_kepler_elements_arg_perigee.png" alt="orb_kepler_elements_arg_perigee"/>
<img src="/docs/orb_kepler_elements_mean_anomaly.png" alt="orb_kepler_elements_mean_anomaly"/>


### orb_air_density
Plots the atmospheric density at the satellite altitude over the simulation time
(logarithmic scale), together with the altitude itself on a second axis. The
density is sampled every epoch from the atmosphere model configured as HPOP
DragModel (NRLMSISE00 or DTM2000 driven by the CSSI space weather data, or the
static HarrisPriester model), at the actual satellite position — so the day/night
bulge and the density increase while the orbit decays are visible. Requires
`<OrbitPropagator>HPOP</OrbitPropagator>`.

The following parameters are needed:
```
<Analysis>
    <Type>orb_air_density</Type>
</Analysis>
```
Optional are, to select one constellation or one satellite:
```
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
```
<img src="/docs/orb_air_density.png" alt="orb_air_density"/>


### orb_disturbance_forces
Plots the magnitude of every perturbation acceleration enabled in the HPOP force
model over the simulation time (logarithmic scale), evaluated each epoch on the
propagated state of the first selected satellite: geopotential harmonics,
atmospheric drag, solar radiation pressure, third-body Sun/Moon/planets, solid
and ocean tides and relativity, with the central gravity term as dashed
reference. This shows which perturbations dominate at the mission altitude and
how they vary along the orbit (e.g. drag peaks at perigee/day side, SRP dropping
to zero in eclipse). Requires `<OrbitPropagator>HPOP</OrbitPropagator>`.

The following parameters are needed:
```
<Analysis>
    <Type>orb_disturbance_forces</Type>
</Analysis>
```
Optional are, to select the satellite that is evaluated (the first match):
```
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
```
<img src="/docs/orb_disturbance_forces.png" alt="orb_disturbance_forces"/>


### orb_pole_wobble
Plots the wobble of the Earth rotation axis (the IERS polar motion xp/yp from
the EOP data of the Orekit archive) over the simulation time: the xp/yp time
series and the trace of the pole on the Earth surface. The Chandler and annual
wobble circle the mean pole by roughly 0.1–0.3 arcsec (several metres on the
ground) in about a year, so longer simulation windows (weeks to months) show a
larger arc of the circle. Requires `<OrbitPropagator>HPOP</OrbitPropagator>`
(for the Orekit EOP data).

The following parameters are needed:
```
<Analysis>
    <Type>orb_pole_wobble</Type>
</Analysis>
```
<img src="/docs/orb_pole_wobble_timeseries.png" alt="orb_pole_wobble_timeseries"/>
<img src="/docs/orb_pole_wobble_track.png" alt="orb_pole_wobble_track"/>


### orb_deltav_element
Estimates the station-keeping delta-v required to hold one orbit element within a
deadband around a target value over the simulation time. The element is the mean
element (the osculating value averaged over one orbital period); it drifts under
the modelled perturbations and whenever the controlled value leaves the deadband
an impulsive correction resets it to the target, costed with the standard
impulsive-maneuver formulas (tangential burn for Altitude/SemiMajorAxis/
Eccentricity, plane change at the node for Inclination/RAAN, apsidal rotation for
ArgOfPerigee). The plot shows the uncontrolled drift, the controlled element
bouncing inside the deadband with the maneuver epochs marked, and the cumulative
delta-v; the log reports the total and the extrapolated m/s per year. Works with
any propagator, but only the HPOP (and to a lesser degree SGP4) propagation
actually drifts — with Keplerian propagation or drift-free elements over short
windows the answer is simply zero, as it should be.

The following parameters are needed:
```
<Analysis>
    <Type>orb_deltav_element</Type>
    <TargetType>Altitude</TargetType>
    <DeadBand>1000</DeadBand>
</Analysis>
```
- TargetType: Altitude, SemiMajorAxis, Eccentricity, Inclination, RAAN or
  ArgOfPerigee.
- DeadBand: half width of the deadband around the target, in meters for
  Altitude/SemiMajorAxis, degrees for the angles, [-] for Eccentricity.

Optional are:
```
    <TargetValue>240000</TargetValue>
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
```
- TargetValue: the target for the mean element (same units as DeadBand). When
  omitted, the element value at the simulation start is held. Note the mean
  value can differ noticeably from the configured osculating element (e.g. the
  mean semi-major axis of a J2-perturbed LEO lies ~10 km below the initial
  osculating value).
- ConstellationID/SatelliteID: select the satellite (the first match is analysed).

<img src="/docs/orb_deltav_element.png" alt="orb_deltav_element"/>

### orb_beta_angle
Solar beta angle (the angle between the Sun direction and the orbit plane)
over the simulation time, together with the analytic eclipse fraction of a
circular orbit at that beta angle. Beta drives the eclipse pattern, the
thermal hot/cold cases and the power sizing, so this analysis ties the
sat_ platform analyses together; run it over months to see the seasonal cycle
(e.g. the classic ~60-day beat of an ISS-type orbit, or the near-constant
beta of a sun-synchronous orbit). The eclipse fraction drops to zero while
|beta| exceeds the shadow half-angle asin(R/r) — the eclipse-free season. The
log reports the beta range and the maximum eclipse minutes per orbit.

Note on propagators: Keplerian elements have no J2 nodal regression, so for
non-sun-synchronous orbits use SGP4 or HPOP (or an LTAN-defined SSO orbit) to
capture the full beta cycle.

All parameters are optional:
```
<Analysis>
    <Type>orb_beta_angle</Type>
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
</Analysis>
```
- ConstellationID/SatelliteID: select the satellite(s); omitted or 0 means all.

<img src="/docs/orb_beta_angle.png" alt="orb_beta_angle"/>

### orb_lifetime
Orbital lifetime under atmospheric drag. The mean semi-major axis at the end
of the simulation window is decayed semi-analytically (da/dt = -rho Cd A/m
sqrt(mu a), circular-orbit approximation, piecewise-exponential atmosphere
scaled by DensityScale for solar activity) until the re-entry altitude or the
MaxYears horizon. The analysis reports compliance with the 25-year
debris-mitigation lifetime rule, the delta-v of an immediate deorbit burn
(perigee lowered to the re-entry altitude) and, when the orbit is not
compliant, the circular disposal altitude with a 25-year lifetime and the
Hohmann delta-v to reach it. Works with any propagator — with HPOP the
projection starts from the actually decayed state at the end of the window.
Note the exponential atmosphere is a static mid-activity model: treat the
lifetime as an order-of-magnitude figure and bracket it with DensityScale
(~0.5 solar minimum, ~2 solar maximum).

The following parameters are needed (both fall back to the constellation's
Mass/FrontalArea when present there):
```
<Analysis>
    <Type>orb_lifetime</Type>
    <Mass>500</Mass>
    <DragArea>3.2</DragArea>
</Analysis>
```
- Mass: satellite mass in kg.
- DragArea: aerodynamic cross section in m^2.

Optional are:
```
    <DragCoefficient>2.2</DragCoefficient>
    <DensityScale>1.0</DensityScale>
    <MaxYears>100</MaxYears>
    <ReentryAltitude>120000</ReentryAltitude>
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
```
- DragCoefficient: drag coefficient Cd (default 2.2).
- DensityScale: multiplier on the atmosphere density for solar activity
  (default 1).
- MaxYears: integration horizon (default 100); an orbit that has not decayed
  by then is reported as "> MaxYears".
- ReentryAltitude: re-entry interface altitude in m (default 120 km).
- ConstellationID/SatelliteID: select the satellite (the first match is analysed).

<img src="/docs/orb_lifetime.png" alt="orb_lifetime"/>

### orb_environment
Space environment along the orbit, as a SPENVIS-style summary in six plots
(orb_environment_trapped_flux.png, _flux_timeseries, _drift_shell, _dose_depth,
_atomic_oxygen, _micrometeoroids):
- Trapped radiation: the geomagnetic field is modelled as an eccentric tilted dipole
  (the offset dipole centre is what creates the South Atlantic Anomaly), giving the
  McIlwain L-shell and field strength B along the orbit. AE8/AP8-style parametrisations
  of the trapped proton (>10 MeV) and electron (>1 MeV) omnidirectional fluxes — a
  Gaussian belt profile in L, an off-equator attenuation in B/B_eq and a drift-loss
  gate that confines the inner belts to the SAA at LEO altitudes — produce the SAA and
  polar horn crossings on the map, the flux time series and the annual fluences.
- Total ionizing dose versus aluminium shielding thickness (dose-depth curve) from the
  annual fluences with fluence-to-dose conversion factors (electrons dominating behind
  thin shielding, protons behind thick).
- Atomic oxygen: number density at the orbit altitude (MSIS-class mean-activity
  profile), ram fluence and kapton-equivalent erosion depth over the mission.
- Micrometeoroids: Gruen (1985) cumulative flux versus particle mass with Earth
  shielding and gravitational focusing, and the expected impact count on the
  spacecraft area over the mission.

All models are first-order engineering estimates that place the SAA/horn geometry
correctly and give order-of-magnitude fluxes and doses — use SPENVIS or IRENE
(AE9/AP9) for design and qualification values.

The following parameters are needed:
```
<Analysis>
    <Type>orb_environment</Type>
</Analysis>
```

Optional in the analysis part are:
```
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
    <MissionYears>5</MissionYears>
    <SurfaceArea>20</SurfaceArea>
```
- ConstellationID/SatelliteID: select the satellite (default: the first satellite).
- MissionYears: accumulation period for fluences, erosion and impact counts (default 5).
- SurfaceArea: area exposed to micrometeoroids in m^2 (default 10).

<img src="/docs/orb_environment_trapped_flux.png" alt="orb_environment_trapped_flux"/>
<img src="/docs/orb_environment_flux_timeseries.png" alt="orb_environment_flux_timeseries"/>
<img src="/docs/orb_environment_drift_shell.png" alt="orb_environment_drift_shell"/>
<img src="/docs/orb_environment_dose_depth.png" alt="orb_environment_dose_depth"/>
<img src="/docs/orb_environment_atomic_oxygen.png" alt="orb_environment_atomic_oxygen"/>
<img src="/docs/orb_environment_micrometeoroids.png" alt="orb_environment_micrometeoroids"/>

### sat_thermal
Single-node spacecraft thermal balance over the orbit. Per time step the heat
inputs — direct solar flux (zero in the Earth shadow, same eclipse model as the
power analyses), Earth albedo (scaled with the sun elevation over the subsatellite
point and the Earth view factor), Earth infrared and the internal electrical
dissipation — are balanced against the radiated heat (epsilon sigma A T^4) and
integrated to a temperature history. The plot shows the temperature (left axis)
and the individual heat flows (right axis); the log reports the temperature range
and the classic hot-case/cold-case equilibrium temperatures. Works with any
propagator.

The following parameters are needed:
```
<Analysis>
    <Type>sat_thermal</Type>
    <SurfaceArea>6.0</SurfaceArea>
    <HeatCapacity>50000</HeatCapacity>
</Analysis>
```
- SurfaceArea: total radiating surface area in m^2.
- HeatCapacity: thermal capacitance in J/K (roughly satellite mass x 900 J/kg/K
  for an aluminium-dominated bus).

Optional are:
```
    <CrossSectionSun>1.5</CrossSectionSun>
    <CrossSectionEarth>1.5</CrossSectionEarth>
    <Absorptivity>0.3</Absorptivity>
    <Emissivity>0.8</Emissivity>
    <InternalPowerW>300</InternalPowerW>
    <InitialTemperature>20</InitialTemperature>
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
```
- CrossSectionSun/CrossSectionEarth: projected areas towards the Sun and the
  Earth in m^2; default SurfaceArea/4 (a sphere presents a quarter of its
  surface to any direction).
- Absorptivity/Emissivity: solar absorptivity alpha (default 0.3) and infrared
  emissivity epsilon (default 0.8) of the outer surface.
- InternalPowerW: dissipated electrical power in W (default 0).
- InitialTemperature: start temperature in degC; when omitted the node starts in
  equilibrium with the first-epoch fluxes.
- ConstellationID/SatelliteID: select the satellite (the first match is analysed).

<img src="/docs/sat_thermal.png" alt="sat_thermal"/>

### sat_aocs
AOCS disturbance torques over the orbit with the standard worst-case models
(SMAD): gravity gradient (at a configurable worst-case attitude deviation),
aerodynamic torque (dynamic pressure on the drag area with the centre-of-pressure
lever arm), solar radiation pressure torque (zero in eclipse) and the magnetic
torque of the residual dipole in a tilted-dipole geomagnetic field. The
atmospheric density is sampled from the HPOP DragModel (NRLMSISE00 etc.) when
that propagator is active and from a built-in piecewise-exponential atmosphere
otherwise. The plot shows the torque components over time (log scale, left axis)
and the momentum buildup (the integral of the total torque — a conservative
reaction wheel sizing figure, right axis); the log reports the worst-case torque
per source and the momentum accumulation per orbit and per day.

The following parameters are needed:
```
<Analysis>
    <Type>sat_aocs</Type>
    <InertiaXX>100</InertiaXX>
    <InertiaYY>120</InertiaYY>
    <InertiaZZ>80</InertiaZZ>
</Analysis>
```
- InertiaXX/InertiaYY/InertiaZZ: principal moments of inertia in kg m^2 (the
  gravity gradient torque uses the largest difference).

Optional are:
```
    <MaxPointingOffset>1.0</MaxPointingOffset>
    <ResidualDipole>1.0</ResidualDipole>
    <DragArea>2.5</DragArea>
    <DragCoefficient>2.2</DragCoefficient>
    <SrpArea>2.5</SrpArea>
    <Reflectivity>0.6</Reflectivity>
    <CopOffset>0.2</CopOffset>
    <WheelMomentum>15</WheelMomentum>
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
```
- MaxPointingOffset: worst-case attitude deviation from the local vertical in
  degrees for the gravity gradient torque (default 1).
- ResidualDipole: residual magnetic dipole in A m^2 (default 1).
- DragArea/DragCoefficient: aerodynamic area in m^2 (default the constellation
  FrontalArea, else 1) and drag coefficient (default 2.2).
- SrpArea/Reflectivity: solar radiation pressure area in m^2 (default DragArea)
  and reflectance factor q (default 0.6).
- CopOffset: centre-of-pressure to centre-of-mass offset in m, the lever arm of
  the aerodynamic and SRP torques (default 0.1).
- WheelMomentum: optional reaction wheel momentum capacity in N m s, drawn as a
  line on the momentum plot.
- ConstellationID/SatelliteID: select the satellite (the first match is analysed).

<img src="/docs/sat_aocs.png" alt="sat_aocs"/>

### sat_drag_coefficient
Estimates the aerodynamic drag coefficient from the satellite geometry with the
Sentman free-molecular panel method — the correct flow regime at orbital altitudes
(the molecular mean free path is kilometres, so no CFD is involved). Every facet of
the satellite mesh contributes drag through the analytic diffuse-reflection
gas-surface interaction, with ray-cast shadowing between panels. The mesh is the
`<SatelliteModelFile>` STL taken at true dimensions in metres (needs pyvista, the
same dependency as the 3D plots), or the built-in bus + solar panel model.

Outputs:
- the drag area CdA over a body-frame attitude sweep of the flow direction
  (map plot with the ram (+x) attitude marked, plus the tumbling average),
- CdA along the orbit for a nadir-fixed spacecraft (+x flight direction; the
  molecular speed ratio follows the altitude through an NRLMSISE-class mean
  composition and temperature profile),
- a recommended `<DragArea>`/`<DragCd>` pair for the HPOP force model and the
  orb_/sat_ analyses, logged and dumped to CSV.

Note that satellites with large surfaces flying edge-on (solar panels) genuinely
show Cd well above the classic 2.2 — the grazing thermal flux on those surfaces is
real free-molecular drag. The energy accommodation coefficient dominates the model
uncertainty (order 10-20% on Cd); surfaces covered by adsorbed atomic oxygen below
~500 km are nearly fully diffuse (accommodation 0.9-1.0).

The following parameters are needed:
```
<Analysis>
    <Type>sat_drag_coefficient</Type>
</Analysis>
```

Optional in the analysis part are:
```
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
    <SatelliteModelFile>../input/satellite.stl</SatelliteModelFile>
    <SatelliteModelScale>1.0</SatelliteModelScale>
    <AccommodationCoefficient>0.93</AccommodationCoefficient>
    <WallTemperature>300</WallTemperature>
    <ExosphericTemperature>1000</ExosphericTemperature>
    <AttitudeStep>15</AttitudeStep>
    <Shadowing>True</Shadowing>
    <MaxFacets>5000</MaxFacets>
```
- ConstellationID/SatelliteID: select the satellite (default: the first satellite).
- SatelliteModelFile: STL mesh in metres (default: the built-in bus + panel model);
  SatelliteModelScale rescales STL units to metres.
- AccommodationCoefficient: energy accommodation of the gas-surface interaction
  (default 0.93).
- WallTemperature: spacecraft surface temperature in K for the re-emission term
  (default 300).
- ExosphericTemperature: atmosphere temperature scale in K (default 1000).
- AttitudeStep: attitude sweep resolution in degrees (default 15).
- Shadowing: ray-cast panel-on-panel shadowing (default True; switch off for speed
  on convex shapes).
- MaxFacets: larger meshes are decimated to this facet count (default 5000).

<img src="/docs/sat_drag_coefficient.png" alt="sat_drag_coefficient"/>
<img src="/docs/sat_drag_coefficient_orbit.png" alt="sat_drag_coefficient_orbit"/>





