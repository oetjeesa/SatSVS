# Satellite Service Volume Simulator
## Open Source satellite, ground station and user tool
### M. Tossaint - 2026 - v3

<img src="/docs/schema.png" alt="schema"/>

## Installation & first run
Download from github, and install the following libraries:
- Numpy, pandas and numba
- Astropy and sgp4
- Cartopy and xarray
- Geopandas and shapely
- Itur
- For the HPOP orbit propagator only: orekit_jpype and jdk4py (bundled JVM), plus the
  Orekit physical data archive saved as __input/orekit-data.zip__
  (download from https://gitlab.orekit.org/orekit/orekit-data)
- For the 3D plots only: pyvista, plus an equirectangular Earth texture saved as
  __input/earth_texture.jpg__ (e.g. the public domain NASA Blue Marble image
  https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57752/land_shallow_topo_2048.jpg)

To run, edit the config.xml file and run: python main.py

## Introduction
Framework takes care of geometry computations, satellite propagation, ground station and user rotation in ECI/ECF.
It will also automatically compute links between stations and satellites, users and satellites, and between satellites.

### Main structure of the tool

<img src="/docs/satsvs_architecture.png" alt="schema"/>

### Configuration of the tool

Configuration of the tool can be done in the config.xml file where satellites, ground stations, users, simulation
parameters and analysis are defined. Analysis can be added as wished, the baseline of analysis available are below 
(and explained further below):

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

### Communication
- __com_gr2sp_budget__: For station-satellite received power, losses and C/N0
- __com_sp2sp_budget__: For satellite-satellite received power, losses and C/N0
- __com_doppler__: For satellite-station elevation and doppler

### Navigation
- __nav_dillution_of_precision__: DOP values for user(s) (also spacecraft user)
- __nav_accuracy__: Navigation accuracy (UERE*DOP) values for user(s) (also spacecraft user)

### Satellite power
- __pow_battery_depth_discharge__: Battery state-of-charge / depth-of-discharge and power generation vs. draw over orbit
- __pow_eclipse_duration__: Eclipse duration per orbit over the simulation

### Data handling
- __dat_storage__: Solid State Recorder (SSR) fill level over orbit (recording vs. downlink)
- __dat_latency__: Data latency statistics from acquisition to ground reception

### Orbit
- __orb_semi_major_axis__: Osculating semi-major axis over time, e.g. the orbital decay
  under atmospheric drag with the HPOP propagator

_To be implemented at a later stage:_

### Satellite
- __sat_thermal__: Thermal conditions over orbit
- __sat_attitude_control__: Satellite attitude control over orbit

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

    <Analysis>
          <Type>9</Type>
          <ConstellationID>1</ConstellationID>
    </Analysis>
</SimulationManager>
```
The following explanations apply for the parameters:
- The Start/Stop time parameters are in UTC time and TimeStep in seconds. 
- The IncludeStation2SpaceLinks, etc. parameters determine whether links between different objects: sat, station and user are computed. Normally leave these to True so that all analysis works. Time could
be saved by disabling some. 
- The OrbitsFromPreviousRun flag (True/False) reuses the satellite ECI orbits cached in 'output/orbits_internal.txt' from a previous run instead of re-propagating, to save time when only the analysis changes.
- The OrbitPropagator determines which propagator to take: 'Keplerian', 'SGP4' or 'HPOP'.

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
Ground stations are always drawn as magenta markers. Near-Earth scenes use a
perspective camera above the (last) satellite; scenes whose orbits reach far
above the Earth (MEO/GEO, e.g. GNSS constellations) are rendered with a fitted
parallel projection instead, so the Earth and the orbits appear at their true
relative size.

The specific parameters for the existing analysis are given here below:

### cov_ground_track
Plots the ground track of one or more satellites over simulation time. The following parameters are needed, to plot the ground track of satellites in a constellation:
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
The user segment is used to define the grid of analysis and defines the granularity of the result.
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
what kind of statistic is displayed per user location.
- SaveOutput: NetCDF or Numpy This flag will enable saving user swath coverage for every timestep.
- Plot3D flag: additionally renders the swath as a semi-transparent ribbon on a
  textured 3D Earth, saved as output/obs_swath_conical_3d.png (see the 3D plots
  section for all Plot3D parameters).

<img src="/docs/obs_swath_conical.png" alt="obs_swath_conical"/>
<img src="/docs/obs_swath_conical_revisit.png" alt="obs_swath_conical_revisit"/>
<img src="/docs/obs_swath_conical_3d.png" alt="obs_swath_conical_3d"/>

### obs_swath_push_broom
Plots the swath coverage for a push broom scanner on one or more satellites defined in the space segment. 
The user segment is used to define the grid of analysis and defines the granularity of the result.
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
what kind of statistic is displayed per user location.
- SaveOutput: NetCDF or Numpy This flag will enable saving user swath coverage for every timestep.
- Plot3D flag: additionally renders the swath as a semi-transparent ribbon on a
  textured 3D Earth, saved as output/obs_swath_push_broom_3d.png (see the 3D plots
  section for all Plot3D parameters).

<img src="/docs/obs_swath_push_broom.png" alt="cov_satellite_push_broom"/>
<img src="/docs/obs_swath_push_broom_revisit.png" alt="cov_satellite_push_broom_revisit"/>
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
```

Some parameters can be entered to get the required CN0:
- ModulationType: 'BPSK' or 'QPSK'
- BitErrorRate: bit error rate required
- DataRateBitPerSec: datarate required

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


### pow_battery_depth_discharge
Models the satellite power subsystem over the simulation: it determines for each epoch whether the
satellite is in eclipse (geometric Earth-shadow check against the Sun direction), computes the solar
power generated, subtracts the bus and payload power draw, and integrates the battery state-of-charge.
The result plots the battery Depth-of-Discharge (DoD) in % together with the generated and drawn power.
The analysis uses the first satellite defined in the space segment.

The following parameters are needed:
```
<Analysis>
    <Type>pow_battery_depth_discharge</Type>
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

<img src="/docs/pow_depth_discharge.png" alt="pow_battery_depth_discharge"/>


### pow_eclipse_duration
Detects eclipse entry/exit for the first satellite using the same geometric Earth-shadow check and
plots the duration in minutes of each eclipse over the simulation (day-of-year on the x-axis).
A small time step is recommended so eclipse transitions are captured accurately.

The following parameters are needed:
```
<Analysis>
    <Type>pow_eclipse_duration</Type>
</Analysis>
```

<img src="/docs/pow_eclipse_duration.png" alt="pow_eclipse_duration"/>


### dat_storage
Models the on-board Solid State Recorder (SSR). At each epoch data is recorded when the satellite is
below the payload latitude limit, and downlinked when a ground station is in view (using the
station-to-space links, so IncludeStation2SpaceLinks must be True). The SSR fill level is integrated
and clipped between zero and the capacity, and plotted over time with downlink-active periods shaded.
The analysis uses the first satellite defined in the space segment.

The following parameters are needed:
```
<Analysis>
    <Type>dat_storage</Type>
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

<img src="/docs/dat_storage.png" alt="dat_storage"/>


### dat_latency
Extends the SSR model of dat_storage with a first-in-first-out data queue to compute the latency
between data acquisition and reception on the ground (orbit time until downlink plus a fixed ground
processing delay). It plots the latency time series and a histogram, and reports mean, 95th
percentile and the percentage of data received within 2 hours.
The analysis uses the first satellite defined in the space segment and requires station-to-space
links to be enabled.

The following parameters are needed:
```
<Analysis>
    <Type>dat_latency</Type>
    <SSRCapacityGbits>2000</SSRCapacityGbits>
    <InitialFillGbits>0</InitialFillGbits>
    <InstrumentRateMbps>500</InstrumentRateMbps>
    <DownlinkRateMbps>1200</DownlinkRateMbps>
</Analysis>
```
Parameters are the same as for dat_storage.

Optional in the analysis part are:
```
    <PayloadLatitudeLimit>60</PayloadLatitudeLimit>
    <GroundProcessingMin>15</GroundProcessingMin>
```
- PayloadLatitudeLimit: The instrument only records when the satellite is below this absolute
  latitude in degrees (default 90).
- GroundProcessingMin: Fixed ground processing delay in minutes added to each latency value (default 0).

<img src="/docs/dat_latency_stats.png" alt="dat_latency"/>


### orb_semi_major_axis
Plots the osculating semi-major axis of the satellite(s) over the simulation time,
computed from the ECI state vector with the vis-viva equation. Run with the HPOP
propagator and drag enabled this shows the orbital decay (the legend reports the
secular change, i.e. the difference of the mean semi-major axis at the start and
the end of the simulation, averaging out the J2 short-period oscillation). It also
works with the other propagators: constant semi-major axis for Keplerian,
mean-element variations for SGP4.

The following parameters are needed:
```
<Analysis>
    <Type>orb_semi_major_axis</Type>
</Analysis>
```
Optional are, to select one constellation or one satellite:
```
    <ConstellationID>1</ConstellationID>
    <SatelliteID>1</SatelliteID>
```
<img src="/docs/orb_semi_major_axis.png" alt="semi_major_axis"/>





