"""
Generates the Config.xml (and copies the input files, e.g. TLE files) for every
analysis test folder under tests/. Run once:  py tests/make_configs.py
Each scenario is chosen to be representative for the analysis under test
(GPS constellation for navigation, LEO SSO satellites for EO/power/data, etc.)
while keeping the runtime of a single test in the minutes range.
"""
import os
import re
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS = os.path.join(ROOT, 'tests')
TLE_DIR = os.path.join(ROOT, 'input', 'example_tle_files')

UERE_GPS = ('1.72,1.72,1.17,1.02,0.92,0.92,0.85,0.85,0.81,0.81,'
            '0.80,0.80,0.79,0.79,0.79,0.79,0.79,0.79')


def gps_constellation(uere=False):
    """24-satellite GPS Walker constellation from the example config blocks."""
    with open(os.path.join(ROOT, 'input', 'example_conf_blocks', 'ConstellationGPS.xml')) as f:
        txt = f.read()
    m = re.search(r'<Constellation>.*</Constellation>', txt, re.DOTALL)
    block = m.group(0)
    if uere:
        block = block.replace('<ReceiverConstellation>1111</ReceiverConstellation>',
                              '<ReceiverConstellation>1111</ReceiverConstellation>\n'
                              f'      <UERE>{UERE_GPS}</UERE>')
    return block


def tle_constellation(name, tle_file, test_name, extra=''):
    """Constellation defined by a TLE file placed inside the test folder."""
    return f"""<Constellation>
      <ConstellationID>1</ConstellationID>
      <NumOfSatellites>1</NumOfSatellites>
      <NumOfPlanes>1</NumOfPlanes>
      <ConstellationName>{name}</ConstellationName>
      <ReceiverConstellation>1</ReceiverConstellation>
      <TLEFileName>../tests/{test_name}/{tle_file}</TLEFileName>
{extra}    </Constellation>"""


def sso_constellation(num_sat=1, extra='', ltan=False, delta_ma_deg=2.0, raan=50.0,
                      sma=7078137, incl=98.19):
    """LEO sun-synchronous satellite(s), Keplerian elements, ~700 km altitude."""
    sats = ''
    for i in range(num_sat):
        if ltan:
            orbit = """<EpochMJD>61072.5</EpochMJD>
        <Altitude>700000</Altitude>
        <Eccentricity>0.001</Eccentricity>
        <LTAN>22.25</LTAN>"""
        else:
            orbit = f"""<EpochMJD>61072.0</EpochMJD>
        <SemiMajorAxis>{sma}</SemiMajorAxis>
        <Eccentricity>0.001</Eccentricity>
        <Inclination>{incl}</Inclination>
        <RAAN>{raan}</RAAN>"""
        sats += f"""      <Satellite>
        <SatelliteID>{i + 1}</SatelliteID>
        <Plane>1</Plane>
        {orbit}
        <ArgOfPerigee>0.0</ArgOfPerigee>
        <MeanAnomaly>{i * delta_ma_deg}</MeanAnomaly>
      </Satellite>
"""
    return f"""<Constellation>
      <ConstellationID>1</ConstellationID>
      <NumOfSatellites>{num_sat}</NumOfSatellites>
      <NumOfPlanes>1</NumOfPlanes>
      <ConstellationName>LEO-SSO</ConstellationName>
      <ReceiverConstellation>1</ReceiverConstellation>
{extra}{sats}    </Constellation>"""


STATIONS = {
    'Svalbard': (77.875, 20.9752),
    'Kiruna': (67.83, 20.42),
    'Inuvik': (68.24, -133.39),
    'Kourou': (5.25, -52.80),
}


def ground_segment(names, mask=5):
    gs = ''
    for i, n in enumerate(names):
        lat, lon = STATIONS[n]
        gs += f"""    <GroundStation>
      <Type>Downlink</Type>
      <ConstellationID>1</ConstellationID>
      <GroundStationID>{i + 1}</GroundStationID>
      <GroundStationName>{n}</GroundStationName>
      <Latitude>{lat}</Latitude>
      <Longitude>{lon}</Longitude>
      <Height>0</Height>
      <ReceiverConstellation>1</ReceiverConstellation>
      <ElevationMask>{mask}</ElevationMask>
    </GroundStation>
"""
    return f"""<GroundSegment>
    <Network>
      <NumStation>{len(names)}</NumStation>
      <NetworkName>Net1</NetworkName>
{gs}    </Network>
  </GroundSegment>"""


def static_users(locations, mask=5):
    us = ''
    for lat, lon in locations:
        us += f"""    <User>
      <Type>Static</Type>
      <Latitude>{lat}</Latitude>
      <Longitude>{lon}</Longitude>
      <Height>0</Height>
      <ReceiverConstellation>1</ReceiverConstellation>
      <ElevationMask>{mask}</ElevationMask>
    </User>
"""
    return f'<UserSegment>\n{us}  </UserSegment>'


def grid_users(step, mask=5):
    return f"""<UserSegment>
    <User>
      <Type>Grid</Type>
      <LatMin>-90</LatMin>
      <LatMax>90</LatMax>
      <LonMin>-180</LonMin>
      <LonMax>180</LonMax>
      <LatStep>{step}</LatStep>
      <LonStep>{step}</LonStep>
      <Height>0</Height>
      <ReceiverConstellation>1</ReceiverConstellation>
      <ElevationMask>{mask}</ElevationMask>
    </User>
  </UserSegment>"""


def polygon_users(name, points, step, mask=5):
    """Grid users clipped to a polygon: the AOI of the obs_aoi_revisit test.
    points is a string of (lon, lat) tuples as documented in readme.md."""
    return f"""<UserSegment>
    <User>
      <Type>Polygon</Type>
      <Name>{name}</Name>
      <PolygonList>{points}</PolygonList>
      <LatStep>{step}</LatStep>
      <LonStep>{step}</LonStep>
      <Height>0</Height>
      <ReceiverConstellation>1</ReceiverConstellation>
      <ElevationMask>{mask}</ElevationMask>
    </User>
  </UserSegment>"""


def config(space, ground, users, start, stop, step, analysis,
           propagator='Keplerian', gr2sp=True, usr2sp=True, sp2sp=False, extra_sim=''):
    if isinstance(analysis, str):  # One or a list of <Analysis> block contents
        analysis = [analysis]
    analysis_blocks = '\n'.join(f'    <Analysis>\n{a}\n    </Analysis>' for a in analysis)
    return f"""<?xml version="1.0" encoding="utf-8"?>
<!--Simulation Scenario-->
<Scenario>
  <SpaceSegment>
    {space}
  </SpaceSegment>
  {ground}
  {users}
  <SimulationManager>
    <StartDate>{start}</StartDate>
    <StopDate>{stop}</StopDate>
    <TimeStep>{step}</TimeStep>
    <IncludeStation2SpaceLinks>{gr2sp}</IncludeStation2SpaceLinks>
    <IncludeUser2SpaceLinks>{usr2sp}</IncludeUser2SpaceLinks>
    <IncludeSpace2SpaceLinks>{sp2sp}</IncludeSpace2SpaceLinks>
    <OrbitsFromPreviousRun>False</OrbitsFromPreviousRun>
    <OrbitPropagator>{propagator}</OrbitPropagator>
{extra_sim}{analysis_blocks}
  </SimulationManager>
</Scenario>
"""


def write_test(name, cfg, tle=None):
    d = os.path.join(TESTS, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, 'Config.xml'), 'w') as f:
        f.write(cfg)
    if tle:
        shutil.copy(os.path.join(TLE_DIR, tle), d)
    print(f'wrote {name}')


DELFT = [(52.0, 4.36)]
FEB26 = ('2026-02-01 00:00:00', '2026-02-02 00:00:00')
APR20 = ('2020-04-07 00:00:00', '2020-04-08 00:00:00')  # Sentinel-1 TLE epoch
SEP18 = ('2018-09-18 00:00:00', '2018-09-19 00:00:00')  # MetOp-A TLE epoch

SWATH_PB = """      <ObsSwathStart>250000.0</ObsSwathStart>
      <ObsSwathStop>650000.0</ObsSwathStop>
"""
SWATH_CON = """      <ObsSwathStop>650000.0</ObsSwathStop>
"""

# ----------------------------------------------------------------- coverage
write_test('cov_ground_track', config(
    tle_constellation('TerraSAR-X', 'terrasarx.txt', 'cov_ground_track'),
    ground_segment(['Kiruna']), static_users(DELFT), *FEB26, 60,
    '      <Type>cov_ground_track</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>31698</SatelliteID>\n'
    '      <EarthImage>True</EarthImage>\n'
    '      <Coastlines>False</Coastlines>\n'
    '      <ShowStations>True</ShowStations>\n'
    '      <ShowUsers>True</ShowUsers>\n'
    '      <Plot3D>True</Plot3D>\n'
    '      <EarthClouds>True</EarthClouds>\n'
    '      <MP4>True</MP4>',
    propagator='SGP4'), tle='terrasarx.txt')

write_test('cov_satellite_pvt', config(
    gps_constellation(), ground_segment(['Kourou']), static_users(DELFT), *FEB26, 300,
    '      <Type>cov_satellite_pvt</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>'))

write_test('cov_satellite_visible', config(
    gps_constellation(), ground_segment(['Kourou']),
    static_users([(52.0, 4.36), (1.35, 103.82)]), *FEB26, 300,
    '      <Type>cov_satellite_visible</Type>'))

write_test('cov_satellite_visible_grid', config(
    gps_constellation(), ground_segment(['Kourou']), grid_users(10), *FEB26, 300,
    '      <Type>cov_satellite_visible_grid</Type>\n'
    '      <Statistic>Min</Statistic>\n'
    '      <Plot3D>True</Plot3D>\n'
    '      <ShowOrbit>False</ShowOrbit>\n'
    '      <SatelliteModelScale>500000</SatelliteModelScale>\n'
    '      <MP4>True</MP4>'))

write_test('cov_satellite_visible_id', config(
    gps_constellation(), ground_segment(['Kourou']), static_users(DELFT), *FEB26, 120,
    '      <Type>cov_satellite_visible_id</Type>\n'
    '      <ConstellationID>1</ConstellationID>'))

write_test('cov_satellite_contour', config(
    gps_constellation(), ground_segment(['Kourou']), static_users(DELFT), *FEB26, 300,
    '      <Type>cov_satellite_contour</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <ElevationMask>10</ElevationMask>\n'
    '      <Plot3D>True</Plot3D>\n'
    '      <SatelliteModelScale>500000</SatelliteModelScale>'))

write_test('cov_satellite_sky_angles', config(
    gps_constellation(), ground_segment(['Kourou']), static_users(DELFT), *FEB26, 120,
    '      <Type>cov_satellite_sky_angles</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>'))

write_test('cov_satellite_highest', config(
    gps_constellation(), ground_segment(['Kourou']), grid_users(10), *FEB26, 300,
    '      <Type>cov_satellite_highest</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <Statistic>Mean</Statistic>\n'
    '      <Plot3D>True</Plot3D>\n'
    '      <SatelliteModelScale>500000</SatelliteModelScale>'))

write_test('cov_depth_of_coverage', config(
    tle_constellation('TerraSAR-X', 'terrasarx.txt', 'cov_depth_of_coverage'),
    ground_segment(['Svalbard', 'Kiruna', 'Inuvik']), static_users(DELFT), *FEB26, 60,
    '      <Type>cov_depth_of_coverage</Type>\n'
    '      <Plot3D>True</Plot3D>',
    propagator='SGP4'), tle='terrasarx.txt')

write_test('cov_pass_time', config(
    gps_constellation(), ground_segment(['Kourou']), grid_users(10), *FEB26, 300,
    '      <Type>cov_pass_time</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <Statistic>Mean</Statistic>\n'
    '      <Plot3D>True</Plot3D>\n'
    '      <ShowSatellite>False</ShowSatellite>'))

# ------------------------------------------------------------ earth observation
write_test('obs_swath_conical', config(
    tle_constellation('MetOp-A', 'metop_a.txt', 'obs_swath_conical', extra=SWATH_CON),
    ground_segment(['Svalbard']), grid_users(2, mask=0), *SEP18, 60,
    '      <Type>obs_swath_conical</Type>\n'
    '      <Revisit>True</Revisit>\n'
    '      <Statistic>Mean</Statistic>\n'
    '      <SaveOutput>Numpy</SaveOutput>\n'
    '      <Plot3D>True</Plot3D>\n'
    '      <MP4>True</MP4>',
    propagator='SGP4', usr2sp=False), tle='metop_a.txt')

write_test('obs_swath_push_broom', config(
    tle_constellation('Sentinel-1', 'sentinel-1.txt', 'obs_swath_push_broom', extra=SWATH_PB),
    ground_segment(['Svalbard']), grid_users(2, mask=0), *APR20, 60,
    '      <Type>obs_swath_push_broom</Type>\n'
    '      <Revisit>True</Revisit>\n'
    '      <Statistic>Mean</Statistic>\n'
    '      <SaveOutput>NetCDF</SaveOutput>\n'
    '      <Plot3D>True</Plot3D>',
    propagator='SGP4', usr2sp=False), tle='sentinel-1.txt')

write_test('obs_sza_push_broom', config(
    tle_constellation('Sentinel-1', 'sentinel-1.txt', 'obs_sza_push_broom', extra=SWATH_PB),
    ground_segment(['Svalbard']), grid_users(5, mask=0),
    '2020-04-07 00:00:00', '2020-04-07 12:00:00', 60,
    '      <Type>obs_sza_push_broom</Type>\n'
    '      <Plot3D>True</Plot3D>',
    propagator='SGP4', usr2sp=False), tle='sentinel-1.txt')

write_test('obs_sza_subsat', config(
    tle_constellation('Sentinel-1', 'sentinel-1.txt', 'obs_sza_subsat'),
    ground_segment(['Svalbard']), static_users(DELFT),
    '2020-04-07 00:00:00', '2020-04-09 00:00:00', 60,
    '      <Type>obs_sza_subsat</Type>\n'
    '      <RangeLatitude>-80,80,10</RangeLatitude>\n'
    '      <SaveOutput>Numpy</SaveOutput>\n'
    '      <Plot3D>True</Plot3D>',
    propagator='SGP4'), tle='sentinel-1.txt')

IBERIA_FRANCE = ('(-10.0, 36.0),(3.5, 36.0),(8.0, 43.0),(8.0, 49.0),'
                 '(2.0, 51.0),(-2.0, 49.0),(-10.0, 44.0)')

write_test('obs_aoi_revisit', config(
    # 3 days so most AOI points collect several push-broom passes (revisit gaps)
    sso_constellation(extra=SWATH_PB), ground_segment(['Svalbard']),
    polygon_users('IberiaFrance', IBERIA_FRANCE, 1.0, mask=0),
    '2026-02-01 00:00:00', '2026-02-04 00:00:00', 60,
    '      <Type>obs_aoi_revisit</Type>\n'
    '      <Statistic>Max</Statistic>',
    usr2sp=False))

write_test('obs_target_imaging', config(
    # RAAN 140 gives a noon/midnight orbit in Feb 2026: the dayside passes
    # clear the MinSunElevation daylight constraint at the mid-latitude
    # targets (Longyearbyen stays in polar night - zero opportunities)
    sso_constellation(raan=140.0), ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-02-03 00:00:00', 30,
    '      <Type>obs_target_imaging</Type>\n'
    '      <ShowGroundTrack>True</ShowGroundTrack>\n'
    '      <MaxOffNadir>45</MaxOffNadir>\n'
    '      <MinSunElevation>10</MinSunElevation>\n'
    '      <Target>Rome, 41.9, 12.5</Target>\n'
    '      <Target>Delft, 52.0, 4.36</Target>\n'
    '      <Target>Longyearbyen, 78.2, 15.6</Target>\n'
    '      <Target>SaoPaulo, -23.5, -46.6</Target>\n'
    '      <Target>Sydney, -33.9, 151.2</Target>'))

# ------------------------------------------------------------- communication
write_test('com_gr2sp_budget', config(
    tle_constellation('TerraSAR-X', 'terrasarx.txt', 'com_gr2sp_budget'),
    ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 60,
    '      <Type>com_gr2sp_budget</Type>\n'
    '      <GroundStationID>1</GroundStationID>\n'
    '      <TransmitterObject>Satellite</TransmitterObject>\n'
    '      <CarrierFrequency>8.025e9</CarrierFrequency>\n'
    '      <TransmitPowerW>20</TransmitPowerW>\n'
    '      <TransmitLossesdB>2</TransmitLossesdB>\n'
    '      <TransmitGaindB>20</TransmitGaindB>\n'
    '      <TransmitAntennaPatternFile>../input/example_antenna_patterns/isoflux_x_band.cut</TransmitAntennaPatternFile>\n'
    '      <ReceiveGaindB>45</ReceiveGaindB>\n'
    '      <ReceiveAntennaPatternFile>../input/example_antenna_patterns/dish_x_band_3m.grd</ReceiveAntennaPatternFile>\n'
    '      <ReceiveLossesdB>3</ReceiveLossesdB>\n'
    '      <ReceiveTempK>200</ReceiveTempK>\n'
    '      <PExceedPerc>0.5</PExceedPerc>\n'
    '      <IncludeGas>True</IncludeGas>\n'
    '      <IncludeRain>True</IncludeRain>\n'
    '      <IncludeScintillation>False</IncludeScintillation>\n'
    '      <IncludeClouds>False</IncludeClouds>\n'
    '      <ModulationType>QPSK</ModulationType>\n'
    '      <BitErrorRate>1e-5</BitErrorRate>\n'
    '      <DataRateBitPerSec>100e6</DataRateBitPerSec>',
    propagator='SGP4'), tle='terrasarx.txt')

write_test('com_gr2sp_budget_interference', config(
    sso_constellation(num_sat=2, delta_ma_deg=1.0),
    ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 60,
    '      <Type>com_gr2sp_budget_interference</Type>\n'
    '      <GroundStationID>1</GroundStationID>\n'
    '      <TransmitterObject>Satellite</TransmitterObject>\n'
    '      <CarrierFrequency>26.25e9</CarrierFrequency>\n'
    '      <BandWidth>675e6</BandWidth>\n'
    '      <TransmitPowerW>50</TransmitPowerW>\n'
    '      <TransmitLossesdB>5</TransmitLossesdB>\n'
    '      <TransmitGaindB>33</TransmitGaindB>\n'
    '      <TransmitAntennaDiameter>0.25</TransmitAntennaDiameter>\n'
    '      <TransmitAntennaPatternFile>../input/example_antenna_patterns/dish_ka_band_25cm.cut</TransmitAntennaPatternFile>\n'
    '      <ReceiveGaindB>64</ReceiveGaindB>\n'
    '      <ReceiveAntennaDiameter>6.8</ReceiveAntennaDiameter>\n'
    '      <ReceiveAntennaPatternFile>../input/example_antenna_patterns/dish_ka_band_6_8m.cut</ReceiveAntennaPatternFile>\n'
    '      <ReceiveLossesdB>3</ReceiveLossesdB>\n'
    '      <ReceiveTempK>290</ReceiveTempK>\n'
    '      <PExceedPerc>0.5</PExceedPerc>\n'
    '      <IncludeGas>True</IncludeGas>\n'
    '      <IncludeRain>True</IncludeRain>\n'
    '      <IncludeScintillation>False</IncludeScintillation>\n'
    '      <IncludeClouds>False</IncludeClouds>'))

write_test('com_sp2sp_budget', config(
    gps_constellation(), ground_segment(['Kourou']), static_users(DELFT), *FEB26, 300,
    '      <Type>com_sp2sp_budget</Type>\n'
    '      <SatelliteID1>1</SatelliteID1>\n'
    '      <SatelliteID2>2</SatelliteID2>\n'
    '      <CarrierFrequency>23e9</CarrierFrequency>\n'
    '      <TransmitPowerW>10</TransmitPowerW>\n'
    '      <TransmitLossesdB>2</TransmitLossesdB>\n'
    '      <TransmitGaindB>30</TransmitGaindB>\n'
    '      <ReceiveGaindB>30</ReceiveGaindB>\n'
    '      <ReceiveLossesdB>2</ReceiveLossesdB>\n'
    '      <ReceiveTempK>290</ReceiveTempK>\n'
    '      <ModulationType>BPSK</ModulationType>\n'
    '      <BitErrorRate>1e-5</BitErrorRate>\n'
    '      <DataRateBitPerSec>1e6</DataRateBitPerSec>',
    sp2sp=True))

write_test('com_doppler', config(
    tle_constellation('TerraSAR-X', 'terrasarx.txt', 'com_doppler'),
    ground_segment(['Kiruna']), static_users(DELFT), *FEB26, 30,
    '      <Type>com_doppler</Type>\n'
    '      <StationID>1</StationID>\n'
    '      <CarrierFrequency>8.025e9</CarrierFrequency>',
    propagator='SGP4'), tle='terrasarx.txt')

write_test('com_contact_plan', config(
    # 2 co-planar SSO satellites 10 deg apart in mean anomaly: their passes
    # overlap at every station, exercising the conflict detection
    sso_constellation(num_sat=2, delta_ma_deg=10.0),
    ground_segment(['Svalbard', 'Inuvik']), static_users(DELFT), *FEB26, 30,
    '      <Type>com_contact_plan</Type>\n'
    '      <MinDuration>60</MinDuration>\n'
    '      <DownlinkRateMbps>1070</DownlinkRateMbps>'))

write_test('com_pfd', config(
    tle_constellation('TerraSAR-X', 'terrasarx.txt', 'com_pfd'),
    ground_segment(['Svalbard'], mask=0), static_users(DELFT), *FEB26, 30,
    '      <Type>com_pfd</Type>\n'
    '      <GroundStationID>1</GroundStationID>\n'
    '      <CarrierFrequency>8025e6</CarrierFrequency>\n'
    '      <TransmitPowerW>20</TransmitPowerW>\n'
    '      <TransmitLossesdB>1</TransmitLossesdB>\n'
    '      <BandWidth>300e6</BandWidth>\n'
    '      <ReferenceBandwidth>4000</ReferenceBandwidth>\n'
    '      <TransmitAntennaPatternFile>../input/example_antenna_patterns/isoflux_x_band.cut</TransmitAntennaPatternFile>\n'
    '      <PfdLimit>-150</PfdLimit>',
    propagator='SGP4'), tle='terrasarx.txt')

# ---------------------------------------------------------------- navigation
write_test('nav_dilution_of_precision', config(
    gps_constellation(), ground_segment(['Kourou']), grid_users(10), *FEB26, 300,
    '      <Type>nav_dilution_of_precision</Type>\n'
    '      <Direction>Ver</Direction>\n'
    '      <Statistic>Max</Statistic>\n'
    '      <Plot3D>True</Plot3D>\n'
    '      <ShowSatellite>False</ShowSatellite>'))

write_test('nav_accuracy', config(
    gps_constellation(uere=True), ground_segment(['Kourou']), grid_users(10), *FEB26, 300,
    '      <Type>nav_accuracy</Type>\n'
    '      <Direction>Hor</Direction>\n'
    '      <Statistic>Mean</Statistic>\n'
    '      <Plot3D>True</Plot3D>\n'
    '      <ShowSatellite>False</ShowSatellite>'))

# -------------------------------------------------------------------- power
write_test('sat_battery_depth_discharge', config(
    sso_constellation(ltan=True), ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 60,
    '      <Type>sat_battery_depth_discharge</Type>\n'
    '      <BatteryCapacityWh>500</BatteryCapacityWh>\n'
    '      <InitialSoC>1.0</InitialSoC>\n'
    '      <SolarPanelArea>5.0</SolarPanelArea>\n'
    '      <PanelEfficiency>0.3</PanelEfficiency>\n'
    '      <BasePowerDrawW>200</BasePowerDrawW>\n'
    '      <InstrumentPowerDrawW>400</InstrumentPowerDrawW>\n'
    '      <PayloadLatitudeLimit>60</PayloadLatitudeLimit>'))

write_test('sat_eclipse_duration', config(
    # RAAN 140 puts the orbit plane close to the Feb sun direction (beta ~ 0),
    # giving maximum-length eclipses every orbit; RAAN 50 would be eclipse-free
    sso_constellation(raan=140.0), ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-02-03 00:00:00', 30,
    '      <Type>sat_eclipse_duration</Type>'))

# ------------------------------------------------------------- data handling
write_test('sat_data_storage', config(
    sso_constellation(), ground_segment(['Svalbard', 'Inuvik']), static_users(DELFT), *FEB26, 30,
    '      <Type>sat_data_storage</Type>\n'
    '      <SSRCapacityGbits>2000</SSRCapacityGbits>\n'
    '      <InitialFillGbits>100</InitialFillGbits>\n'
    '      <InstrumentRateMbps>60</InstrumentRateMbps>\n'
    '      <DownlinkRateMbps>1070</DownlinkRateMbps>\n'
    '      <PayloadLatitudeLimit>60</PayloadLatitudeLimit>'))

# --------------------------------------------------------------------- orbit
def hpop_block(mass, area):
    return f"""    <HPOP>
      <IntegratorMinStep>0.001</IntegratorMinStep>
      <IntegratorMaxStep>300</IntegratorMaxStep>
      <IntegratorPositionTolerance>1.0</IntegratorPositionTolerance>
      <Mass>{mass}</Mass>
      <Geopotential>True</Geopotential>
      <GeopotentialDegree>21</GeopotentialDegree>
      <GeopotentialOrder>21</GeopotentialOrder>
      <EarthPoleRotation>True</EarthPoleRotation>
      <Drag>True</Drag>
      <DragArea>{area}</DragArea>
      <DragCd>2.2</DragCd>
      <DragModel>NRLMSISE00</DragModel>
      <SolarRadiationPressure>True</SolarRadiationPressure>
      <SRPArea>{area}</SRPArea>
      <SRPCr>1.5</SRPCr>
      <ThirdBodySun>True</ThirdBodySun>
      <ThirdBodyMoon>True</ThirdBodyMoon>
      <ThirdBodyPlanets>False</ThirdBodyPlanets>
      <SolidTides>True</SolidTides>
      <OceanTides>False</OceanTides>
      <Relativity>False</Relativity>
    </HPOP>
"""


write_test('orb_air_density', config(
    # 250 km LEO: day/night density variation and decay clearly visible
    sso_constellation(sma=6628137, incl=96.5, raan=140.0),
    ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 60,
    '      <Type>orb_air_density</Type>\n'
    '      <ConstellationID>1</ConstellationID>',
    propagator='HPOP', gr2sp=False, usr2sp=False,
    extra_sim=hpop_block(mass=500.0, area=3.2)))

write_test('orb_disturbance_forces', config(
    sso_constellation(sma=6628137, incl=96.5, raan=140.0),
    ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 60,
    '      <Type>orb_disturbance_forces</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>',
    propagator='HPOP', gr2sp=False, usr2sp=False,
    extra_sim=hpop_block(mass=500.0, area=3.2)))

write_test('orb_deltav_element', config(
    # 250 km LEO with strong drag: several altitude-keeping maneuvers in 2 days
    sso_constellation(sma=6628137, incl=96.5, raan=140.0),
    ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-02-03 00:00:00', 60,
    '      <Type>orb_deltav_element</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <TargetType>Altitude</TargetType>\n'
    '      <TargetValue>240000</TargetValue>\n'  # Mean altitude (below osculating 250 km)
    '      <DeadBand>1000</DeadBand>',
    propagator='HPOP', gr2sp=False, usr2sp=False,
    extra_sim=hpop_block(mass=500.0, area=3.2)))

# Pole wobble: 60 days so a visible arc of the Chandler+annual circle is
# traced; light force model to keep the long HPOP integration fast
HPOP_LIGHT = """    <HPOP>
      <IntegratorMaxStep>900</IntegratorMaxStep>
      <IntegratorPositionTolerance>10.0</IntegratorPositionTolerance>
      <Mass>500.0</Mass>
      <Geopotential>True</Geopotential>
      <GeopotentialDegree>8</GeopotentialDegree>
      <GeopotentialOrder>8</GeopotentialOrder>
      <EarthPoleRotation>True</EarthPoleRotation>
      <Drag>False</Drag>
      <SolarRadiationPressure>False</SolarRadiationPressure>
      <ThirdBodySun>False</ThirdBodySun>
      <ThirdBodyMoon>False</ThirdBodyMoon>
      <SolidTides>False</SolidTides>
      <OceanTides>False</OceanTides>
      <Relativity>False</Relativity>
    </HPOP>
"""
write_test('orb_pole_wobble', config(
    sso_constellation(), ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-04-02 00:00:00', 1800,
    '      <Type>orb_pole_wobble</Type>',
    propagator='HPOP', gr2sp=False, usr2sp=False, extra_sim=HPOP_LIGHT))

write_test('orb_kepler_elements', config(
    # 250 km LEO: strong drag makes the semi-major axis decay clearly visible
    sso_constellation(sma=6628137, incl=96.5, raan=140.0),
    ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-02-04 00:00:00', 60,
    '      <Type>orb_kepler_elements</Type>\n'
    '      <ConstellationID>1</ConstellationID>',
    propagator='HPOP', gr2sp=False, usr2sp=False,
    extra_sim=hpop_block(mass=500.0, area=3.2)))

write_test('orb_beta_angle', config(
    # ISS-like inclination propagated with SGP4 (J2 nodal regression) over
    # 4 months: the classic ~60-day beta beat cycle with eclipse-free peaks
    sso_constellation(sma=6798137, incl=51.6, raan=0.0,
                      extra='      <FrontalArea>3.2</FrontalArea>\n'
                            '      <Mass>500</Mass>\n'),
    ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-06-01 00:00:00', 3600,
    '      <Type>orb_beta_angle</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>',
    propagator='SGP4'))

write_test('orb_lifetime', config(
    # 500 km SSO: a few years of drag decay with the accelerating knee visible
    sso_constellation(sma=6878137, incl=97.4, raan=140.0),
    ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 60,
    '      <Type>orb_lifetime</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <Mass>500</Mass>\n'
    '      <DragArea>3.2</DragArea>\n'
    '      <DragCoefficient>2.2</DragCoefficient>\n'
    '      <DensityScale>1.0</DensityScale>\n'
    '      <MaxYears>100</MaxYears>'))

# Conjunction screening of Sentinel-1A against a frozen 265-object CelesTrak
# snapshot (both committed in the test folder, July 2026 epochs) - via
# CelestrakGroupFile the test needs no network and stays deterministic
write_test('orb_collision_check', config(
    tle_constellation('Sentinel-1A', 'mission_tle.txt', 'orb_collision_check'),
    ground_segment(['Svalbard']), static_users(DELFT),
    '2026-07-15 00:00:00', '2026-07-16 00:00:00', 300,
    '      <Type>orb_collision_check</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <CelestrakGroupFile>../tests/orb_collision_check/catalog_snapshot.txt</CelestrakGroupFile>\n'
    '      <ScreeningDistance>10000</ScreeningDistance>\n'
    '      <ScreeningStep>30</ScreeningStep>',
    propagator='SGP4'))

# Altitude-band neighbours only (the static sieve stage), reusing the frozen
# snapshot and mission TLE of the orb_collision_check test
write_test('orb_collision_alt_check', config(
    tle_constellation('Sentinel-1A', '../orb_collision_check/mission_tle.txt',
                      'orb_collision_alt_check'),
    ground_segment(['Svalbard']), static_users(DELFT),
    '2026-07-15 00:00:00', '2026-07-15 06:00:00', 600,
    '      <Type>orb_collision_alt_check</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <CelestrakGroupFile>../tests/orb_collision_check/catalog_snapshot.txt</CelestrakGroupFile>\n'
    '      <AltitudeMargin>10000</AltitudeMargin>'))

# Impulsive delta-v budget calculators on the satellite-block orbit: the time
# loop is irrelevant, so a short coarse window keeps the tests fast
write_test('orb_deltav_injection', config(
    sso_constellation(), ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-02-01 06:00:00', 600,
    '      <Type>orb_deltav_injection</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <Launcher>Ariane 62, 5000, 0.04</Launcher>\n'
    '      <Launcher>Vega-C, 15000, 0.15</Launcher>\n'
    '      <Launcher>Falcon 9, 15000, 0.10</Launcher>'))

write_test('orb_deltav_reentry', config(
    sso_constellation(), ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-02-01 06:00:00', 600,
    '      <Type>orb_deltav_reentry</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <IntermediatePerigee>250000</IntermediatePerigee>\n'
    '      <FinalPerigee>50000</FinalPerigee>'))

write_test('orb_deltav_collision', config(
    sso_constellation(), ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-02-01 06:00:00', 600,
    '      <Type>orb_deltav_collision</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <AvoidanceAltitude>10000</AvoidanceAltitude>'))

write_test('sat_drag_coefficient', config(
    # 500 km SSO with the built-in bus + solar panel model: free-molecular
    # Sentman panel drag with shadowing, attitude sweep at 30 deg for speed
    sso_constellation(sma=6878137, incl=97.4, raan=140.0),
    ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 300,
    '      <Type>sat_drag_coefficient</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <AccommodationCoefficient>0.93</AccommodationCoefficient>\n'
    '      <WallTemperature>300</WallTemperature>\n'
    '      <AttitudeStep>30</AttitudeStep>'))

write_test('orb_environment', config(
    # 800 km SSO: SAA proton crossings, outer-belt electron horns at high
    # latitude, and a representative AO/micrometeoroid environment
    sso_constellation(sma=7178137, incl=98.6, raan=140.0),
    ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 60,
    '      <Type>orb_environment</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <MissionYears>5</MissionYears>\n'
    '      <SurfaceArea>20</SurfaceArea>'))

write_test('hpop_benchmark', config(
    tle_constellation('TerraSAR-X', 'terrasarx.txt', 'hpop_benchmark'),
    ground_segment(['Kiruna']), static_users(DELFT), *FEB26, 60,
    '      <Type>cov_satellite_pvt</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>31698</SatelliteID>',
    propagator='HPOP', extra_sim=hpop_block(mass=1230.0, area=3.2)), tle='terrasarx.txt')

# Same scenario propagated with SGP4: generates the reference ("known") orbit
# file for the benchmark. Written next to the main Config.xml of the test.
_ref = config(
    tle_constellation('TerraSAR-X', 'terrasarx.txt', 'hpop_benchmark'),
    ground_segment(['Kiruna']), static_users(DELFT), *FEB26, 60,
    '      <Type>cov_satellite_pvt</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>31698</SatelliteID>',
    propagator='SGP4')
with open(os.path.join(TESTS, 'hpop_benchmark', 'Config_sgp4_reference.xml'), 'w') as f:
    f.write(_ref)
print('wrote hpop_benchmark/Config_sgp4_reference.xml')

# ------------------------------------------------------------- data handling
write_test('sat_data_latency', config(
    sso_constellation(), ground_segment(['Svalbard', 'Inuvik']), static_users(DELFT), *FEB26, 30,
    '      <Type>sat_data_latency</Type>\n'
    '      <SSRCapacityGbits>2000</SSRCapacityGbits>\n'
    '      <InitialFillGbits>100</InitialFillGbits>\n'
    '      <InstrumentRateMbps>60</InstrumentRateMbps>\n'
    '      <DownlinkRateMbps>1070</DownlinkRateMbps>\n'
    '      <PayloadLatitudeLimit>60</PayloadLatitudeLimit>\n'
    '      <GroundProcessingMin>15</GroundProcessingMin>'))

# --------------------------------------------------------- satellite platform
write_test('sat_thermal', config(
    # Same SSO/eclipse scenario as sat_eclipse_duration: the eclipses drive
    # the per-orbit temperature saw-tooth
    sso_constellation(raan=140.0), ground_segment(['Svalbard']), static_users(DELFT),
    *FEB26, 60,
    '      <Type>sat_thermal</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <SurfaceArea>6.0</SurfaceArea>\n'
    '      <CrossSectionSun>1.5</CrossSectionSun>\n'
    '      <CrossSectionEarth>1.5</CrossSectionEarth>\n'
    '      <Absorptivity>0.3</Absorptivity>\n'
    '      <Emissivity>0.8</Emissivity>\n'
    '      <InternalPowerW>300</InternalPowerW>\n'
    '      <HeatCapacity>50000</HeatCapacity>'))

write_test('sat_aocs', config(
    sso_constellation(raan=140.0), ground_segment(['Svalbard']), static_users(DELFT),
    *FEB26, 60,
    '      <Type>sat_aocs</Type>\n'
    '      <ConstellationID>1</ConstellationID>\n'
    '      <SatelliteID>1</SatelliteID>\n'
    '      <InertiaXX>100</InertiaXX>\n'
    '      <InertiaYY>120</InertiaYY>\n'
    '      <InertiaZZ>80</InertiaZZ>\n'
    '      <MaxPointingOffset>1.0</MaxPointingOffset>\n'
    '      <ResidualDipole>1.0</ResidualDipole>\n'
    '      <DragArea>2.5</DragArea>\n'
    '      <DragCoefficient>2.2</DragCoefficient>\n'
    '      <SrpArea>2.5</SrpArea>\n'
    '      <Reflectivity>0.6</Reflectivity>\n'
    '      <CopOffset>0.2</CopOffset>'))

# ------------------------------------------------- multiple analyses per run
# Three analyses in one simulation, incl. a repeated type (sky angles for two
# satellites) whose second output must come out numbered *_2.png; also the
# in-run mission report (<Report>) collecting all results in report.html
write_test('multi_analysis', config(
    gps_constellation(), ground_segment(['Kourou']),
    static_users([(52.0, 4.36), (1.35, 103.82)]), *FEB26, 300,
    ['      <Type>cov_satellite_visible</Type>',
     '      <Type>cov_satellite_sky_angles</Type>\n'
     '      <ConstellationID>1</ConstellationID>\n'
     '      <SatelliteID>1</SatelliteID>',
     '      <Type>cov_satellite_sky_angles</Type>\n'
     '      <ConstellationID>1</ConstellationID>\n'
     '      <SatelliteID>7</SatelliteID>'],
    extra_sim='    <Report>True</Report>\n'))

print('done')
