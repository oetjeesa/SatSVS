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
    '      <Plot3D>True</Plot3D>',
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
    '      <SatelliteModelScale>500000</SatelliteModelScale>'))

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
    '      <Plot3D>True</Plot3D>',
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

# ------------------------------------------------------------- communication
write_test('com_gr2sp_budget', config(
    tle_constellation('TerraSAR-X', 'terrasarx.txt', 'com_gr2sp_budget'),
    ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 60,
    '      <Type>com_gr2sp_budget</Type>\n'
    '      <GroundStationID>1</GroundStationID>\n'
    '      <TransmitterObject>Satellite</TransmitterObject>\n'
    '      <CarrierFrequency>8.025e9</CarrierFrequency>\n'
    '      <TransmitPowerW>10</TransmitPowerW>\n'
    '      <TransmitLossesdB>2</TransmitLossesdB>\n'
    '      <TransmitGaindB>20</TransmitGaindB>\n'
    '      <ReceiveGaindB>45</ReceiveGaindB>\n'
    '      <ReceiveLossesdB>3</ReceiveLossesdB>\n'
    '      <ReceiveTempK>200</ReceiveTempK>\n'
    '      <PExceedPerc>0.5</PExceedPerc>\n'
    '      <IncludeGas>True</IncludeGas>\n'
    '      <IncludeRain>True</IncludeRain>\n'
    '      <IncludeScintillation>False</IncludeScintillation>\n'
    '      <IncludeClouds>False</IncludeClouds>\n'
    '      <ModulationType>QPSK</ModulationType>\n'
    '      <BitErrorRate>1e-5</BitErrorRate>\n'
    '      <DataRateBitPerSec>300e6</DataRateBitPerSec>',
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
    '      <ReceiveGaindB>64</ReceiveGaindB>\n'
    '      <ReceiveAntennaDiameter>6.8</ReceiveAntennaDiameter>\n'
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
write_test('pow_battery_depth_discharge', config(
    sso_constellation(ltan=True), ground_segment(['Svalbard']), static_users(DELFT), *FEB26, 60,
    '      <Type>pow_battery_depth_discharge</Type>\n'
    '      <BatteryCapacityWh>500</BatteryCapacityWh>\n'
    '      <InitialSoC>1.0</InitialSoC>\n'
    '      <SolarPanelArea>5.0</SolarPanelArea>\n'
    '      <PanelEfficiency>0.3</PanelEfficiency>\n'
    '      <BasePowerDrawW>200</BasePowerDrawW>\n'
    '      <InstrumentPowerDrawW>400</InstrumentPowerDrawW>\n'
    '      <PayloadLatitudeLimit>60</PayloadLatitudeLimit>'))

write_test('pow_eclipse_duration', config(
    # RAAN 140 puts the orbit plane close to the Feb sun direction (beta ~ 0),
    # giving maximum-length eclipses every orbit; RAAN 50 would be eclipse-free
    sso_constellation(raan=140.0), ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-02-03 00:00:00', 30,
    '      <Type>pow_eclipse_duration</Type>'))

# ------------------------------------------------------------- data handling
write_test('dat_storage', config(
    sso_constellation(), ground_segment(['Svalbard', 'Inuvik']), static_users(DELFT), *FEB26, 30,
    '      <Type>dat_storage</Type>\n'
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


write_test('orb_kepler_elements', config(
    # 250 km LEO: strong drag makes the semi-major axis decay clearly visible
    sso_constellation(sma=6628137, incl=96.5, raan=140.0),
    ground_segment(['Svalbard']), static_users(DELFT),
    '2026-02-01 00:00:00', '2026-02-04 00:00:00', 60,
    '      <Type>orb_kepler_elements</Type>\n'
    '      <ConstellationID>1</ConstellationID>',
    propagator='HPOP', gr2sp=False, usr2sp=False,
    extra_sim=hpop_block(mass=500.0, area=3.2)))

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
write_test('dat_latency', config(
    sso_constellation(), ground_segment(['Svalbard', 'Inuvik']), static_users(DELFT), *FEB26, 30,
    '      <Type>dat_latency</Type>\n'
    '      <SSRCapacityGbits>2000</SSRCapacityGbits>\n'
    '      <InitialFillGbits>100</InitialFillGbits>\n'
    '      <InstrumentRateMbps>60</InstrumentRateMbps>\n'
    '      <DownlinkRateMbps>1070</DownlinkRateMbps>\n'
    '      <PayloadLatitudeLimit>60</PayloadLatitudeLimit>\n'
    '      <GroundProcessingMin>15</GroundProcessingMin>'))

# ------------------------------------------------- multiple analyses per run
# Three analyses in one simulation, incl. a repeated type (sky angles for two
# satellites) whose second output must come out numbered *_2.png
write_test('multi_analysis', config(
    gps_constellation(), ground_segment(['Kourou']),
    static_users([(52.0, 4.36), (1.35, 103.82)]), *FEB26, 300,
    ['      <Type>cov_satellite_visible</Type>',
     '      <Type>cov_satellite_sky_angles</Type>\n'
     '      <ConstellationID>1</ConstellationID>\n'
     '      <SatelliteID>1</SatelliteID>',
     '      <Type>cov_satellite_sky_angles</Type>\n'
     '      <ConstellationID>1</ConstellationID>\n'
     '      <SatelliteID>7</SatelliteID>']))

print('done')
