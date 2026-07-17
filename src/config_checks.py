"""
Validation of the Config.xml scenario file before anything is loaded: clear
error messages for missing/mistyped parameters instead of a mid-run crash,
and "did you mean" warnings for unknown tags (a misspelled optional tag would
otherwise be silently ignored). Called by main.load_configuration; errors are
fatal, warnings only logged.
"""
import difflib
import xml.etree.ElementTree as ET

from astropy.time import Time

import logging_svs as ls

# 3D-plot options shared by every world-map analysis (AnalysisPlot3D mixin)
PLOT3D_PARAMS = {'Plot3D', 'ShowSatellite', 'ShowOrbit', 'SatelliteModelFile',
                 'SatelliteModelScale', 'ModelRamAxis', 'ModelNadirAxis',
                 'MP4', 'EarthClouds', 'StationCones', 'SatelliteCone'}

# 2D-map decorations shared by every world-map analysis (AnalysisMap2D mixin)
MAP2D_PARAMS = {'ShowStations', 'ShowUsers', 'EarthImage', 'ShowGroundTrack',
                'Coastlines'}

COM_GR2SP_PARAMS = {
    'GroundStationID', 'TransmitterObject', 'CarrierFrequency', 'TransmitPowerW',
    'TransmitLossesdB', 'TransmitGaindB', 'TransmitAntennaPatternFile',
    'ReceiveAntennaPatternFile', 'SatelliteAntennaPointing', 'PExceedPerc',
    'IncludeRain', 'IncludeGas', 'IncludeScintillation', 'IncludeClouds',
    'ReceiveGaindB', 'ReceiveLossesdB', 'ReceiveTempK', 'ModulationType',
    'BitErrorRate', 'DataRateBitPerSec'}

# Known parameters per analysis type (the authoritative reference is readme.md)
ANALYSIS_PARAMS = {
    'cov_ground_track': {'ConstellationID', 'SatelliteID'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'cov_depth_of_coverage': PLOT3D_PARAMS,
    'cov_pass_time': {'ConstellationID', 'Statistic'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'cov_satellite_contour': {'ConstellationID', 'SatelliteID', 'ElevationMask'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'cov_satellite_highest': {'ConstellationID', 'Statistic'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'cov_satellite_pvt': {'ConstellationID', 'SatelliteID'},
    'cov_satellite_sky_angles': {'ConstellationID', 'SatelliteID'},
    'cov_satellite_visible': set(),
    'cov_satellite_visible_grid': {'Statistic'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'cov_satellite_visible_id': {'ConstellationID'},
    'obs_swath_conical': {'PolarView', 'Revisit', 'Statistic', 'SaveOutput'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'obs_swath_push_broom': {'PolarView', 'Revisit', 'Statistic', 'SaveOutput'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'obs_sza_push_broom': {'PolarView', 'Statistic', 'SaveOutput'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'obs_sza_subsat': {'PolarView', 'SaveOutput', 'RangeLatitude'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'obs_aoi_revisit': {'Statistic'} | MAP2D_PARAMS,
    'obs_target_imaging': {'Target', 'TargetFile', 'MaxOffNadir', 'MinSunElevation',
                           'ConstellationID'} | MAP2D_PARAMS,
    'com_gr2sp_budget': COM_GR2SP_PARAMS,
    'com_gr2sp_budget_interference': COM_GR2SP_PARAMS | {
        'BandWidth', 'TransmitGainManualdB', 'TransmitAntennaDiameter',
        'ReceiveAntennaDiameter'},
    'com_sp2sp_budget': {
        'SatelliteID1', 'SatelliteID2', 'CarrierFrequency', 'TransmitPowerW',
        'TransmitLossesdB', 'TransmitGaindB', 'ReceiveGaindB', 'ReceiveLossesdB',
        'ReceiveTempK', 'ModulationType', 'BitErrorRate', 'DataRateBitPerSec',
        'TransmitAntennaPatternFile', 'ReceiveAntennaPatternFile'},
    'com_doppler': {'StationID', 'CarrierFrequency'},
    'com_contact_plan': {'GroundStationID', 'ConstellationID', 'MinDuration',
                         'DownlinkRateMbps'},
    'com_pfd': {'GroundStationID', 'CarrierFrequency', 'TransmitPowerW',
                'TransmitLossesdB', 'BandWidth', 'ReferenceBandwidth',
                'TransmitGaindB', 'TransmitAntennaPatternFile', 'PfdLimit'},
    'nav_dilution_of_precision': {'Direction', 'Statistic'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'nav_accuracy': {'Direction', 'Statistic'} | PLOT3D_PARAMS | MAP2D_PARAMS,
    'sat_battery_depth_discharge': {
        'BatteryCapacityWh', 'InitialSoC', 'SolarPanelArea', 'PanelEfficiency',
        'BasePowerDrawW', 'InstrumentPowerDrawW', 'PayloadLatitudeLimit'},
    'sat_eclipse_duration': set(),
    'sat_data_storage': {'SSRCapacityGbits', 'InitialFillGbits', 'InstrumentRateMbps',
                         'DownlinkRateMbps', 'PayloadLatitudeLimit'},
    'sat_data_latency': {'SSRCapacityGbits', 'InitialFillGbits', 'InstrumentRateMbps',
                         'DownlinkRateMbps', 'PayloadLatitudeLimit', 'GroundProcessingMin'},
    'orb_altitude': {'ConstellationID', 'SatelliteID'},
    'orb_kepler_elements': {'ConstellationID', 'SatelliteID'},
    'orb_air_density': {'ConstellationID', 'SatelliteID'},
    'orb_disturbance_forces': {'ConstellationID', 'SatelliteID'},
    'orb_pole_wobble': set(),
    'orb_deltav_element': {'ConstellationID', 'SatelliteID', 'TargetType',
                           'TargetValue', 'DeadBand'},
    'orb_beta_angle': {'ConstellationID', 'SatelliteID'},
    'orb_lifetime': {'ConstellationID', 'SatelliteID', 'Mass', 'DragArea',
                     'DragCoefficient', 'DensityScale', 'MaxYears', 'ReentryAltitude'},
    'orb_environment': {'ConstellationID', 'SatelliteID', 'MissionYears',
                        'SurfaceArea'} | MAP2D_PARAMS,
    'orb_deltav_injection': {'ConstellationID', 'SatelliteID', 'Launcher'},
    'orb_deltav_reentry': {'ConstellationID', 'SatelliteID',
                           'IntermediatePerigee', 'FinalPerigee'},
    'orb_deltav_collision': {'ConstellationID', 'SatelliteID',
                             'AvoidanceAltitude'},
    'orb_collision_check': {'ConstellationID', 'SatelliteID', 'CelestrakGroup',
                            'CelestrakGroupFile', 'ScreeningDistance',
                            'ScreeningStep', 'HardBodyRadius', 'CovarianceRadial',
                            'CovarianceAlongTrack', 'CovarianceCrossTrack'},
    'orb_collision_alt_check': {'ConstellationID', 'SatelliteID', 'CelestrakGroup',
                                'CelestrakGroupFile', 'AltitudeMargin'},
    'sat_thermal': {'ConstellationID', 'SatelliteID', 'SurfaceArea', 'CrossSectionSun',
                    'CrossSectionEarth', 'Absorptivity', 'Emissivity', 'InternalPowerW',
                    'HeatCapacity', 'InitialTemperature'},
    'sat_aocs': {'ConstellationID', 'SatelliteID', 'InertiaXX', 'InertiaYY', 'InertiaZZ',
                 'MaxPointingOffset', 'ResidualDipole', 'DragArea', 'DragCoefficient',
                 'SrpArea', 'Reflectivity', 'CopOffset', 'WheelMomentum'},
    'sat_drag_coefficient': {'ConstellationID', 'SatelliteID', 'SatelliteModelFile',
                             'SatelliteModelScale', 'ModelRamAxis', 'ModelNadirAxis',
                             'AccommodationCoefficient',
                             'WallTemperature', 'ExosphericTemperature',
                             'AttitudeStep', 'Shadowing', 'MaxFacets'},
}

# Parameters that must be present per analysis type (their absence crashes the
# run or silently produces an empty result)
ANALYSIS_REQUIRED = {
    'cov_ground_track': ['ConstellationID'],
    'cov_pass_time': ['ConstellationID', 'Statistic'],
    'cov_satellite_contour': ['ConstellationID', 'ElevationMask'],
    'cov_satellite_highest': ['ConstellationID', 'Statistic'],
    'cov_satellite_pvt': ['ConstellationID'],
    'cov_satellite_sky_angles': ['ConstellationID', 'SatelliteID'],
    'cov_satellite_visible_grid': ['Statistic'],
    'cov_satellite_visible_id': ['ConstellationID'],
    'com_gr2sp_budget': ['GroundStationID', 'TransmitterObject', 'CarrierFrequency',
                         'TransmitPowerW', 'TransmitLossesdB', 'PExceedPerc',
                         'IncludeRain', 'IncludeGas', 'IncludeScintillation',
                         'IncludeClouds', 'ReceiveLossesdB', 'ReceiveTempK'],
    'com_gr2sp_budget_interference': ['GroundStationID', 'TransmitterObject',
                                      'CarrierFrequency', 'BandWidth', 'TransmitPowerW',
                                      'TransmitLossesdB', 'PExceedPerc', 'IncludeRain',
                                      'IncludeGas', 'IncludeScintillation',
                                      'IncludeClouds', 'ReceiveLossesdB', 'ReceiveTempK'],
    'com_sp2sp_budget': ['SatelliteID1', 'SatelliteID2', 'CarrierFrequency',
                         'TransmitPowerW', 'TransmitLossesdB', 'ReceiveLossesdB',
                         'ReceiveTempK'],
    'com_doppler': ['StationID', 'CarrierFrequency'],
    'com_pfd': ['CarrierFrequency', 'TransmitPowerW', 'BandWidth'],
    'nav_dilution_of_precision': ['Direction', 'Statistic'],
    'nav_accuracy': ['Direction', 'Statistic'],
    'obs_target_imaging': ['MaxOffNadir'],
    'sat_battery_depth_discharge': ['BatteryCapacityWh', 'InitialSoC', 'SolarPanelArea',
                                    'PanelEfficiency', 'BasePowerDrawW',
                                    'InstrumentPowerDrawW'],
    'sat_data_storage': ['SSRCapacityGbits', 'InitialFillGbits', 'InstrumentRateMbps',
                         'DownlinkRateMbps'],
    'sat_data_latency': ['SSRCapacityGbits', 'InitialFillGbits', 'InstrumentRateMbps',
                         'DownlinkRateMbps'],
    'orb_deltav_element': ['DeadBand'],
    'sat_thermal': ['SurfaceArea', 'HeatCapacity'],
    'sat_aocs': ['InertiaXX', 'InertiaYY', 'InertiaZZ'],
}

# One-of parameter groups per analysis type: at least one must be present
ANALYSIS_ANY_OF = {
    'com_gr2sp_budget': [('TransmitGaindB', 'TransmitAntennaPatternFile'),
                         ('ReceiveGaindB', 'ReceiveAntennaPatternFile')],
    'com_sp2sp_budget': [('TransmitGaindB', 'TransmitAntennaPatternFile'),
                         ('ReceiveGaindB', 'ReceiveAntennaPatternFile')],
    'com_gr2sp_budget_interference': [
        ('TransmitGaindB', 'TransmitGainManualdB', 'TransmitAntennaPatternFile'),
        ('ReceiveGaindB', 'ReceiveAntennaPatternFile'),
        ('TransmitAntennaDiameter', 'TransmitGainManualdB', 'TransmitAntennaPatternFile'),
        ('ReceiveAntennaDiameter', 'ReceiveAntennaPatternFile')],
    'com_pfd': [('TransmitGaindB', 'TransmitAntennaPatternFile')],
    'obs_target_imaging': [('Target', 'TargetFile')],
}

# Allowed children per structural block
SCHEMA = {
    'Scenario': {'SpaceSegment', 'GroundSegment', 'UserSegment', 'SimulationManager'},
    'SpaceSegment': {'Constellation'},
    'Constellation': {'ConstellationID', 'NumOfSatellites', 'NumOfPlanes',
                      'ConstellationName', 'ReceiverConstellation', 'TLEFileName',
                      'TLEFromCelestrak',
                      'ObsIncidenceAngleStart', 'ObsIncidenceAngleStop', 'ObsSwathStart',
                      'ObsSwathStop', 'ElevationMask', 'ElevationMaskMaximum', 'UERE',
                      'FrontalArea', 'Mass', 'Satellite'},
    'Satellite': {'SatelliteID', 'Plane', 'EpochMJD', 'Altitude', 'SemiMajorAxis',
                  'Eccentricity', 'Inclination', 'RAAN', 'LTAN', 'ArgOfPerigee',
                  'MeanAnomaly'},
    'GroundSegment': {'Network'},
    'Network': {'NumStation', 'NetworkName', 'GroundStation'},
    'GroundStation': {'Type', 'ConstellationID', 'GroundStationID', 'GroundStationName',
                      'Latitude', 'Longitude', 'Height', 'ReceiverConstellation',
                      'ElevationMask', 'ElevationMaskMaximum'},
    'UserSegment': {'User'},
    'User': {'Type', 'Latitude', 'Longitude', 'Height', 'ReceiverConstellation',
             'ElevationMask', 'ElevationMaskMaximum', 'LatMin', 'LatMax', 'LonMin',
             'LonMax', 'LatStep', 'LonStep', 'Name', 'PolygonList', 'PolygonFile',
             'TLEFileName'},
    'SimulationManager': {'StartDate', 'StopDate', 'TimeStep',
                          'IncludeStation2SpaceLinks', 'IncludeUser2SpaceLinks',
                          'IncludeSpace2SpaceLinks', 'OrbitsFromPreviousRun',
                          'OrbitPropagator', 'HPOP', 'Analysis', 'Report'},
    'HPOP': {'IntegratorMinStep', 'IntegratorMaxStep', 'IntegratorPositionTolerance',
             'Mass', 'Geopotential', 'GeopotentialDegree', 'GeopotentialOrder',
             'EarthPoleRotation', 'Drag', 'DragArea', 'DragCd', 'DragModel',
             'SolarRadiationPressure', 'SRPArea', 'SRPCr', 'ThirdBodySun',
             'ThirdBodyMoon', 'ThirdBodyPlanets', 'SolidTides', 'OceanTides',
             'OceanTidesDegree', 'OceanTidesOrder', 'Relativity'},
}

# Value type per tag name, for "is not a number"-style messages. floatlist is
# a comma-separated list of numbers; bits a string of 0/1 flags per constellation.
INT_TAGS = {'ConstellationID', 'NumOfSatellites', 'NumOfPlanes', 'SatelliteID',
            'Plane', 'GroundStationID', 'NumStation', 'GeopotentialDegree',
            'GeopotentialOrder', 'OceanTidesDegree', 'OceanTidesOrder',
            'SatelliteID1', 'SatelliteID2', 'StationID'}
FLOAT_TAGS = {'EpochMJD', 'Altitude', 'SemiMajorAxis', 'Eccentricity', 'Inclination',
              'RAAN', 'LTAN', 'ArgOfPerigee', 'MeanAnomaly', 'Latitude', 'Longitude',
              'Height', 'TimeStep', 'FrontalArea', 'Mass', 'ObsIncidenceAngleStart',
              'ObsIncidenceAngleStop', 'ObsSwathStart', 'ObsSwathStop', 'LatMin',
              'LatMax', 'LonMin', 'LonMax', 'LatStep', 'LonStep', 'PolarView',
              'CarrierFrequency', 'TransmitPowerW', 'TransmitLossesdB', 'TransmitGaindB',
              'ReceiveGaindB', 'ReceiveLossesdB', 'ReceiveTempK', 'PExceedPerc',
              'BitErrorRate', 'DataRateBitPerSec', 'BandWidth', 'TransmitAntennaDiameter',
              'ReceiveAntennaDiameter', 'BatteryCapacityWh', 'InitialSoC',
              'SolarPanelArea', 'PanelEfficiency', 'BasePowerDrawW',
              'InstrumentPowerDrawW', 'PayloadLatitudeLimit', 'SSRCapacityGbits',
              'InitialFillGbits', 'InstrumentRateMbps', 'DownlinkRateMbps',
              'GroundProcessingMin', 'TargetValue', 'DeadBand', 'SatelliteModelScale',
              'IntegratorMinStep', 'IntegratorMaxStep', 'IntegratorPositionTolerance',
              'DragArea', 'DragCd', 'SRPArea', 'SRPCr', 'SurfaceArea', 'CrossSectionSun',
              'CrossSectionEarth', 'Absorptivity', 'Emissivity', 'InternalPowerW',
              'HeatCapacity', 'InitialTemperature', 'InertiaXX', 'InertiaYY', 'InertiaZZ',
              'MaxPointingOffset', 'ResidualDipole', 'DragCoefficient', 'SrpArea',
              'Reflectivity', 'CopOffset', 'WheelMomentum', 'DensityScale', 'MaxYears',
              'ReentryAltitude', 'MinDuration', 'ReferenceBandwidth', 'PfdLimit',
              'MaxOffNadir', 'MinSunElevation', 'MissionYears',
              'AccommodationCoefficient', 'WallTemperature', 'ExosphericTemperature',
              'AttitudeStep', 'MaxFacets', 'IntermediatePerigee', 'FinalPerigee',
              'AvoidanceAltitude', 'ScreeningDistance', 'ScreeningStep',
              'AltitudeMargin', 'HardBodyRadius', 'CovarianceRadial',
              'CovarianceAlongTrack', 'CovarianceCrossTrack'}
BOOL_TAGS = {'IncludeStation2SpaceLinks', 'IncludeUser2SpaceLinks',
             'IncludeSpace2SpaceLinks', 'OrbitsFromPreviousRun', 'IncludeRain',
             'IncludeGas', 'IncludeScintillation', 'IncludeClouds', 'Revisit',
             'Plot3D', 'ShowSatellite', 'ShowOrbit', 'MP4', 'EarthClouds',
             'ShowStations', 'ShowUsers', 'EarthImage', 'ShowGroundTrack',
             'Coastlines', 'Report', 'Shadowing', 'StationCones', 'SatelliteCone',
             'Geopotential', 'EarthPoleRotation', 'Drag', 'SolarRadiationPressure',
             'ThirdBodySun', 'ThirdBodyMoon', 'ThirdBodyPlanets', 'SolidTides',
             'OceanTides', 'Relativity'}
FLOATLIST_TAGS = {'ElevationMask', 'ElevationMaskMaximum', 'UERE', 'RangeLatitude'}
BOOL_VALUES = {'yes', 'no', 'true', 'false', 't', 'f', '0', '1'}


def _text(node, tag):
    child = node.find(tag)
    return child.text.strip() if child is not None and child.text else None


def _suggest(tag, known):
    match = difflib.get_close_matches(tag, known, n=1)
    return f" (did you mean <{match[0]}>?)" if match else ''


def _check_children(node, path, allowed, warnings):
    for child in node:
        if child.tag not in allowed:
            warnings.append(f'{path}: unknown tag <{child.tag}> is ignored'
                            + _suggest(child.tag, allowed))


def _check_values(node, path, errors):
    """Type check of every known leaf tag inside node (one level deep)."""
    for child in node:
        text = (child.text or '').strip()
        where = f'{path}/{child.tag}'
        try:
            if child.tag in INT_TAGS:
                int(text)
            elif child.tag in FLOAT_TAGS:
                float(text)
            elif child.tag in FLOATLIST_TAGS:
                [float(v) for v in text.split(',')]
        except ValueError:
            errors.append(f"{where}: value '{text}' is not a number")
        if child.tag in BOOL_TAGS and text.lower() not in BOOL_VALUES:
            errors.append(f"{where}: value '{text}' is not a boolean (use True/False)")
        if child.tag == 'ReceiverConstellation' and (not text or set(text) - {'0', '1'}):
            errors.append(f"{where}: value '{text}' must be a string of 0/1 flags, "
                          f"one per constellation")


def _check_required(node, path, required, errors, context=''):
    for tag in required:
        if _text(node, tag) is None:
            errors.append(f'{path}: missing required <{tag}>{context}')


def validate_config(file_name):
    """Validate the scenario file; log warnings, exit on errors."""
    try:
        tree = ET.parse(file_name)
    except FileNotFoundError:
        ls.logger.error(f'Configuration file {file_name} not found')
        exit(1)
    except ET.ParseError as err:
        ls.logger.error(f'Configuration file {file_name} is not valid XML: {err}')
        exit(1)
    root = tree.getroot()
    errors, warnings = [], []

    if root.tag != 'Scenario':
        warnings.append(f'root element is <{root.tag}>, expected <Scenario>')
    _check_children(root, root.tag, SCHEMA['Scenario'], warnings)

    # --- Space segment ---------------------------------------------------
    constellation_ids = []
    constellations = list(root.iter('Constellation'))
    if not constellations:
        errors.append('SpaceSegment: no <Constellation> found - nothing to simulate')
    for i, const in enumerate(constellations):
        path = f'Constellation[{i + 1}]'
        _check_children(const, path, SCHEMA['Constellation'], warnings)
        _check_values(const, path, errors)
        _check_required(const, path, ['ConstellationID', 'NumOfSatellites',
                                      'NumOfPlanes', 'ConstellationName',
                                      'ReceiverConstellation'], errors)
        if _text(const, 'ConstellationID') is not None:
            try:
                constellation_ids.append(int(_text(const, 'ConstellationID')))
            except ValueError:
                pass
        satellites = list(const.iter('Satellite'))
        if not satellites and _text(const, 'TLEFileName') is None \
                and _text(const, 'TLEFromCelestrak') is None:
            errors.append(f'{path}: needs either <Satellite> block(s), a '
                          f'<TLEFileName> or a <TLEFromCelestrak>')
        for j, sat in enumerate(satellites):
            sat_path = f'{path}/Satellite[{j + 1}]'
            _check_children(sat, sat_path, SCHEMA['Satellite'], warnings)
            _check_values(sat, sat_path, errors)
            _check_required(sat, sat_path, ['SatelliteID', 'Plane', 'EpochMJD',
                                            'Eccentricity', 'ArgOfPerigee',
                                            'MeanAnomaly'], errors)
            if _text(sat, 'Altitude') is None and _text(sat, 'SemiMajorAxis') is None:
                errors.append(f'{sat_path}: needs <Altitude> or <SemiMajorAxis>')
            if _text(sat, 'LTAN') is None:  # SSO orbits derive inclination/RAAN from LTAN
                for tag in ('Inclination', 'RAAN'):
                    if _text(sat, tag) is None:
                        errors.append(f'{sat_path}: needs <{tag}> (or <LTAN> for an SSO orbit)')
            ecc = _text(sat, 'Eccentricity')
            try:
                if ecc is not None and not 0.0 <= float(ecc) < 1.0:
                    errors.append(f'{sat_path}: Eccentricity {ecc} outside [0, 1)')
            except ValueError:
                pass
    num_const = max(constellation_ids) if constellation_ids else 0

    def check_rx_constellation(node, path):
        rx = _text(node, 'ReceiverConstellation')
        if rx is not None and not set(rx) - {'0', '1'} and len(rx) < num_const:
            errors.append(f'{path}: ReceiverConstellation "{rx}" has {len(rx)} flag(s) '
                          f'but the space segment defines ConstellationID up to {num_const}')

    for i, const in enumerate(constellations):
        check_rx_constellation(const, f'Constellation[{i + 1}]')

    # --- Ground segment --------------------------------------------------
    for i, network in enumerate(root.iter('Network')):
        _check_children(network, f'Network[{i + 1}]', SCHEMA['Network'], warnings)
    for i, station in enumerate(root.iter('GroundStation')):
        path = f'GroundStation[{i + 1}]'
        _check_children(station, path, SCHEMA['GroundStation'], warnings)
        _check_values(station, path, errors)
        _check_required(station, path, ['ConstellationID', 'GroundStationID',
                                        'GroundStationName', 'Latitude', 'Longitude',
                                        'Height', 'ReceiverConstellation',
                                        'ElevationMask'], errors)
        check_rx_constellation(station, path)

    # --- User segment ----------------------------------------------------
    USER_REQUIRED = {  # Per user Type
        'Static': ['Latitude', 'Longitude', 'Height', 'ReceiverConstellation',
                   'ElevationMask'],
        'Grid': ['LatMin', 'LatMax', 'LonMin', 'LonMax', 'LatStep', 'LonStep',
                 'Height', 'ReceiverConstellation', 'ElevationMask'],
        'Polygon': ['Name', 'LatStep', 'LonStep', 'Height', 'ReceiverConstellation',
                    'ElevationMask'],
        'Spacecraft': ['TLEFileName', 'ReceiverConstellation', 'ElevationMask'],
    }
    for i, user in enumerate(root.iter('User')):
        path = f'User[{i + 1}]'
        _check_children(user, path, SCHEMA['User'], warnings)
        _check_values(user, path, errors)
        user_type = _text(user, 'Type')
        if user_type is None:
            errors.append(f'{path}: missing required <Type> '
                          f'(Static, Grid, Polygon or Spacecraft)')
        elif user_type not in USER_REQUIRED:
            errors.append(f"{path}: unknown Type '{user_type}' "
                          f"(use Static, Grid, Polygon or Spacecraft)")
        else:
            _check_required(user, path, USER_REQUIRED[user_type], errors,
                            context=f' (Type {user_type})')
            if user_type == 'Polygon' and _text(user, 'PolygonList') is None \
                    and _text(user, 'PolygonFile') is None:
                errors.append(f'{path}: Polygon user needs <PolygonList> or <PolygonFile>')
        check_rx_constellation(user, path)

    # --- Simulation manager ----------------------------------------------
    sims = list(root.iter('SimulationManager'))
    if not sims:
        errors.append('missing <SimulationManager> block')
    for sim in sims:
        path = 'SimulationManager'
        _check_children(sim, path, SCHEMA['SimulationManager'], warnings)
        _check_values(sim, path, errors)
        _check_required(sim, path, ['StartDate', 'StopDate', 'TimeStep',
                                    'IncludeStation2SpaceLinks', 'IncludeUser2SpaceLinks',
                                    'IncludeSpace2SpaceLinks', 'OrbitsFromPreviousRun',
                                    'OrbitPropagator'], errors)
        times = {}
        for tag in ('StartDate', 'StopDate'):
            text = _text(sim, tag)
            if text is not None:
                try:
                    times[tag] = Time(text, scale='utc').mjd
                except ValueError:
                    errors.append(f"{path}/{tag}: '{text}' is not a valid date "
                                  f"(use e.g. 2026-02-01 00:00:00)")
        if len(times) == 2 and times['StopDate'] <= times['StartDate']:
            errors.append(f'{path}: StopDate must be after StartDate')
        step = _text(sim, 'TimeStep')
        try:
            if step is not None and float(step) <= 0:
                errors.append(f'{path}: TimeStep must be positive, got {step}')
        except ValueError:
            pass
        propagator = _text(sim, 'OrbitPropagator')
        if propagator is not None and propagator not in ('Keplerian', 'SGP4', 'HPOP'):
            errors.append(f"{path}: unknown OrbitPropagator '{propagator}' "
                          f"(use Keplerian, SGP4 or HPOP)")
        hpop = sim.find('HPOP')
        if hpop is not None:
            _check_children(hpop, f'{path}/HPOP', SCHEMA['HPOP'], warnings)
            _check_values(hpop, f'{path}/HPOP', errors)

    # --- Analyses ----------------------------------------------------------
    for i, analysis in enumerate(root.iter('Analysis')):
        type_str = _text(analysis, 'Type')
        path = f'Analysis[{i + 1}]'
        if type_str is None:
            errors.append(f'{path}: missing required <Type>')
            continue
        path = f'{path} ({type_str})'
        if type_str not in ANALYSIS_PARAMS:
            errors.append(f'{path}: unknown analysis type'
                          + _suggest(type_str, list(ANALYSIS_PARAMS)))
            continue
        _check_children(analysis, path, ANALYSIS_PARAMS[type_str] | {'Type'}, warnings)
        _check_values(analysis, path, errors)
        _check_required(analysis, path, ANALYSIS_REQUIRED.get(type_str, []), errors)
        for group in ANALYSIS_ANY_OF.get(type_str, []):
            if all(_text(analysis, tag) is None for tag in group):
                errors.append(f'{path}: needs one of ' +
                              ' / '.join(f'<{tag}>' for tag in group))
        # Analysis ConstellationID must reference a defined constellation
        const_id = _text(analysis, 'ConstellationID')
        try:
            if const_id is not None and constellation_ids and \
                    int(const_id) not in constellation_ids:
                errors.append(f'{path}: ConstellationID {const_id} not defined in the '
                              f'space segment (defined: {sorted(constellation_ids)})')
        except ValueError:
            pass
        if type_str in ('obs_swath_conical', 'obs_swath_push_broom') and \
                misc_str2bool(_text(analysis, 'Revisit')) and \
                _text(analysis, 'Statistic') is None:
            errors.append(f'{path}: <Revisit>True</Revisit> needs a <Statistic> '
                          f'(min, mean, max, std or median)')

    for warning in warnings:
        ls.logger.warning(f'Config check: {warning}')
    if errors:
        for error in errors:
            ls.logger.error(f'Config check: {error}')
        ls.logger.error(f'Configuration {file_name} invalid: {len(errors)} error(s), '
                        f'see messages above (readme.md documents all parameters)')
        exit(1)
    ls.logger.info(f'Configuration {file_name} validated: '
                   f'{len(warnings)} warning(s), no errors')


def misc_str2bool(text):
    return text is not None and text.lower() in ('yes', 'true', 't', '1')
