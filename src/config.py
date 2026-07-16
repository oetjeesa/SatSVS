import xml.etree.ElementTree as ET
from math import ceil, radians
from astropy.time import Time
from sgp4.api import Satrec, WGS84
import os
import matplotlib.pyplot as plt
from geopandas import GeoSeries, GeoDataFrame
from shapely.geometry import Point, Polygon
import numpy as np
import ast

# Project modules
from constants import *
from analysis_cov import *
from analysis_obs import *
from analysis_com import *
from analysis_nav import *
from analysis_orb import *
from analysis_sat import *  # sat_ platform subsystem analyses (thermal, AOCS, power, data)

from segments import Constellation, Satellite, Station, User, Ground2SpaceLink, User2SpaceLink, Space2SpaceLink
import logging_svs as ls
import misc_fn
import copy


fetch_celestrak_tle = misc_fn.fetch_celestrak_tle  # moved to misc_fn (also used by orb_collision_check)


class AppConfig:
    
    def __init__(self, file_name=None, output_dir='../output'):

        self.output_dir = output_dir  # Directory for plots, data dumps and main.log
        if file_name is not None:  # File references inside the config may be
            # relative to the config file itself (see misc_fn.resolve_path)
            misc_fn.set_config_dir(os.path.dirname(os.path.abspath(file_name)))

        self.constellations = []
        self.satellites = []
        self.stations = []
        self.users = []
        self.gr2sp = []
        self.usr2sp = []
        self.sp2sp = []

        self.analyses = []  # One analysis instance per <Analysis> block, each with its own metric memory
        self.file_name = file_name

        self.start_time = 0  # MJD start time simulation
        self.stop_time = 0  # MJD stop time simulation
        self.time_step = 0  # in seconds
        self.num_epoch = 0  # Number of epochs in simulation
        self.cnt_epoch = 0
        self.time_gmst = 0  # Loop time
        self.time_mjd = 0  # Loop time
        self.time_str = ''  # Loop time
        self.time_datetime = None  # Loop time

        self.include_gr2sp = True
        self.include_usr2sp = True
        self.include_sp2sp = True
        self.orbits_from_previous_run = False
        self.report = False  # <Report>: write report.html at the end of the run
        self.data_orbits = []  # Orbits from previous run
        self.orbit_propagator = ''  # Text with either Keplerian / SGP4 / HPOP
        self.hpop_config = None  # Parsed <HPOP> block (HPOP propagator only)
        self.hpop = None  # propagation_hpop.HpopPropagation instance (HPOP only)
        self.num_constellation = 0
        self.num_sat = 0
        self.num_station = 0
        self.num_user = 0
        self.num_sp2sp = 0
        self.num_st2sp = 0
        self.num_usr2sp = 0
        self.user_latitudes = None
        self.user_longitudes = None

        self.file_orbits = None  # Handle of the orbit cache file, open during the time loop
        self.sat_pos_ecf = None  # (num_sat, 3) satellite ECF positions for the batched link kernels
        self.time_jd = 0  # Loop time, julian date split in day + fraction (for sgp4)
        self.time_fr = 0
        self.times_mjd_pre = None  # Precomputed per-epoch times, see main.precompute_times
        self.times_gmst_pre = None
        self.times_jd_pre = None
        self.times_fr_pre = None
        self.times_str_pre = None
        self.times_datetime_pre = None
        self.times_f_doy_pre = None

    def output_path(self, file_name):
        """Full path of an output file (plot, data dump, cache) in the run's
        output directory (--output-dir, default ../output)."""
        return os.path.join(self.output_dir, file_name)

    def load_satellites(self):

        # Get the list of constellations
        tree = ET.parse(self.file_name)
        root = tree.getroot()
        for constellation in root.iter('Constellation'):
            const = Constellation()
            # Mandatory
            const.constellation_id = int(constellation.find('ConstellationID').text)
            const.num_sat = int(constellation.find('NumOfSatellites').text)
            const.num_planes = int(constellation.find('NumOfPlanes').text)
            const.constellation_name = constellation.find('ConstellationName').text
            const.rx_constellation = constellation.find('ReceiverConstellation').text

            # Now some optional parameters
            if constellation.find('TLEFileName') is not None:
                const.tle_file_name = constellation.find('TLEFileName').text
            if constellation.find('TLEFromCelestrak') is not None:
                # Force-download the latest TLE from CelesTrak before the run;
                # the identifier is a NORAD catalog number or a satellite name
                identifier = constellation.find('TLEFromCelestrak').text.strip()
                if const.tle_file_name is not None:
                    dest = misc_fn.resolve_path(const.tle_file_name)
                else:  # No TLE file configured: cache next to the Config.xml
                    safe = ''.join(c if c.isalnum() or c in '-_' else '_'
                                   for c in identifier)
                    dest = os.path.join(misc_fn.config_dir(),
                                        f'tle_celestrak_{safe}.txt')
                    const.tle_file_name = dest
                if fetch_celestrak_tle(identifier, dest):
                    ls.logger.info(f'Downloaded latest TLE for "{identifier}" '
                                   f'from CelesTrak into {dest}')
                elif os.path.isfile(dest):
                    ls.logger.warning(f'Using the existing TLE file {dest} instead')
                else:
                    ls.logger.error(f'No TLE available for "{identifier}": the '
                                    f'CelesTrak download failed and there is no '
                                    f'previously downloaded file {dest}')
                    exit()
            if constellation.find('ObsIncidenceAngleStart') is not None:
                const.obs_inci_angle_start = radians(float(constellation.find('ObsIncidenceAngleStart').text))
            if constellation.find('ObsIncidenceAngleStop') is not None:
                const.obs_inci_angle_stop = radians(float(constellation.find('ObsIncidenceAngleStop').text))
            if constellation.find('ObsSwathStart') is not None:
                const.obs_swath_start = float(constellation.find('ObsSwathStart').text)
            if constellation.find('ObsSwathStop') is not None:
                const.obs_swath_stop = float(constellation.find('ObsSwathStop').text)
            if constellation.find('ElevationMask') is not None:
                mask_values = constellation.find('ElevationMask').text.split(',')
                for mask_value in mask_values:
                    const.elevation_mask.append(radians(float(mask_value)))
                if constellation.find('ElevationMaskMaximum') is not None: # Optional
                    mask_max_values = constellation.find('ElevationMaskMaximum').text.split(',')
                    for mask_max_value in mask_max_values:
                        const.el_mask_max.append(radians(mask_max_value))
                else:  # If no maximum set to 90 degrees
                    const.el_mask_max = len(mask_values)*[radians(90.0)]
            if constellation.find('UERE') is not None:
                uere_values = constellation.find('UERE').text.split(',')
                for uere_value in uere_values:
                    const.uere_list.append(float(uere_value))
            if constellation.find('FrontalArea') is not None:
                const.frontal_area = float(constellation.find('FrontalArea').text)
            if constellation.find('Mass') is not None:
                const.mass = float(constellation.find('Mass').text)
                # B* = Cd*A/m*rho_0/2 , which is an SGP4-type drag coefficient, according to Celestrak
                # if rho_0 is give as 0.1570kg/m2/earth_radii then B* has units of (earth radii)-1 as it should be in sgp4.
                # typically Cd is between 2 and 2.5 for different shapes, boxes, plates, etc. lets say 2.2
                # However the NASA website says that B* is more a radiation pressure coefficient
                # It is better to set it to 0
                # const.bstar = 2.2 * const.frontal_area/const.mass*0.1570/2
                const.bstar = 0.0

            ls.logger.info(const.__dict__)
            self.constellations.append(const)

            for satellite in constellation.iter('Satellite'):
                sat = Satellite()
                sat.sat_id = int(satellite.find('SatelliteID').text)
                sat.plane = int(satellite.find('Plane').text)
                sat.constellation_id = const.constellation_id
                sat.constellation_name = const.constellation_name
                sat.rx_constellation = const.rx_constellation
                sat.elevation_mask = const.elevation_mask
                sat.el_mask_max = const.el_mask_max
                sat.uere_list = const.uere_list
                sat.kepler.epoch_mjd = float(satellite.find('EpochMJD').text)
                if satellite.find('Altitude') is not None:
                    sat.kepler.semi_major_axis = float(satellite.find('Altitude').text)+R_EARTH
                else:
                    sat.kepler.semi_major_axis = float(satellite.find('SemiMajorAxis').text)
                sat.kepler.eccentricity = float(satellite.find('Eccentricity').text)
                if satellite.find('LTAN') is not None:
                    sat.kepler.inclination =  math.acos(-pow(sat.kepler.semi_major_axis/12352000,7/2))
                    sat.ltan = float(satellite.find('LTAN').text)
                else:
                    sat.kepler.inclination = radians(float(satellite.find('Inclination').text))
                    sat.kepler.right_ascension = radians(float(satellite.find('RAAN').text))
                sat.kepler.arg_perigee = radians(float(satellite.find('ArgOfPerigee').text))
                sat.kepler.mean_anomaly = radians(float(satellite.find('MeanAnomaly').text))
                if const.frontal_area is not None:
                    sat.frontal_area = const.frontal_area
                    sat.mass = const.mass
                    # sgp4 v2 accelerated (C) propagator, initialised from the Kepler elements
                    sat.satrec = Satrec()
                    sat.satrec.sgp4init(
                        WGS84, 'i', sat.sat_id,
                        sat.kepler.epoch_mjd + 2400000.5 - 2433281.5,  # epoch as days since 1949 December 31 00:00 UT
                        const.bstar, 0.0, 0.0,  # bstar, ndot, nddot
                        sat.kepler.eccentricity, sat.kepler.arg_perigee, sat.kepler.inclination,
                        sat.kepler.mean_anomaly,
                        np.sqrt(GM_EARTH/np.power(sat.kepler.semi_major_axis, 3))*60.0,  # no_kozai in [rad/min]
                        sat.kepler.right_ascension)
                ls.logger.info(sat.__dict__)
                self.satellites.append(sat)
            # TLE files to load: the <TLEFileName> element(s), or the file
            # downloaded from CelesTrak when only <TLEFromCelestrak> is given
            tle_files = [item.text for item in constellation.iter('TLEFileName')]
            if not tle_files and const.tle_file_name is not None:
                tle_files = [const.tle_file_name]
            for tle_file in tle_files:
                with open(misc_fn.resolve_path(tle_file)) as file:
                    cnt = 0
                    for line in file:
                        if cnt == 0:
                            sat = Satellite()
                            sat.Name = line
                        if cnt == 1:
                            sat.tle_line1 = line
                            sat.sat_id = int(line[2:7])  # NORAD Catalog Number
                            sat.constellation_id = const.constellation_id
                            sat.constellation_name = const.constellation_name
                            full_year = "20" + line[18:20]
                            epoch_yr = int(full_year)
                            epoch_doy = float(line[20:31])
                            sat.kepler.epoch_mjd = misc_fn.yyyy_doy2mjd(epoch_yr, epoch_doy)
                        if cnt == 2:
                            sat.tle_line2 = line
                            mean_motion = float(line[52:59]) / 86400 * 2 * PI  # rad / s
                            sat.kepler.semi_major_axis = pow((GM_EARTH/(pow(mean_motion, 2))), (1.0 / 3.0))
                            eccentricity = "0." + line[26:33]
                            sat.kepler.eccentricity = float(eccentricity)
                            sat.kepler.inclination = radians(float(line[8:15]))
                            sat.kepler.right_ascension = radians(float(line[17:24]))
                            sat.kepler.arg_perigee = radians(float(line[34:41]))
                            sat.kepler.mean_anomaly = radians(float(line[43:50]))
                            ls.logger.info(['Found satellite:', str(sat.sat_id), 'Kepler Elements', str(sat.kepler.__dict__)])
                            sat.rx_constellation = const.rx_constellation
                            sat.elevation_mask = const.elevation_mask
                            sat.el_mask_max = const.el_mask_max
                            sat.satrec = Satrec.twoline2rv(sat.tle_line1, sat.tle_line2, WGS84)
                            self.satellites.append(sat)
                            cnt = -1  # reset so that next tle set starts at cnt=0
                        cnt += 1
        self.num_sat = len(self.satellites)
        self.num_constellation = len(self.constellations)

    def load_stations(self):
        # Get the list of constellations
        tree = ET.parse(self.file_name)
        root = tree.getroot()
        for station_el in root.iter('GroundStation'):
            station = Station()
            station.constellation_id = int(station_el.find('ConstellationID').text)
            station.station_id = int(station_el.find('GroundStationID').text)
            station.station_name = station_el.find('GroundStationName').text
            station.lla[0] = radians(float(station_el.find('Latitude').text))
            station.lla[1] = radians(float(station_el.find('Longitude').text))
            station.lla[2] = float(station_el.find('Height').text)
            station.rx_constellation = station_el.find('ReceiverConstellation').text

            mask_values = station_el.find('ElevationMask').text.split(',')
            for mask_value in mask_values:
                station.elevation_mask.append(radians(float(mask_value)))
            if station_el.find('ElevationMaskMaximum') is not None: # Optional
                mask_max_values = station_el.find('ElevationMaskMaximum').text.split(',')
                for mask_max_value in mask_max_values:
                    station.el_mask_max.append(radians(mask_max_value))
            else:  # If no maximum set to 90 degrees
                station.el_mask_max = len(mask_values)*[radians(90.0)]

            station.det_posvel_ecf()  # do it once to establish ECF coordinates
            ls.logger.info(station.__dict__)
            self.stations.append(station)

        self.num_station = len(self.stations)

    def load_users(self):
        # Get the list of constellations
        tree = ET.parse(self.file_name)
        root = tree.getroot()
        for user_element in root.iter('User'):

            if user_element.find('Type').text == 'Static':
                user = User()
                user.type = user_element.find('Type').text
                user.user_id = len(self.users)
                user.lla[0] = radians(float(user_element.find('Latitude').text))
                user.lla[1] = radians(float(user_element.find('Longitude').text))
                user.lla[2] = float(user_element.find('Height').text)
                user.rx_constellation = user_element.find('ReceiverConstellation').text
                mask_values = user_element.find('ElevationMask').text.split(',')
                user.elevation_mask = [radians(float(n)) for n in mask_values]
                if user_element.find('ElevationMaskMaximum') is not None:
                    mask_max_values = user_element.find('ElevationMaskMaximum').text.split(',')
                    user.el_mask_max = [radians(float(n))for n in mask_max_values]
                else:
                    user.el_mask_max = len(mask_values) * [radians(90.0)]
                user.det_posvel_ecf()
                self.users.append(user)

            if user_element.find('Type').text == 'Grid':
                lat_min = float(user_element.find('LatMin').text)
                lat_max = float(user_element.find('LatMax').text)
                lon_min = float(user_element.find('LonMin').text)
                lon_max = float(user_element.find('LonMax').text)
                lat_step = float(user_element.find('LatStep').text)
                lon_step = float(user_element.find('LonStep').text)
                num_lat = int((lat_max - lat_min) / lat_step) + 1
                num_lon = int((lon_max - lon_min) / lon_step) + 1
                self.user_latitudes = np.linspace(lat_min,lat_max,num_lat)
                self.user_longitudes = np.linspace(lon_min,lon_max,num_lon)
                height = float(user_element.find('Height').text)
                rx_constellation = user_element.find('ReceiverConstellation').text
                mask_values = user_element.find('ElevationMask').text.split(',')
                mask_max_values = []
                if user_element.find('ElevationMaskMaximum') is not None:
                    mask_max_values = user_element.find('ElevationMaskMaximum').text.split(',')

                # Now make the full list of users
                user = User()
                user.type = "Grid"
                user.UserName = "Grid"
                user.num_lat = num_lat
                user.num_lon = num_lon
                user.rx_constellation = rx_constellation
                user.elevation_mask = [radians(float(n)) for n in mask_values]
                user.lla[2] = height
                if user_element.find('ElevationMaskMaximum') is not None:
                    user.el_mask_max = [radians(float(n)) for n in mask_max_values]
                else:
                    user.el_mask_max = len(mask_values) * [radians(90.0)]
                for i in range(num_lat):
                    latitude = lat_min + lat_step * i
                    for j in range(num_lon):
                        longitude = lon_min + lon_step * j
                        user.lla[0] = radians(latitude)
                        user.lla[1] = radians(longitude)
                        user.det_posvel_ecf()  # Do it once to initialise
                        user.user_id = len(self.users)
                        self.users.append(copy.deepcopy(user))

            if user_element.find('Type').text == 'Polygon':
                name = user_element.find('Name').text
                lat_step = float(user_element.find('LatStep').text)
                lon_step = float(user_element.find('LonStep').text)
                height = float(user_element.find('Height').text)
                rx_constellation = user_element.find('ReceiverConstellation').text
                mask_values = user_element.find('ElevationMask').text.split(',')
                mask_max_values = []
                if user_element.find('ElevationMaskMaximum') is not None:
                    mask_max_values = user_element.find('ElevationMaskMaximum').text.split(',')
                if user_element.find('PolygonList') is not None:
                    points = list(ast.literal_eval(user_element.find('PolygonList').text))
                    poly = Polygon(points)
                else:
                    poly = GeoDataFrame.from_file(misc_fn.resolve_path(user_element.find('PolygonFile').text))
                    poly = poly['geometry'].iloc[0]
                    #poly = poly[poly.ADMIN == 'Denmark']['geometry'].iloc[0]
                # Setup a rough grid
                xmin, xmax, ymin, ymax = poly.bounds[0], poly.bounds[2], poly.bounds[1], poly.bounds[3]
                ls.logger.info(f'User polygon {name} bounds: lon: {round(xmin,1)},{round(xmax,1)} [deg] lat: {round(ymin,1)},{round(ymax,1)} [deg]')
                ls.logger.info(f'User polygon {name} area: {round(poly.area,1)} [deg^2] percentage: {round(poly.area/360/180*100,1)} [%] Earth surface')
                xx, yy = np.meshgrid(np.arange(xmin, xmax, lon_step), np.arange(ymin, ymax, lat_step))
                xc = xx.flatten(); yc = yy.flatten()
                # Check the ones within the polygon
                pts = GeoSeries([Point(x, y) for x, y in zip(xc, yc)])
                in_map = np.array([pts.within(poly)]).sum(axis=0)
                pts = GeoSeries([val for pos, val in enumerate(pts) if in_map[pos]])

                # Now make the full list of users
                for i, pt in enumerate(pts):
                    user = User()
                    user.type = "Polygon"
                    user.name = name
                    user.user_id = i
                    user.UserName = "Polygon"
                    user.lla[0] = radians(pt.y)
                    user.lla[1] = radians(pt.x)
                    user.lla[2] = height
                    user.num_lat = len(pts)
                    user.num_lon = len(pts)
                    user.rx_constellation = rx_constellation
                    user.elevation_mask = [radians(float(n)) for n in mask_values]
                    if user_element.find('ElevationMaskMaximum') is not None:
                        user.el_mask_max = [radians(float(n)) for n in mask_max_values]
                    else:
                        user.el_mask_max = len(mask_values) * [radians(90.0)]
                    user.det_posvel_ecf()  # Do it once to initialise
                    self.users.append(user)

            if user_element.find('Type').text == 'Spacecraft':
                user = User()
                user.user_id = 0
                user.type = 'Spacecraft'
                user.rx_constellation = user_element.find('ReceiverConstellation').text
                mask_values = user_element.find('ElevationMask').text.split(',')
                user.elevation_mask = [radians(float(n)) for n in mask_values]
                if user_element.find('ElevationMaskMaximum') is not None:
                    mask_max_values = user_element.find('ElevationMaskMaximum').text.split(',')
                    user.el_mask_max = [radians(float(n)) for n in mask_max_values]
                else:
                    user.el_mask_max = len(mask_values) * [radians(90.0)]

                user.tle_file_name = misc_fn.resolve_path(user_element.find('TLEFileName').text)
                with open(user.tle_file_name) as file:
                    cnt = 0
                    for line in file:
                        if cnt == 0:
                            user.UserName = line
                        if cnt == 1:
                            user.user_id = int(line[2:7])  # NORAD Catalog Number
                            full_year = "20" + line[18:20]
                            epoch_yr = int(full_year)
                            epoch_doy = float(line[20:31])
                            user.kepler.epoch_mjd = misc_fn.yyyy_doy2mjd(epoch_yr, epoch_doy)
                        if cnt == 2:
                            mu = GM_EARTH
                            mean_motion = float(line[52:59]) / 86400 * 2 * PI  # rad / s
                            user.kepler.semi_major_axis = pow((mu / (pow(mean_motion, 2))), (1.0 / 3.0))
                            eccentricity = "0." + line[26:33]
                            user.kepler.eccentricity = float(eccentricity)
                            user.kepler.inclination = radians(float(line[8:15]))
                            user.kepler.right_ascension = radians(float(line[17:24]))
                            user.kepler.arg_perigee = radians(float(line[34:41]))
                            user.kepler.mean_anomaly = radians(float(line[43:50]))
                            ls.logger.info(['Found user satellite:', str(user.user_id), 'Kepler Elements', str(user.kepler.__dict__)])
                            cnt = -1
                        cnt += 1
                gmst_requested = misc_fn.mjd2gmst(user.kepler.epoch_mjd)
                user.det_posvel_tle(gmst_requested, user.kepler.epoch_mjd)
                user.det_posvel_ecf()
                self.users.append(user)
        ls.logger.info(user.__dict__)
        self.num_user = len(self.users)

    def setup_links(self):

        # Setup the list of objects properly indexed in 2d, eg. [0][0]
        # It is important to have new instance every time, python does by reference!!!
        # STATION <->SATELLITE LINKS
        self.gr2sp = [[Ground2SpaceLink() for j in range(self.num_sat)] for i in range(self.num_station)]
        for idx_station in range(self.num_station):
            for idx_sat in range(self.num_sat):
                if self.stations[idx_station].rx_constellation[self.satellites[idx_sat].constellation_id - 1] == '1':
                    self.gr2sp[idx_station][idx_sat].link_in_use = True
                else:
                    self.gr2sp[idx_station][idx_sat].link_in_use = False
        ls.logger.info(f'Loaded {len(self.gr2sp)*self.num_sat} number of Station to Spacecraft links')

        # USER <->SATELLITE LINKS
        self.usr2sp = [[User2SpaceLink() for j in range(self.num_sat)] for i in range(self.num_user)]
        for idx_user in range(self.num_user):
            for idx_sat in range(self.num_sat):
                if self.users[idx_user].rx_constellation[self.satellites[idx_sat].constellation_id - 1] == '1':
                    self.usr2sp[idx_user][idx_sat].link_in_use = True
                else:
                    self.usr2sp[idx_user][idx_sat].link_in_use = False
        ls.logger.info(f'Loaded {len(self.usr2sp)*self.num_sat} number of User to Spacecraft links')

        # SATELLITE <->SATELLITE LINKS
        self.sp2sp = [[Space2SpaceLink() for j in range(self.num_sat)] for i in range(self.num_sat)]
        for idx_sat1 in range(self.num_sat):
            for idx_sat2 in range(self.num_sat):
                if self.satellites[idx_sat1].rx_constellation[self.satellites[idx_sat2].constellation_id - 1] == '1'\
                        and idx_sat1 != idx_sat2:  # avoid link to itself
                    self.sp2sp[idx_sat1][idx_sat2].link_in_use = True
                else:
                    self.sp2sp[idx_sat1][idx_sat2].link_in_use = False
        ls.logger.info(f'Loaded {len(self.sp2sp)*self.num_sat} number of Spacecraft to Spacecraft links')

    def load_simulation(self):
        # Load the simulation parameters
        tree = ET.parse(self.file_name)
        root = tree.getroot()
        for sim in root.iter('SimulationManager'):
            self.start_time = Time(sim.find('StartDate').text, scale='utc').mjd
            self.time_mjd = self.start_time  # initiate the time loop
            self.stop_time = Time(sim.find('StopDate').text, scale='utc').mjd
            self.time_step = float(sim.find('TimeStep').text)
            self.num_epoch = ceil(86400.0 * (self.stop_time - self.start_time) / self.time_step)

            self.include_gr2sp = misc_fn.str2bool(sim.find('IncludeStation2SpaceLinks').text)
            self.include_usr2sp = misc_fn.str2bool(sim.find('IncludeUser2SpaceLinks').text)
            self.include_sp2sp = misc_fn.str2bool(sim.find('IncludeSpace2SpaceLinks').text)
            self.orbits_from_previous_run = misc_fn.str2bool(sim.find('OrbitsFromPreviousRun').text)
            self.orbit_propagator = sim.find('OrbitPropagator').text
            if sim.find('Report') is not None:  # report.html at the end of the run
                self.report = misc_fn.str2bool(sim.find('Report').text)

            if self.orbit_propagator == 'HPOP':
                # Deferred import: orekit/JVM dependencies only load for HPOP runs
                from propagation_hpop import HpopConfig
                self.hpop_config = HpopConfig()
                if sim.find('HPOP') is not None:
                    self.hpop_config.read_config(sim.find('HPOP'))
                else:
                    ls.logger.warning('OrbitPropagator is HPOP but no <HPOP> block was found; '
                                      'using the default full force model')

            ls.logger.info(f'Loaded simulation, start MJD: {self.start_time}, stop MJD: {self.stop_time},' +
                           f' size time steps in sec: {self.time_step}')
            type_counts = {}  # Occurrences per analysis type, to keep output file names unique
            for analysis_node in root.iter('Analysis'):  # One or more analyses can be performed per run
                type_str = analysis_node.find('Type').text
                analysis = None
                if type_str == 'cov_ground_track':
                    analysis = AnalysisCovGroundTrack()
                if type_str == 'cov_depth_of_coverage':
                    analysis = AnalysisCovDepthOfCoverage()
                if type_str == 'cov_pass_time':
                    analysis = AnalysisCovPassTime()
                if type_str == 'cov_satellite_contour':
                    analysis = AnalysisCovSatelliteContour()
                if type_str == 'cov_satellite_highest':
                    analysis = AnalysisCovSatelliteHighest()
                if type_str == 'cov_satellite_pvt':
                    analysis = AnalysisCovSatellitePvt()
                if type_str == 'cov_satellite_sky_angles':
                    analysis = AnalysisCovSatelliteSkyAngles()
                if type_str == 'cov_satellite_visible':
                    analysis = AnalysisCovSatelliteVisible()
                if type_str == 'cov_satellite_visible_grid':
                    analysis = AnalysisCovSatelliteVisibleGrid()
                if type_str == 'cov_satellite_visible_id':
                    analysis = AnalysisCovSatelliteVisibleId()
                if type_str == 'obs_swath_conical':
                    analysis = AnalysisObsSwathConical()
                if type_str == 'obs_swath_push_broom':
                    analysis = AnalysisObsSwathPushBroom()
                if type_str == 'obs_sza_push_broom':
                    analysis = AnalysisObsSzaPushBroom()
                if type_str == 'obs_sza_subsat':
                    analysis = AnalysisObsSzaSubSat()
                if type_str == 'obs_aoi_revisit':
                    analysis = AnalysisObsAoiRevisit()
                if type_str == 'obs_target_imaging':
                    analysis = AnalysisObsTargetImaging()
                if type_str == 'com_gr2sp_budget':
                    analysis = AnalysisComGr2SpBudget()
                if type_str == 'com_sp2sp_budget':
                    analysis = AnalysisComSp2SpBudget()
                if type_str == 'com_doppler':
                    analysis = AnalysisComDoppler()
                if type_str == 'com_gr2sp_budget_interference':
                    analysis = AnalysisComGr2SpBudgetInterference()
                if type_str == 'com_contact_plan':
                    analysis = AnalysisComContactPlan()
                if type_str == 'com_pfd':
                    analysis = AnalysisComPfd()
                if type_str == 'nav_dilution_of_precision':
                    analysis = AnalysisNavDOP()
                if type_str == 'nav_accuracy':
                    analysis = AnalysisNavAccuracy()
                if type_str == 'sat_battery_depth_discharge':
                    analysis = AnalysisSatBatteryDepthDischarge()
                if type_str == 'sat_eclipse_duration':
                    analysis = AnalysisSatEclipseDuration()
                if type_str == 'sat_data_storage':
                    analysis = AnalysisSatDataStorage()
                if type_str == 'sat_data_latency':
                    analysis = AnalysisSatDataLatency()
                if type_str == 'orb_kepler_elements':
                    analysis = AnalysisOrbKeplerElements()
                if type_str == 'orb_air_density':
                    analysis = AnalysisOrbAirDensity()
                if type_str == 'orb_disturbance_forces':
                    analysis = AnalysisOrbDisturbanceForces()
                if type_str == 'orb_pole_wobble':
                    analysis = AnalysisOrbPoleWobble()
                if type_str == 'orb_deltav_element':
                    analysis = AnalysisOrbDeltaVElement()
                if type_str == 'orb_beta_angle':
                    analysis = AnalysisOrbBetaAngle()
                if type_str == 'orb_lifetime':
                    analysis = AnalysisOrbLifetime()
                if type_str == 'orb_environment':
                    analysis = AnalysisOrbEnvironment()
                if type_str == 'orb_deltav_injection':
                    analysis = AnalysisOrbDeltaVInjection()
                if type_str == 'orb_deltav_reentry':
                    analysis = AnalysisOrbDeltaVReentry()
                if type_str == 'orb_deltav_collision':
                    analysis = AnalysisOrbDeltaVCollision()
                if type_str == 'orb_collision_check':
                    analysis = AnalysisOrbCollisionCheck()
                if type_str == 'orb_collision_alt_check':
                    analysis = AnalysisOrbCollisionAltCheck()
                if type_str == 'sat_thermal':
                    analysis = AnalysisSatThermal()
                if type_str == 'sat_aocs':
                    analysis = AnalysisSatAocs()
                if type_str == 'sat_drag_coefficient':
                    analysis = AnalysisSatDragCoefficient()
                if analysis is None:
                    ls.logger.error(f'Unknown analysis type: {type_str}. '
                                    f'See readme.md for the available analyses.')
                    exit()
                # Output files are named after the analysis type; number repeats so a
                # second analysis of the same type writes e.g. cov_ground_track_2.png
                type_counts[type_str] = type_counts.get(type_str, 0) + 1
                analysis.type = type_str if type_counts[type_str] == 1 else \
                    f'{type_str}_{type_counts[type_str]}'
                analysis.read_config(analysis_node)  # Read the configuration for the specific analysis
                analysis.read_config_map2d(analysis_node)  # Shared 2D-map decoration flags
                self.analyses.append(analysis)
            if not self.analyses:
                ls.logger.warning('No <Analysis> block found; the run will only propagate the orbits')

