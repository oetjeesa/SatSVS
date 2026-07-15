# Import python modules
import argparse
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

# Import project modules
import config
import config_checks
import logging_svs as ls
import misc_fn
from misc_fn import benchmark


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='SatSVS - Satellite Service Volume Simulator. Runs the analyses '
                    'defined in the configuration file (see readme.md).')
    parser.add_argument('config', nargs='?', default='../input/Config.xml',
                        help='Path to the Config.xml scenario file '
                             '(default: ../input/Config.xml)')
    parser.add_argument('-o', '--output-dir', default='../output',
                        help='Directory for plots, data dumps and main.log '
                             '(created if missing, default: ../output)')
    return parser.parse_args(argv)


def load_configuration(config_file='../input/Config.xml', output_dir='../output'):
    config_checks.validate_config(config_file)  # Exits with clear errors on a bad config
    sm = config.AppConfig(config_file, output_dir)
    sm.load_satellites()
    sm.load_stations()
    sm.load_users()
    sm.load_simulation()
    sm.setup_links()
    if sm.orbit_propagator == 'HPOP' and not sm.orbits_from_previous_run:
        import propagation_hpop  # deferred: orekit/JVM only needed for HPOP runs
        sm.hpop = propagation_hpop.HpopPropagation(sm)
    return sm  # Configuration is used as state machine


def run_before_time_loop(sm):
    for analysis in sm.analyses:
        analysis.before_loop(sm)  # Run analyses which are needed before time loop
        analysis.before_loop_map2d(sm)  # Shared 2D-map decorations (ground track memory)
    ls.logger.info('Finished reading configuration file')


def run_time_loop(sm):
    clear_load_orbit_file(sm)  # Clear or load previous orbit file
    precompute_times(sm)  # Precompute per-epoch time conversions in one go
    sm.sat_pos_ecf = np.zeros((sm.num_sat, 3))  # Per-epoch satellite ECF positions for batched link geometry

    for sm.cnt_epoch in range(sm.num_epoch):  # Loop over simulation time window

        convert_times(sm)  # Look up times in gmst, fdoy and string format for this epoch
        if sm.cnt_epoch % 100 == 0:  # Throttled: a log line per epoch is measurable at large num_epoch
            ls.logger.info(f'Simulation time: {sm.time_str}, time step: {sm.cnt_epoch}')

        update_satellites(sm)  # Update pvt on satellites and links
        update_stations(sm)  # Update pvt on ground stations and links
        update_users(sm)  # Update pvt on users and links

        for analysis in sm.analyses:
            analysis.in_loop(sm)  # Run analyses which are needed in time loop
            analysis.in_loop_map2d(sm)  # Shared 2D-map decorations (ground track)

    if sm.file_orbits is not None:  # Close the orbit cache file
        sm.file_orbits.close()
        sm.file_orbits = None


def run_after_time_loop(sm):
    for analysis in sm.analyses:
        ls.logger.info(f'Plotting analysis {analysis.type}')
        analysis.after_loop(sm)  # Run analysis after time loop
        plt.close('all')  # Fresh matplotlib state so analyses cannot touch each other's figures


def precompute_times(sm):
    # All per-epoch time conversions done once, vectorised. Creating an astropy Time
    # plus strptime every epoch used to cost ~1 ms per step; convert_times now only indexes.
    steps = np.arange(sm.num_epoch)
    sm.times_mjd_pre = sm.time_mjd + steps * (sm.time_step / 86400.0)
    sm.times_gmst_pre = np.array([misc_fn.mjd2gmst(mjd) for mjd in sm.times_mjd_pre])
    sm.times_jd_pre = np.floor(sm.times_mjd_pre) + 2400000.5  # Split julian date for sgp4
    sm.times_fr_pre = sm.times_mjd_pre - np.floor(sm.times_mjd_pre)
    iso_strings = Time(sm.times_mjd_pre, format='mjd').iso  # One vectorised astropy call
    sm.times_str_pre = [s[:-4] for s in iso_strings]
    sm.times_datetime_pre = [datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S') for s in sm.times_str_pre]
    sm.times_f_doy_pre = [d.timetuple().tm_yday + d.hour / 24 + d.minute / 60 / 24 + d.second / 3600 / 24
                          for d in sm.times_datetime_pre]


def convert_times(sm):
    cnt = sm.cnt_epoch
    sm.time_mjd = float(sm.times_mjd_pre[cnt])
    sm.time_gmst = float(sm.times_gmst_pre[cnt])
    sm.time_jd = float(sm.times_jd_pre[cnt])
    sm.time_fr = float(sm.times_fr_pre[cnt])
    sm.time_str = sm.times_str_pre[cnt]
    sm.time_datetime = sm.times_datetime_pre[cnt]
    for analysis in sm.analyses:  # Each analysis keeps its own time lists for plotting
        analysis.times_mjd.append(sm.time_mjd)
        analysis.times_f_doy.append(sm.times_f_doy_pre[cnt])


def update_satellites(sm):
    # Compute satellite positions in ECF/ECI and compute the links
    # and remember which ones are in view
    if sm.orbits_from_previous_run:  # Read the ECI posvel from all satellites for one epoch from file
        for idx_sat, satellite in enumerate(sm.satellites):
            satellite.pos_eci = sm.data_orbits[sm.cnt_epoch * sm.num_sat + idx_sat, 0:3]
            satellite.vel_eci = sm.data_orbits[sm.cnt_epoch * sm.num_sat + idx_sat, 3:6]
            satellite.det_posvel_ecf(sm.time_gmst)
            satellite.idx_stat_in_view = []  # Reset before loop
    else:
        for idx_sat, satellite in enumerate(sm.satellites):
            if sm.orbit_propagator == 'Keplerian':
                satellite.det_posvel_eci_keplerian(sm.time_mjd)
            if sm.orbit_propagator == 'SGP4':
                satellite.det_posvel_eci_sgp4(sm.time_jd, sm.time_fr)
            if sm.orbit_propagator == 'HPOP':
                sm.hpop.update_satellite(satellite, idx_sat, sm.time_mjd, sm.time_gmst)
            satellite.det_posvel_ecf(sm.time_gmst)
            satellite.idx_stat_in_view = []  # Reset before loop
        write_posvel_satellites(sm)

    for idx_sat, satellite in enumerate(sm.satellites):  # Gather for the batched link kernels
        sm.sat_pos_ecf[idx_sat, :] = satellite.pos_ecf

    # Compute satellite to satellite links
    if sm.include_sp2sp:  # Only when links needed
        for idx_sat, satellite in enumerate(sm.satellites):
            satellite.idx_sat_in_view = []  # Reset before loop
            # az/el/dist for all satellites from this one in a single numba call
            az, el, az2, el2, dist, los = misc_fn.calc_az_el_dist_batch(satellite.pos_ecf, sm.sat_pos_ecf, True)
            sp2sp_row = sm.sp2sp[idx_sat]
            for idx_sat2, satellite2 in enumerate(sm.satellites):
                link = sp2sp_row[idx_sat2]
                if link.link_in_use:
                    # Same convention as Space2SpaceLink.compute_link: azimuth/elevation
                    # are satellite 1 as seen from satellite 2
                    link.azimuth = az2[idx_sat2]
                    link.elevation = el2[idx_sat2]
                    link.azimuth2 = az[idx_sat2]
                    link.elevation2 = el[idx_sat2]
                    link.distance = dist[idx_sat2]
                    link.sp2sp_ecf = los[idx_sat2]
                    # Check if above elevation mask and not through Earth
                    if link.check_masking(satellite, satellite2):
                        satellite.idx_sat_in_view.append(idx_sat2)


def update_stations(sm):
    # Compute ground station positions/velocities in ECI, compute connection to satellite,
    # and remember which ones are in view
    if sm.include_gr2sp:  # Only when links needed
        for idx_station, station in enumerate(sm.stations):
            station.idx_sat_in_view = []  # Reset before loop
            station.det_posvel_eci(sm.time_gmst)

            # Station to satellite links for all satellites in a single numba call
            az, el, az2, el2, dist, los = misc_fn.calc_az_el_dist_batch(station.pos_ecf, sm.sat_pos_ecf, True)
            gr2sp_row = sm.gr2sp[idx_station]
            for idx_sat, satellite in enumerate(sm.satellites):
                link = gr2sp_row[idx_sat]
                if link.link_in_use:  # from receiver constellation
                    link.azimuth = az[idx_sat]
                    link.elevation = el[idx_sat]
                    link.azimuth2 = az2[idx_sat]
                    link.elevation2 = el2[idx_sat]
                    link.distance = dist[idx_sat]
                    link.gr2sp_ecf = los[idx_sat]
                    if link.check_masking(station):  # Above elevation mask
                        station.idx_sat_in_view.append(idx_sat)
                        # Compute which stations are in view from this satellite (DOC)
                        satellite.idx_stat_in_view.append(idx_station)


def update_users(sm):
    # Compute user positions/velocities in ECI, compute connection to satellite,
    # and remember which ones are in view
    if sm.include_usr2sp:  # Only when links needed
        for idx_user, user in enumerate(sm.users):
            user.idx_sat_in_view = []  # Reset before loop
            if user.type == "Static" or user.type == "Grid":
                user.det_posvel_eci(sm.time_gmst)  # Compute position/velocity in ECI
            if user.type == "Spacecraft":
                user.det_posvel_tle(sm.time_gmst, sm.time_mjd)  # Spacecraft position from TLE

            # User to satellite links for all satellites in a single numba call
            az, el, az2, el2, dist, los = misc_fn.calc_az_el_dist_batch(user.pos_ecf, sm.sat_pos_ecf, False)
            usr2sp_row = sm.usr2sp[idx_user]
            for idx_sat, satellite in enumerate(sm.satellites):
                link = usr2sp_row[idx_sat]
                if link.link_in_use:  # From Receiver Constellation of User
                    link.azimuth = az[idx_sat]
                    link.elevation = el[idx_sat]
                    link.distance = dist[idx_sat]
                    link.usr2sp_ecf = los[idx_sat]
                    if link.check_masking(user):  # Above elevation mask
                        user.idx_sat_in_view.append(idx_sat)


def write_posvel_satellites(sm):  # Cache the propagated ECI orbits for OrbitsFromPreviousRun
    f = sm.file_orbits
    for satellite in sm.satellites:
        f.write("%13.6f,%13.6f,%13.6f,%13.6f,%13.6f,%13.6f\n"
                % (satellite.pos_eci[0], satellite.pos_eci[1], satellite.pos_eci[2],
                   satellite.vel_eci[0], satellite.vel_eci[1], satellite.vel_eci[2]))


def clear_load_orbit_file(sm):  # Clear or load previous orbit file
    sm.file_orbits = None
    if sm.orbits_from_previous_run:
        sm.data_orbits = np.genfromtxt(sm.output_path('orbits_internal.txt'), delimiter=',')
    else:
        # One file handle for the whole run; opening per satellite per epoch dominated the loop
        sm.file_orbits = open(sm.output_path('orbits_internal.txt'), 'w')


if __name__ == '__main__':

    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    ls.init_file(os.path.join(args.output_dir, 'main.log'))

    sm = load_configuration(args.config, args.output_dir)  # Load config into sm status machine holds status of sat, station, user and links

    run_before_time_loop(sm)  # Run the before and during time loop

    run_time_loop(sm)  # Run the before and during time loop

    run_after_time_loop(sm)  # Run the after time loop

    if sm.report:  # <Report>True</Report>: one HTML page with all the results
        import report
        ls.logger.info(f'Written mission report {report.write_report_from_sm(sm)}')

# TODO analysis SZA_pushbroom faster...
# TODO incorporate datashader or geoviews
