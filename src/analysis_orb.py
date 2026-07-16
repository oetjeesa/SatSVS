from math import cos, degrees, radians, sin

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from astropy.coordinates import get_sun
from astropy.time import Time
import astropy.units as u

# Project modules
from constants import GM_EARTH, OMEGA_EARTH, R_EARTH
from analysis import AnalysisBase, make_map_cyl
from analysis_sat import air_density_exponential
import misc_fn
import logging_svs as ls


def _running_mean(x, win):
    """Centred moving average with proper handling of the partial windows at
    the edges (used to turn osculating into per-orbit mean elements)."""
    if win <= 1:
        return np.asarray(x, dtype=float).copy()
    kernel = np.ones(win) / win
    summed = np.convolve(x, kernel, mode='same')
    norm = np.convolve(np.ones(len(x)), kernel, mode='same')
    return summed / norm


def _require_hpop(sm, analysis_type):
    """The orb_ environment analyses sample the Orekit models behind the HPOP
    propagator (atmosphere, force models, EOP data)."""
    if sm.hpop is None:
        ls.logger.error(f'Analysis {analysis_type} needs <OrbitPropagator>HPOP</OrbitPropagator> '
                        f'(and OrbitsFromPreviousRun False). Analysis skipped.')
        return False
    return True


def rv2kepler(pos, vel):
    """Osculating Kepler elements from ECI state vector histories.

    :param pos: (n, 3) positions [m]
    :param vel: (n, 3) velocities [m/s]
    :return: dict of (n,) arrays: sma [m], ecc [-], incl, raan, arg_perigee,
             true_anomaly, mean_anomaly [rad]
    Standard conventions (e.g. Vallado); for near-circular orbits the argument
    of perigee is returned as 0 with the anomaly measured from the ascending
    node, for near-equatorial orbits the RAAN is returned as 0.
    """
    eps = 1e-11
    two_pi = 2.0 * np.pi

    r = np.linalg.norm(pos, axis=1)
    v2 = np.einsum('ij,ij->i', vel, vel)
    h_vec = np.cross(pos, vel)
    n_vec = np.column_stack((-h_vec[:, 1], h_vec[:, 0], np.zeros(len(h_vec))))  # z x h
    n = np.linalg.norm(n_vec, axis=1)
    e_vec = (np.cross(vel, h_vec) - GM_EARTH * pos / r[:, None]) / GM_EARTH
    ecc = np.linalg.norm(e_vec, axis=1)

    sma = 1.0 / (2.0 / r - v2 / GM_EARTH)
    h_norm = np.linalg.norm(h_vec, axis=1)
    incl = np.arccos(np.clip(h_vec[:, 2] / np.where(h_norm < eps, 1.0, h_norm), -1.0, 1.0))
    raan = np.where(n < eps, 0.0, np.arctan2(n_vec[:, 1], n_vec[:, 0]) % two_pi)

    n_safe = np.where(n < eps, 1.0, n)
    e_safe = np.where(ecc < eps, 1.0, ecc)
    arg_perigee = np.arccos(np.clip(np.einsum('ij,ij->i', n_vec, e_vec) /
                                    (n_safe * e_safe), -1.0, 1.0))
    arg_perigee = np.where(e_vec[:, 2] < 0.0, two_pi - arg_perigee, arg_perigee)
    arg_perigee = np.where((n < eps) | (ecc < eps), 0.0, arg_perigee)

    true_anomaly = np.arccos(np.clip(np.einsum('ij,ij->i', e_vec, pos) /
                                     (e_safe * r), -1.0, 1.0))
    true_anomaly = np.where(np.einsum('ij,ij->i', pos, vel) < 0.0,
                            two_pi - true_anomaly, true_anomaly)
    # Near-circular: anomaly measured from the ascending node (argument of latitude)
    arg_lat = np.arccos(np.clip(np.einsum('ij,ij->i', n_vec, pos) / (n_safe * r), -1.0, 1.0))
    arg_lat = np.where(pos[:, 2] < 0.0, two_pi - arg_lat, arg_lat)
    true_anomaly = np.where(ecc < eps, arg_lat, true_anomaly)

    # Eccentric and mean anomaly (elliptical orbits)
    ecc_ell = np.minimum(ecc, 1.0 - 1e-12)
    ecc_anomaly = 2.0 * np.arctan2(np.sqrt(1.0 - ecc_ell) * np.sin(true_anomaly / 2.0),
                                   np.sqrt(1.0 + ecc_ell) * np.cos(true_anomaly / 2.0))
    mean_anomaly = (ecc_anomaly - ecc_ell * np.sin(ecc_anomaly)) % two_pi

    return {'sma': sma, 'ecc': ecc, 'incl': incl, 'raan': raan,
            'arg_perigee': arg_perigee, 'true_anomaly': true_anomaly,
            'mean_anomaly': mean_anomaly}


class AnalysisOrbKeplerElements(AnalysisBase):
    """Evolution of the osculating Kepler elements over the simulation time,
    computed each epoch from the ECI state vector: semi-major axis,
    eccentricity, inclination, RAAN, argument of perigee and mean anomaly.
    With the HPOP propagator this shows the perturbation effects (drag decay of
    the semi-major axis, J2 RAAN drift, etc.); it also works with the other
    propagators (constant elements for Keplerian, mean-element variations for
    SGP4)."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection
        self.sat_metric = None  # Per-satellite metric memory (num_sat, num_epoch, 6)

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)

    def _selected(self, satellite):
        if self.constellation_id > 0 and satellite.constellation_id != self.constellation_id:
            return False
        if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:
            return False
        return True

    def before_loop(self, sm):
        self.sat_metric = np.full((sm.num_sat, sm.num_epoch, 6), np.nan)  # ECI pos + vel

    def in_loop(self, sm):
        for idx_sat, satellite in enumerate(sm.satellites):
            if self._selected(satellite):
                self.sat_metric[idx_sat, sm.cnt_epoch, 0:3] = satellite.pos_eci
                vel = np.asarray(satellite.vel_eci, dtype=float)
                if sm.orbit_propagator == 'HPOP':
                    # The HPOP tool-frame velocity is the Earth-relative (ITRF)
                    # velocity (see propagation_hpop docstring); add the omega x r
                    # transport term to recover the inertial velocity the
                    # osculating element computation needs
                    pos = satellite.pos_eci
                    vel = vel + OMEGA_EARTH * np.array([-pos[1], pos[0], 0.0])
                self.sat_metric[idx_sat, sm.cnt_epoch, 3:6] = vel

    def after_loop(self, sm):
        panels = [  # (key, scale factor, label, file name suffix)
            ('sma', 1e-3, 'Semi-major axis [km]', 'semi_major_axis'),
            ('ecc', 1.0, 'Eccentricity [-]', 'eccentricity'),
            ('incl', np.degrees(1.0), 'Inclination [deg]', 'inclination'),
            ('raan', np.degrees(1.0), 'RAAN [deg]', 'raan'),
            ('arg_perigee', np.degrees(1.0), 'Argument of perigee [deg]', 'arg_perigee'),
            ('mean_anomaly', np.degrees(1.0), 'Mean anomaly [deg]', 'mean_anomaly'),
        ]
        times = np.asarray(self.times_f_doy)
        results = []  # (satellite, used epochs, osculating elements)
        csv_rows = []
        for idx_sat, satellite in enumerate(sm.satellites):
            if not self._selected(satellite):
                continue
            used = ~np.isnan(self.sat_metric[idx_sat, :, 0])
            if not used.any():
                continue
            elements = rv2kepler(self.sat_metric[idx_sat, used, 0:3], self.sat_metric[idx_sat, used, 3:6])
            csv_rows.append(np.column_stack(
                [times[used], np.full(used.sum(), satellite.sat_id), elements['sma'],
                 elements['ecc'], np.degrees(elements['incl']), np.degrees(elements['raan']),
                 np.degrees(elements['arg_perigee']), np.degrees(elements['mean_anomaly'])]))
            sma_km = elements['sma'] / 1000.0
            # Secular change: difference of the mean over the first and last
            # ~100 epochs, which averages out the J2 short-period oscillation
            n_avg = min(100, len(sma_km))
            first, last = np.mean(sma_km[:n_avg]), np.mean(sma_km[-n_avg:])
            ls.logger.info(f'Satellite {satellite.sat_id}: mean SMA first epochs {first:.3f} km, '
                           f'last epochs {last:.3f} km, change {(last-first)*1000:.1f} m')
            results.append((satellite, used, elements))
        if not results:
            ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                            f'SatelliteID {self.satellite_id}. No plot produced.')
            return
        # One plot per element (<type>_semi_major_axis.png, ...) instead of a
        # single six-panel figure
        for key, scale, label, suffix in panels:
            fig, ax = plt.subplots(figsize=(10, 6))
            for satellite, used, elements in results:
                ax.plot(times[used], elements[key] * scale, '-', linewidth=0.9,
                        label=f'Sat {satellite.sat_id}')
            ax.set_xlabel('DOY [-]')
            ax.set_ylabel(label)
            ax.set_title(f'Osculating Kepler elements: {label}')
            ax.grid(True)
            ax.legend(fontsize=8)
            fig.tight_layout()
            plt.savefig(sm.output_path(f'{self.type}_{suffix}.png'))
            plt.show()

        self.write_csv(sm, ['doy', 'sat_id', 'sma_m', 'eccentricity', 'inclination_deg',
                            'raan_deg', 'arg_perigee_deg', 'mean_anomaly_deg'],
                       np.vstack(csv_rows))


class AnalysisOrbAirDensity(AnalysisBase):
    """Atmospheric density at the satellite altitude over the simulation time,
    sampled per epoch from the atmosphere model configured as HPOP DragModel
    (NRLMSISE00, DTM2000 or HarrisPriester), together with the altitude
    itself. Needs the HPOP propagator."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection
        self.sat_metric = None  # (num_sat, num_epoch, 2): altitude [m], density [kg/m3]
        self.enabled = True

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)

    def _selected(self, satellite):
        if self.constellation_id > 0 and satellite.constellation_id != self.constellation_id:
            return False
        if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:
            return False
        return True

    def before_loop(self, sm):
        self.enabled = _require_hpop(sm, self.type)
        self.sat_metric = np.full((sm.num_sat, sm.num_epoch, 2), np.nan)

    def in_loop(self, sm):
        if not self.enabled:
            return
        for idx_sat, satellite in enumerate(sm.satellites):
            if self._selected(satellite):
                satellite.det_lla()
                r_earth = misc_fn.earth_radius_lat(satellite.lla[0])
                self.sat_metric[idx_sat, sm.cnt_epoch, 0] = \
                    np.linalg.norm(satellite.pos_ecf) - r_earth
                self.sat_metric[idx_sat, sm.cnt_epoch, 1] = \
                    sm.hpop.air_density(satellite.pos_ecf, sm.time_mjd)

    def after_loop(self, sm):
        if not self.enabled:
            return
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax2 = ax1.twinx()
        times = np.asarray(self.times_f_doy)
        plotted = False
        for idx_sat, satellite in enumerate(sm.satellites):
            if not self._selected(satellite):
                continue
            used = ~np.isnan(self.sat_metric[idx_sat, :, 1])
            if not used.any():
                continue
            ax1.semilogy(times[used], self.sat_metric[idx_sat, used, 1], '-',
                         linewidth=0.9, label=f'Density Sat {satellite.sat_id}')
            ax2.plot(times[used], self.sat_metric[idx_sat, used, 0] / 1000.0, '--',
                     linewidth=0.9, label=f'Altitude Sat {satellite.sat_id}')
            plotted = True
        if not plotted:
            ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                            f'SatelliteID {self.satellite_id}. No plot produced.')
            return
        ax1.set_xlabel('DOY [-]')
        ax1.set_ylabel('Air density [kg/m$^3$]')
        ax2.set_ylabel('Altitude [km]')
        ax1.grid(True)
        ax1.legend(loc='upper left', fontsize=8)
        ax2.legend(loc='upper right', fontsize=8)
        fig.suptitle(f'Atmospheric density at satellite altitude ({sm.hpop.cfg.drag_model})')
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        times = np.asarray(self.times_f_doy)
        self.write_csv(sm, ['doy', 'sat_id', 'altitude_m', 'density_kgm3'],
                       np.vstack([np.column_stack([times, np.full(sm.num_epoch, satellite.sat_id),
                                                   self.sat_metric[idx_sat]])
                                  for idx_sat, satellite in enumerate(sm.satellites)
                                  if self._selected(satellite)]))


class AnalysisOrbDisturbanceForces(AnalysisBase):
    """Magnitude of every enabled HPOP perturbation acceleration on the first
    selected satellite over the simulation time, evaluated per epoch on the
    propagated state (geopotential harmonics, drag, solar radiation pressure,
    third bodies, tides, relativity), with the central gravity as reference.
    Needs the HPOP propagator."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection
        self.idx_found_satellite = None
        self.forces = None  # (name, ForceModel) pairs
        self.force_names = []
        self.metric = None  # (num_epoch, 1 + num_forces) accelerations [m/s2]
        self.enabled = True

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)

    def _selected(self, satellite):
        if self.constellation_id > 0 and satellite.constellation_id != self.constellation_id:
            return False
        if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:
            return False
        return True

    def before_loop(self, sm):
        self.enabled = _require_hpop(sm, self.type)
        if not self.enabled:
            return
        for idx_sat, satellite in enumerate(sm.satellites):
            if self._selected(satellite):
                self.idx_found_satellite = idx_sat
                break
        if self.idx_found_satellite is None:
            ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                            f'SatelliteID {self.satellite_id}. Analysis skipped.')
            self.enabled = False
            return
        self.forces = sm.hpop.named_forces()
        self.force_names = ['Central gravity'] + [name for name, force in self.forces]
        self.metric = np.full((sm.num_epoch, len(self.force_names)), np.nan)
        ls.logger.info(f'Disturbance forces evaluated for satellite '
                       f'{sm.satellites[self.idx_found_satellite].sat_id}: ' +
                       ', '.join(self.force_names))

    def in_loop(self, sm):
        if not self.enabled:
            return
        state = sm.hpop.sample_state(self.idx_found_satellite, sm.time_mjd)
        r = state.getPVCoordinates().getPosition().getNorm()
        self.metric[sm.cnt_epoch, 0] = sm.hpop.mu / r ** 2
        for i, (name, force) in enumerate(self.forces):
            acc = force.acceleration(state, force.getParameters())
            self.metric[sm.cnt_epoch, i + 1] = acc.getNorm()

    def after_loop(self, sm):
        if not self.enabled:
            return
        times = np.asarray(self.times_f_doy)
        fig = plt.figure(figsize=(10, 6))
        for i, name in enumerate(self.force_names):
            style = 'k--' if i == 0 else '-'
            plt.semilogy(times, self.metric[:, i], style, linewidth=0.9, label=name)
        plt.xlabel('DOY [-]')
        plt.ylabel('Acceleration [m/s$^2$]')
        plt.title(f'Perturbation accelerations, satellite '
                  f'{sm.satellites[self.idx_found_satellite].sat_id}')
        plt.grid(True)
        plt.legend(fontsize=8)
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        columns = ['doy'] + [name.lower().replace(' ', '_') for name in self.force_names]
        self.write_csv(sm, columns, np.column_stack([times, self.metric]))


class AnalysisOrbPoleWobble(AnalysisBase):
    """Wobble of the Earth rotation axis (IERS polar motion xp/yp from the EOP
    data) over the simulation time: time series and the xp/yp trace. The
    Chandler + annual wobble circles ~0.1-0.3 arcsec in roughly a year, so
    longer simulation windows show more of the circle. Needs the HPOP
    propagator (for the Orekit EOP data)."""

    def __init__(self):
        super().__init__()
        self.metric = None  # (num_epoch, 2): xp, yp [arcsec]
        self.enabled = True

    def read_config(self, node):
        pass

    def before_loop(self, sm):
        self.enabled = _require_hpop(sm, self.type)
        self.metric = np.full((sm.num_epoch, 2), np.nan)

    def in_loop(self, sm):
        if not self.enabled:
            return
        xp, yp = sm.hpop.pole_correction(sm.time_mjd)
        self.metric[sm.cnt_epoch] = [degrees(xp) * 3600.0, degrees(yp) * 3600.0]

    def after_loop(self, sm):
        if not self.enabled:
            return
        times = np.asarray(self.times_f_doy)
        # Two separate plots (<type>_timeseries.png and <type>_track.png)
        # instead of a single two-panel figure
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(times, self.metric[:, 0], 'r-', label='xp')
        ax1.plot(times, self.metric[:, 1], 'b-', label='yp')
        ax1.set_xlabel('DOY [-]')
        ax1.set_ylabel('Polar motion [arcsec]')
        ax1.set_title('Earth pole wobble (IERS polar motion)')
        ax1.grid(True)
        ax1.legend()
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '_timeseries.png'))
        plt.show()

        fig, ax2 = plt.subplots(figsize=(8, 6))
        ax2.plot(self.metric[:, 0], self.metric[:, 1], 'g-')
        ax2.plot(self.metric[0, 0], self.metric[0, 1], 'go', label='start')
        ax2.plot(self.metric[-1, 0], self.metric[-1, 1], 'rs', label='end')
        ax2.set_xlabel('xp [arcsec]')
        ax2.set_ylabel('yp [arcsec]')
        ax2.set_title('Earth pole wobble track (IERS polar motion)')
        ax2.set_aspect('equal', adjustable='datalim')
        ax2.grid(True)
        ax2.legend()
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '_track.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'xp_arcsec', 'yp_arcsec'],
                       np.column_stack([times, self.metric]))


class AnalysisOrbDeltaVElement(AnalysisBase):
    """Station-keeping delta-v estimate for one orbit element of one satellite.

    The element (mean, i.e. the osculating value averaged over one orbital
    period) drifts under the modelled perturbations; whenever the controlled
    value leaves the deadband around the target an impulsive correction resets
    it to the target, costed with the standard impulsive-maneuver formulas:
    tangential burn for altitude/semi-major axis/eccentricity, plane change
    for inclination/RAAN, apsidal rotation for the argument of perigee. The
    drift rate is taken from the uncontrolled propagation (valid to first
    order for the small offsets inside a deadband). Works with any propagator;
    with drift-free propagation (Keplerian) or short windows the required
    delta-v is simply zero."""

    KEYS = ('altitude', 'semimajoraxis', 'eccentricity', 'inclination', 'raan',
            'argofperigee')

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection (first match is analysed)
        self.idx_found_satellite = None
        self.target_type = 'Altitude'
        self.target_value = None  # None: the element value at simulation start
        self.dead_band = None  # Half width of the deadband (m / - / rad)
        self.sat_metric = None  # (num_epoch, 6): ECI pos + vel of the satellite
        self.enabled = True

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)
        if node.find('TargetType') is not None:
            self.target_type = node.find('TargetType').text
        angular = self.target_type.lower() in ('inclination', 'raan', 'argofperigee')
        if node.find('TargetValue') is not None:
            self.target_value = float(node.find('TargetValue').text)
            if angular:
                self.target_value = radians(self.target_value)
        if node.find('DeadBand') is not None:
            self.dead_band = float(node.find('DeadBand').text)
            if angular:
                self.dead_band = radians(self.dead_band)

    def _selected(self, satellite):
        if self.constellation_id > 0 and satellite.constellation_id != self.constellation_id:
            return False
        if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:
            return False
        return True

    def before_loop(self, sm):
        if self.target_type.lower() not in self.KEYS:
            ls.logger.error(f'Unknown TargetType {self.target_type} (use Altitude, '
                            f'SemiMajorAxis, Eccentricity, Inclination, RAAN or '
                            f'ArgOfPerigee). Analysis skipped.')
            self.enabled = False
            return
        if self.dead_band is None or self.dead_band <= 0.0:
            ls.logger.error(f'Analysis {self.type} needs a positive <DeadBand> '
                            f'(m for Altitude/SemiMajorAxis, deg for the angles, '
                            f'- for Eccentricity). Analysis skipped.')
            self.enabled = False
            return
        for idx_sat, satellite in enumerate(sm.satellites):
            if self._selected(satellite):
                self.idx_found_satellite = idx_sat
                break
        if self.idx_found_satellite is None:
            ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                            f'SatelliteID {self.satellite_id}. Analysis skipped.')
            self.enabled = False
            return
        self.sat_metric = np.full((sm.num_epoch, 6), np.nan)

    def in_loop(self, sm):
        if not self.enabled:
            return
        satellite = sm.satellites[self.idx_found_satellite]
        self.sat_metric[sm.cnt_epoch, 0:3] = satellite.pos_eci
        vel = np.asarray(satellite.vel_eci, dtype=float)
        if sm.orbit_propagator == 'HPOP':
            # Earth-relative tool-frame velocity -> inertial (see orb_kepler_elements)
            pos = satellite.pos_eci
            vel = vel + OMEGA_EARTH * np.array([-pos[1], pos[0], 0.0])
        self.sat_metric[sm.cnt_epoch, 3:6] = vel

    def _dv_for_correction(self, key, delta, a, e, i):
        """Delta-v [m/s] of an impulsive correction of the element by delta,
        with the circular-speed approximation v = sqrt(mu/a)."""
        v = np.sqrt(GM_EARTH / a)
        if key in ('altitude', 'semimajoraxis'):
            return v * delta / (2.0 * a)  # Tangential burn: da = 2 a dv / v
        if key == 'eccentricity':
            return v * delta / 2.0  # Tangential burn: de = 2 dv / v
        if key == 'inclination':
            return 2.0 * v * np.sin(delta / 2.0)  # Plane change at the node
        if key == 'raan':
            theta = delta * np.sin(i)  # Plane separation of a small RAAN change
            return 2.0 * v * np.sin(theta / 2.0)
        # argofperigee: in-plane rotation of the line of apsides
        p = a * max(1.0 - e ** 2, 1e-12)
        return 2.0 * np.sqrt(GM_EARTH / p) * e * np.sin(delta / 2.0)

    def after_loop(self, sm):
        if not self.enabled:
            return
        elements = rv2kepler(self.sat_metric[:, 0:3], self.sat_metric[:, 3:6])
        # Mean elements: average the osculating values over one orbital period,
        # otherwise the J2 short-period oscillation triggers spurious maneuvers
        period = 2.0 * np.pi * np.sqrt(np.nanmean(elements['sma']) ** 3 / GM_EARTH)
        win = min(max(1, int(round(period / sm.time_step))), len(self.sat_metric))
        sma = _running_mean(elements['sma'], win)
        ecc = _running_mean(elements['ecc'], win)
        incl = _running_mean(np.unwrap(elements['incl']), win)
        raan = _running_mean(np.unwrap(elements['raan']), win)
        argp = _running_mean(np.unwrap(elements['arg_perigee']), win)
        # Drop the half-window edges: their partial averages still carry the
        # short-period oscillation and would trigger spurious maneuvers
        half = win // 2
        span = slice(half, len(sma) - half) if len(sma) > 3 * half > 0 else slice(None)
        sma, ecc, incl = sma[span], ecc[span], incl[span]
        raan, argp = raan[span], argp[span]
        times = np.asarray(self.times_f_doy)[span]

        key = self.target_type.lower()
        series = {'altitude': sma - R_EARTH, 'semimajoraxis': sma,
                  'eccentricity': ecc, 'inclination': incl, 'raan': raan,
                  'argofperigee': argp}[key]
        angular = key in ('inclination', 'raan', 'argofperigee')

        target = self.target_value
        if target is None:
            target = series[0]
            ls.logger.info(f'{self.type}: no TargetValue given, using the initial '
                           f'{self.target_type} of the simulation')
        elif angular:  # Bring the target into the branch of the unwrapped series
            target = series[0] + (target - series[0] + np.pi) % (2.0 * np.pi) - np.pi

        # Deadband control emulation on the drifting mean element
        controlled = np.empty(len(series))
        dv_cum = np.zeros(len(series))
        maneuver_epochs = []
        offset, dv_total = 0.0, 0.0
        for k in range(len(series)):
            value = series[k] + offset
            if abs(value - target) > self.dead_band:
                delta = target - value
                dv_total += self._dv_for_correction(key, abs(delta),
                                                    sma[k], ecc[k], incl[k])
                offset += delta
                value = target
                maneuver_epochs.append(k)
            controlled[k] = value
            dv_cum[k] = dv_total

        days = len(series) * sm.time_step / 86400.0
        ls.logger.info(f'{self.type}: {self.target_type} kept within +/- '
                       f'{self.dead_band} of {target:.6g} with '
                       f'{len(maneuver_epochs)} maneuver(s), total delta-v '
                       f'{dv_total:.3f} m/s ({dv_total / days * 365.25:.1f} m/s/year)')

        if key in ('altitude', 'semimajoraxis'):
            unit = lambda x: np.asarray(x) / 1000.0
            ylabel = self.target_type + ' [km]'
        elif angular:
            unit = lambda x: np.degrees(x)
            ylabel = self.target_type + ' [deg]'
        else:
            unit = lambda x: np.asarray(x)
            ylabel = 'Eccentricity [-]'
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(times, unit(series), color='gray', linestyle='--', linewidth=0.9,
                 label='Uncontrolled drift')
        ax1.plot(times, unit(controlled), 'b-', linewidth=0.9, label='Controlled')
        ax1.axhline(unit(target), color='green', linewidth=1.0, label='Target')
        ax1.axhline(unit(target + self.dead_band), color='red', linestyle=':',
                    linewidth=1.0, label='Deadband')
        ax1.axhline(unit(target - self.dead_band), color='red', linestyle=':',
                    linewidth=1.0)
        for cnt, k in enumerate(maneuver_epochs):
            ax1.axvline(times[k], color='orange', linewidth=0.6, alpha=0.6,
                        label='Maneuver' if cnt == 0 else None)
        ax1.set_xlabel('DOY [-]')
        ax1.set_ylabel(ylabel)
        ax1.grid(True)
        ax1.legend(loc='upper left', fontsize=8)
        ax2 = ax1.twinx()  # Cumulative delta-v on the right axis
        ax2.step(times, dv_cum, 'r-', linewidth=1.2, where='post',
                 label='Cumulative delta-v')
        ax2.set_ylabel('Cumulative delta-v [m/s]', color='red')
        ax2.tick_params(axis='y', colors='red')
        ax2.set_ylim(bottom=0)
        ax2.legend(loc='upper right', fontsize=8)
        sat_id = sm.satellites[self.idx_found_satellite].sat_id
        fig.suptitle(f'{self.target_type} station keeping, satellite {sat_id}: '
                     f'{len(maneuver_epochs)} maneuvers, {dv_total:.2f} m/s '
                     f'({dv_total / days * 365.25:.1f} m/s/year)')
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['doy', f'uncontrolled_{key}', f'controlled_{key}',
                            'cumulative_deltav_ms'],
                       np.column_stack([times, series, controlled, dv_cum]))


class AnalysisOrbBetaAngle(AnalysisBase):
    """Solar beta angle (angle between the Sun direction and the orbit plane)
    over the simulation time, together with the analytic eclipse fraction of a
    circular orbit at that beta angle. Beta drives the eclipse pattern, the
    thermal hot/cold cases and the power sizing, so this analysis ties the
    sat_ platform analyses together; run it over months to see the seasonal cycle.
    Works with any propagator - note that Keplerian elements have no J2 nodal
    regression, so for non-sun-synchronous orbits use SGP4 or HPOP (or an
    LTAN-defined SSO orbit) to capture the full beta cycle."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection
        self.sat_metric = None  # (num_sat, num_epoch, 2): beta [rad], eclipse fraction

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)

    def _selected(self, satellite):
        if self.constellation_id > 0 and satellite.constellation_id != self.constellation_id:
            return False
        if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:
            return False
        return True

    def before_loop(self, sm):
        self.sat_metric = np.full((sm.num_sat, sm.num_epoch, 2), np.nan)

    def in_loop(self, sm):
        # One Sun direction per epoch, shared by all satellites
        sun_pos_eci = get_sun(Time(sm.time_mjd, format='mjd')).cartesian.xyz.to(u.m).value
        sun_dir = sun_pos_eci / np.linalg.norm(sun_pos_eci)
        for idx_sat, satellite in enumerate(sm.satellites):
            if not self._selected(satellite):
                continue
            pos = np.asarray(satellite.pos_eci, dtype=float)
            vel = np.asarray(satellite.vel_eci, dtype=float)
            if sm.orbit_propagator == 'HPOP':
                # Earth-relative tool-frame velocity -> inertial (see orb_kepler_elements)
                vel = vel + OMEGA_EARTH * np.array([-pos[1], pos[0], 0.0])
            h_vec = np.cross(pos, vel)
            beta = np.arcsin(np.clip(np.dot(h_vec / np.linalg.norm(h_vec), sun_dir),
                                     -1.0, 1.0))
            # Analytic eclipse fraction of a circular orbit (cylindrical shadow):
            # eclipse exists while |beta| < asin(R/r), with arc half-angle from
            # the shadow cylinder cross section
            r = np.linalg.norm(pos)
            if abs(beta) < np.arcsin(min(R_EARTH / r, 1.0)):
                fraction = np.arccos(np.clip(np.sqrt(r ** 2 - R_EARTH ** 2) /
                                             (r * np.cos(beta)), -1.0, 1.0)) / np.pi
            else:
                fraction = 0.0
            self.sat_metric[idx_sat, sm.cnt_epoch] = [beta, fraction]

    def after_loop(self, sm):
        times = np.asarray(self.times_f_doy)
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax2 = ax1.twinx()
        plotted = False
        cnt_plotted = 0
        csv_rows = []
        for idx_sat, satellite in enumerate(sm.satellites):
            if not self._selected(satellite):
                continue
            used = ~np.isnan(self.sat_metric[idx_sat, :, 0])
            if not used.any():
                continue
            beta_deg = np.degrees(self.sat_metric[idx_sat, used, 0])
            fraction = self.sat_metric[idx_sat, used, 1]
            period = 2.0 * np.pi * np.sqrt(np.linalg.norm(satellite.pos_eci) ** 3 / GM_EARTH)
            ls.logger.info(f'{self.type}: satellite {satellite.sat_id} beta '
                           f'{beta_deg.min():.1f} .. {beta_deg.max():.1f} deg, max eclipse '
                           f'{fraction.max() * period / 60.0:.1f} min/orbit '
                           f'({fraction.max() * 100:.1f}% of the orbit)')
            # Both twin axes start their colour cycle at C0, so pair the
            # colours explicitly: beta and eclipse of the same satellite in
            # contrasting cycle colours (C0/C1, C2/C3, ...)
            ax1.plot(times[used], beta_deg, '-', linewidth=1.0,
                     color=f'C{2 * cnt_plotted % 10}',
                     label=f'Beta Sat {satellite.sat_id}')
            ax2.plot(times[used], fraction * 100.0, '--', linewidth=0.9,
                     color=f'C{(2 * cnt_plotted + 1) % 10}',
                     label=f'Eclipse Sat {satellite.sat_id}')
            cnt_plotted += 1
            csv_rows.append(np.column_stack([times[used],
                                             np.full(used.sum(), satellite.sat_id),
                                             beta_deg, fraction]))
            plotted = True
        if not plotted:
            ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                            f'SatelliteID {self.satellite_id}. No plot produced.')
            return
        ax1.axhline(0.0, color='gray', linewidth=0.6)
        ax1.set_xlabel('DOY [-]')
        # Colour each y-axis like its (first) curve: beta C0, eclipse C1
        ax1.set_ylabel('Solar beta angle [deg]', color='C0')
        ax1.tick_params(axis='y', labelcolor='C0')
        ax2.set_ylabel('Eclipse fraction of the orbit [%]', color='C1')
        ax2.tick_params(axis='y', labelcolor='C1')
        ax2.set_ylim(bottom=0)
        ax1.grid(True)
        ax1.legend(loc='upper left', fontsize=8)
        ax2.legend(loc='upper right', fontsize=8)
        fig.suptitle('Solar beta angle and analytic eclipse fraction')
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'sat_id', 'beta_deg', 'eclipse_fraction'],
                       np.vstack(csv_rows))


class AnalysisOrbLifetime(AnalysisBase):
    """Orbital lifetime under atmospheric drag: the mean semi-major axis at
    the end of the simulation window is decayed semi-analytically
    (da/dt = -rho * Cd*A/m * sqrt(mu*a), circular orbit approximation, with a
    piecewise-exponential atmosphere scaled by DensityScale for solar
    activity) until the re-entry altitude or the MaxYears horizon. Reports
    compliance with the 25-year debris-mitigation rule, the delta-v of an
    immediate deorbit burn (perigee lowered to the re-entry altitude) and,
    when non-compliant, the delta-v to move to a 25-year compliant circular
    orbit. Works with any propagator - with HPOP the projection starts from
    the actually decayed state at the end of the window."""

    REENTRY_ALTITUDE = 120e3  # Re-entry interface altitude [m]
    RULE_YEARS = 25.0  # Debris-mitigation lifetime rule

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection (first match is analysed)
        self.idx_found_satellite = None
        self.mass = None  # [kg]; default constellation Mass
        self.drag_area = None  # [m2]; default constellation FrontalArea
        self.drag_cd = 2.2
        self.density_scale = 1.0  # ~0.5 solar minimum, ~2 solar maximum
        self.max_years = 100.0  # Integration horizon
        self.reentry_altitude = self.REENTRY_ALTITUDE
        self.sat_metric = None  # (num_epoch, 6): ECI pos + vel of the satellite
        self.enabled = True

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)
        if node.find('Mass') is not None:
            self.mass = float(node.find('Mass').text)
        if node.find('DragArea') is not None:
            self.drag_area = float(node.find('DragArea').text)
        if node.find('DragCoefficient') is not None:
            self.drag_cd = float(node.find('DragCoefficient').text)
        if node.find('DensityScale') is not None:
            self.density_scale = float(node.find('DensityScale').text)
        if node.find('MaxYears') is not None:
            self.max_years = float(node.find('MaxYears').text)
        if node.find('ReentryAltitude') is not None:
            self.reentry_altitude = float(node.find('ReentryAltitude').text)

    def _selected(self, satellite):
        if self.constellation_id > 0 and satellite.constellation_id != self.constellation_id:
            return False
        if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:
            return False
        return True

    def before_loop(self, sm):
        for idx_sat, satellite in enumerate(sm.satellites):
            if self._selected(satellite):
                self.idx_found_satellite = idx_sat
                break
        if self.idx_found_satellite is None:
            ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                            f'SatelliteID {self.satellite_id}. Analysis skipped.')
            self.enabled = False
            return
        satellite = sm.satellites[self.idx_found_satellite]
        if self.mass is None:
            self.mass = satellite.mass
        if self.drag_area is None:
            self.drag_area = satellite.frontal_area
        if self.mass is None or self.drag_area is None:
            ls.logger.error(f'{self.type} needs <Mass> and <DragArea> (in the Analysis '
                            f'block or as Mass/FrontalArea of the constellation). '
                            f'Analysis skipped.')
            self.enabled = False
            return
        self.sat_metric = np.full((sm.num_epoch, 6), np.nan)

    def in_loop(self, sm):
        if not self.enabled:
            return
        satellite = sm.satellites[self.idx_found_satellite]
        self.sat_metric[sm.cnt_epoch, 0:3] = satellite.pos_eci
        vel = np.asarray(satellite.vel_eci, dtype=float)
        if sm.orbit_propagator == 'HPOP':
            # Earth-relative tool-frame velocity -> inertial (see orb_kepler_elements)
            pos = satellite.pos_eci
            vel = vel + OMEGA_EARTH * np.array([-pos[1], pos[0], 0.0])
        self.sat_metric[sm.cnt_epoch, 3:6] = vel

    def _decay_profile(self, sma_start):
        """Semi-analytic decay of the circular-orbit radius from sma_start
        down to the re-entry altitude, capped at max_years. Returns times
        [years] and altitudes [m]. Step size adapts to ~500 m of decay per
        step, bounded to [60 s, 10 days]."""
        ballistic = self.drag_cd * self.drag_area / self.mass  # Cd*A/m [m2/kg]
        reentry_r = R_EARTH + self.reentry_altitude
        max_seconds = self.max_years * 365.25 * 86400.0
        a, t = float(sma_start), 0.0
        times, alts = [0.0], [a - R_EARTH]
        while a > reentry_r and t < max_seconds and len(times) < 2_000_000:
            rho = air_density_exponential(a - R_EARTH) * self.density_scale
            da_dt = rho * ballistic * np.sqrt(GM_EARTH * a)
            dt = min(max(500.0 / da_dt if da_dt > 0 else 10 * 86400.0, 60.0), 10 * 86400.0)
            a -= da_dt * dt
            t += dt
            times.append(t / (365.25 * 86400.0))
            alts.append(max(a - R_EARTH, self.reentry_altitude))
        return np.asarray(times), np.asarray(alts)

    def _lifetime_years(self, sma_start):
        times, alts = self._decay_profile(sma_start)
        return times[-1] if alts[-1] <= self.reentry_altitude else np.inf

    def _hohmann_dv(self, a_from, a_to):
        """Delta-v [m/s] of a two-burn Hohmann transfer between circular orbits."""
        dv1 = np.sqrt(GM_EARTH / a_from) * abs(np.sqrt(2 * a_to / (a_from + a_to)) - 1.0)
        dv2 = np.sqrt(GM_EARTH / a_to) * abs(1.0 - np.sqrt(2 * a_from / (a_from + a_to)))
        return dv1 + dv2

    def after_loop(self, sm):
        if not self.enabled:
            return
        # Start from the per-orbit mean semi-major axis at the end of the window
        elements = rv2kepler(self.sat_metric[:, 0:3], self.sat_metric[:, 3:6])
        period = 2.0 * np.pi * np.sqrt(np.nanmean(elements['sma']) ** 3 / GM_EARTH)
        win = min(max(1, int(round(period / sm.time_step))), len(elements['sma']))
        sma_start = float(np.nanmean(elements['sma'][-win:]))

        times, alts = self._decay_profile(sma_start)
        decayed = alts[-1] <= self.reentry_altitude
        lifetime = times[-1] if decayed else np.inf
        compliant = lifetime <= self.RULE_YEARS

        # Immediate deorbit: perigee lowered to the re-entry altitude
        reentry_r = R_EARTH + self.reentry_altitude
        v_c = np.sqrt(GM_EARTH / sma_start)
        dv_deorbit = v_c * (1.0 - np.sqrt(2.0 * reentry_r / (sma_start + reentry_r)))
        lifetime_str = f'{lifetime:.2f} years' if decayed else f'> {self.max_years:.0f} years'
        ls.logger.info(f'{self.type}: start altitude {(sma_start - R_EARTH) / 1000:.1f} km, '
                       f'Cd*A/m {self.drag_cd * self.drag_area / self.mass:.4f} m2/kg, '
                       f'density scale {self.density_scale}: lifetime {lifetime_str} '
                       f'-> {self.RULE_YEARS:.0f}-year rule '
                       f'{"COMPLIANT" if compliant else "NOT met"}; immediate deorbit '
                       f'(perigee {self.reentry_altitude / 1000:.0f} km) delta-v '
                       f'{dv_deorbit:.1f} m/s')

        dv_25y, sma_25y = 0.0, sma_start
        if not compliant:
            # Circular disposal altitude with a 25-year lifetime, by bisection
            lo, hi = reentry_r, sma_start
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if self._lifetime_years(mid) > self.RULE_YEARS:
                    hi = mid
                else:
                    lo = mid
            sma_25y = lo
            dv_25y = self._hohmann_dv(sma_start, sma_25y)
            ls.logger.info(f'{self.type}: {self.RULE_YEARS:.0f}-year compliant circular '
                           f'altitude {(sma_25y - R_EARTH) / 1000:.1f} km, transfer '
                           f'delta-v {dv_25y:.1f} m/s')

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(times, alts / 1000.0, 'b-', linewidth=1.2, label='Predicted decay')
        ax.axhline(self.reentry_altitude / 1000.0, color='red', linestyle=':',
                   linewidth=1.0, label=f'Re-entry ({self.reentry_altitude / 1000:.0f} km)')
        if not compliant:
            ax.axhline((sma_25y - R_EARTH) / 1000.0, color='green', linestyle='--',
                       linewidth=1.0,
                       label=f'{self.RULE_YEARS:.0f}-year altitude '
                             f'({(sma_25y - R_EARTH) / 1000:.0f} km, {dv_25y:.0f} m/s)')
        ax.set_xlabel('Time after simulation start [years]')
        ax.set_ylabel('Circular-orbit altitude [km]')
        ax.grid(True)
        ax.legend(fontsize=9)
        sat_id = sm.satellites[self.idx_found_satellite].sat_id
        fig.suptitle(f'Orbital lifetime, satellite {sat_id}: {lifetime_str} '
                     f'({self.RULE_YEARS:.0f}-year rule '
                     f'{"compliant" if compliant else "NOT met"}), '
                     f'deorbit delta-v {dv_deorbit:.0f} m/s')
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['years', 'altitude_m'], np.column_stack([times, alts]))


# ---------------------------------------------------------------------------
# Impulsive maneuver delta-v budgets (orb_deltav_injection / _reentry /
# _collision): vis-viva calculators on the orbit of the satellite block
# ---------------------------------------------------------------------------

def _vis_viva(radius, sma):
    """Orbital speed [m/s] at radius on an orbit with semi-major axis sma."""
    return np.sqrt(GM_EARTH * (2.0 / radius - 1.0 / sma))


class _AnalysisDeltaVBase(AnalysisBase):
    """Shared satellite selection of the impulsive delta-v budget analyses:
    the target/nominal orbit is the one given in the <Satellite> block (or
    the TLE) of the selected satellite."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection
        self.idx_found_satellite = 0

    def _read_selection(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)

    def before_loop(self, sm):
        for idx_sat, satellite in enumerate(sm.satellites):
            if self.constellation_id > 0 and \
                    satellite.constellation_id != self.constellation_id:
                continue
            if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:
                continue
            self.idx_found_satellite = idx_sat
            break

    def _orbit(self, sm):
        """(semi-major axis, apogee radius, perigee radius) of the selected
        satellite's nominal orbit."""
        kepler = sm.satellites[self.idx_found_satellite].kepler
        sma = kepler.semi_major_axis
        return sma, sma * (1.0 + kepler.eccentricity), \
            sma * (1.0 - kepler.eccentricity)


# Typical 3-sigma LEO/SSO injection accuracies (launcher user manuals, first
# order): (name, altitude error [m], inclination error [deg]). Overridable
# with <Launcher>Name, delta_alt_m, delta_inc_deg</Launcher> config tags.
_LAUNCHER_PRESETS = [
    ('Ariane 62', 5000.0, 0.04),
    ('Vega-C', 15000.0, 0.15),
    ('Falcon 9', 15000.0, 0.10),
]


class AnalysisOrbDeltaVInjection(_AnalysisDeltaVBase):
    """Delta-v to correct typical launcher injection errors, with the orbit
    of the <Satellite> block as the target: per launcher the 3-sigma
    altitude error (corrected with tangential burns at both apsides,
    dv = v*dh/(2a)) and inclination error (plane change dv = 2 v sin(di/2))
    are costed and stacked. Reported as the conservative arithmetic sum and
    the root-sum-square (a combined burn corrects both cheaper). Defaults:
    typical Ariane 62 / Vega-C / Falcon 9 LEO accuracies; override with
    <Launcher>Name, delta_alt_m, delta_inc_deg</Launcher> tags."""

    def __init__(self):
        super().__init__()
        self.launchers = []  # (name, altitude error [m], inclination error [rad])

    def read_config(self, node):
        self._read_selection(node)
        for launcher in node.findall('Launcher'):
            name, d_alt, d_inc = [v.strip() for v in launcher.text.split(',')]
            self.launchers.append((name, float(d_alt), radians(float(d_inc))))
        if not self.launchers:
            self.launchers = [(name, d_alt, radians(d_inc))
                              for name, d_alt, d_inc in _LAUNCHER_PRESETS]

    def after_loop(self, sm):
        sma, _, _ = self._orbit(sm)
        v_circ = np.sqrt(GM_EARTH / sma)
        rows = []
        for idx, (name, d_alt, d_inc) in enumerate(self.launchers):
            dv_alt = v_circ * d_alt / (2.0 * sma)  # Both apsides corrected
            dv_inc = 2.0 * v_circ * np.sin(d_inc / 2.0)
            dv_sum = dv_alt + dv_inc
            dv_rss = np.hypot(dv_alt, dv_inc)
            rows.append([idx + 1, d_alt, degrees(d_inc), dv_alt, dv_inc,
                         dv_sum, dv_rss])
            ls.logger.info(f'{self.type}: {name}: altitude +/-{d_alt / 1000:.0f} km '
                           f'-> {dv_alt:.1f} m/s, inclination '
                           f'+/-{degrees(d_inc):.2f} deg -> {dv_inc:.1f} m/s, '
                           f'total {dv_sum:.1f} m/s (RSS {dv_rss:.1f} m/s)')
        rows = np.array(rows)

        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(self.launchers))
        width = 0.28
        ax.bar(x - width, rows[:, 3], width, label='Altitude correction', color='C0')
        ax.bar(x, rows[:, 4], width, label='Inclination correction', color='C1')
        ax.bar(x + width, rows[:, 5], width, label='Total (sum)', color='C3')
        ax.set_xticks(x)
        ax.set_xticklabels([f'{name}\n(±{d_alt / 1000:.0f} km, '
                            f'±{degrees(d_inc):.2f} deg)'
                            for name, d_alt, d_inc in self.launchers], fontsize=9)
        ax.set_ylabel('Delta-v [m/s]')
        ax.set_title(f'Injection error correction to the target orbit '
                     f'(altitude {(sma - R_EARTH) / 1000:.0f} km)')
        ax.grid(True, axis='y', alpha=0.4)
        ax.legend()
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['launcher_id', 'delta_alt_m', 'delta_inc_deg',
                            'dv_altitude_ms', 'dv_inclination_ms',
                            'dv_total_ms', 'dv_rss_ms'], rows)


class AnalysisOrbDeltaVReentry(_AnalysisDeltaVBase):
    """Delta-v of a typical ESA controlled re-entry in two steps: an apogee
    burn lowering the perigee to a safe intermediate altitude
    (<IntermediatePerigee>, default 250 km) and the final apogee burn
    lowering the perigee into the atmosphere for a targeted entry
    (<FinalPerigee>, default 50 km). Vis-viva at the (fixed) apogee of the
    orbit given in the <Satellite> block."""

    def __init__(self):
        super().__init__()
        self.intermediate_perigee = 250e3  # [m] safe intermediate perigee
        self.final_perigee = 50e3  # [m] targeted-entry perigee

    def read_config(self, node):
        self._read_selection(node)
        if node.find('IntermediatePerigee') is not None:
            self.intermediate_perigee = float(node.find('IntermediatePerigee').text)
        if node.find('FinalPerigee') is not None:
            self.final_perigee = float(node.find('FinalPerigee').text)

    def after_loop(self, sm):
        sma, r_apo, r_peri = self._orbit(sm)
        steps = [('nominal orbit', r_peri, 0.0)]
        for label, target in (('intermediate perigee', self.intermediate_perigee),
                              ('final entry perigee', self.final_perigee)):
            sma_before = (r_apo + steps[-1][1]) / 2.0
            r_target = R_EARTH + target
            sma_after = (r_apo + r_target) / 2.0
            dv = _vis_viva(r_apo, sma_before) - _vis_viva(r_apo, sma_after)
            steps.append((label, r_target, dv))
        total_dv = sum(dv for _, _, dv in steps)
        ls.logger.info(f'{self.type}: apogee {(r_apo - R_EARTH) / 1000:.0f} km, '
                       + ', '.join(f'{label} {(r - R_EARTH) / 1000:.0f} km '
                                   f'({dv:.1f} m/s)'
                                   for label, r, dv in steps[1:])
                       + f', total {total_dv:.1f} m/s')

        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(steps))
        apogees = np.full(len(steps), (r_apo - R_EARTH) / 1000.0)
        perigees = np.array([(r - R_EARTH) / 1000.0 for _, r, _ in steps])
        ax.plot(x, apogees, 'o-', color='C0', label='Apogee altitude')
        ax.plot(x, perigees, 'o-', color='C1', label='Perigee altitude')
        for i, (label, _, dv) in enumerate(steps):
            if dv > 0:
                ax.annotate(f'burn {i}: {dv:.1f} m/s', (i, perigees[i]),
                            xytext=(0, 12), textcoords='offset points',
                            ha='center', fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels([label for label, _, _ in steps], fontsize=9)
        ax.set_ylabel('Altitude [km]')
        ax.set_title(f'Controlled re-entry: two apogee burns, '
                     f'total {total_dv:.1f} m/s')
        ax.grid(True, alpha=0.4)
        ax.legend()
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['step', 'perigee_km', 'apogee_km', 'dv_ms'],
                       [[i, (r - R_EARTH) / 1000.0, (r_apo - R_EARTH) / 1000.0, dv]
                        for i, (_, r, dv) in enumerate(steps)])


class AnalysisOrbDeltaVCollision(_AnalysisDeltaVBase):
    """Delta-v of a short-term collision avoidance maneuver: a tangential
    burn half an orbit before the conjunction raises the orbit on the far
    side (the apogee over the conjunction point) by <AvoidanceAltitude>
    (default 10 km), and an equal burn afterwards brings the orbit back to
    nominal - the budget is twice the raise burn. Computed with vis-viva on
    the orbit of the <Satellite> block; the plot shows the cost versus the
    raise altitude with the configured value marked."""

    def __init__(self):
        super().__init__()
        self.avoidance_altitude = 10e3  # [m] apogee raise over the conjunction

    def read_config(self, node):
        self._read_selection(node)
        if node.find('AvoidanceAltitude') is not None:
            self.avoidance_altitude = float(node.find('AvoidanceAltitude').text)

    def _dv_raise(self, r_burn, sma, raise_m):
        sma_raised = sma + raise_m / 2.0
        return _vis_viva(r_burn, sma_raised) - _vis_viva(r_burn, sma)

    def after_loop(self, sm):
        sma, _, r_peri = self._orbit(sm)
        raises = np.linspace(0.0, 2.0 * self.avoidance_altitude, 41)
        dv_up = np.array([self._dv_raise(r_peri, sma, r) for r in raises])
        dv_config = self._dv_raise(r_peri, sma, self.avoidance_altitude)
        ls.logger.info(f'{self.type}: raising the orbit over the conjunction by '
                       f'{self.avoidance_altitude / 1000:.1f} km costs '
                       f'{dv_config:.2f} m/s, return to nominal the same: '
                       f'total {2 * dv_config:.2f} m/s per avoidance maneuver')

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(raises / 1000.0, 2.0 * dv_up, '-', color='C0',
                label='Total (raise + return)')
        ax.plot(raises / 1000.0, dv_up, '--', color='C1', label='Raise burn only')
        ax.plot(self.avoidance_altitude / 1000.0, 2.0 * dv_config, 'r^',
                markersize=10)
        ax.annotate(f'{2 * dv_config:.2f} m/s at '
                    f'{self.avoidance_altitude / 1000:.1f} km',
                    (self.avoidance_altitude / 1000.0, 2.0 * dv_config),
                    xytext=(8, -12), textcoords='offset points', fontsize=9)
        ax.set_xlabel('Altitude raise over the conjunction point [km]')
        ax.set_ylabel('Delta-v [m/s]')
        ax.set_title(f'Collision avoidance maneuver '
                     f'(orbit altitude {(sma - R_EARTH) / 1000:.0f} km)')
        ax.grid(True, alpha=0.4)
        ax.legend()
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['raise_m', 'dv_up_ms', 'dv_total_ms'],
                       np.column_stack([raises, dv_up, 2.0 * dv_up]))


# ---------------------------------------------------------------------------
# Space environment models (orb_environment). First-order engineering
# parametrisations calibrated to published model outputs (AE8/AP8 LEO maps,
# MSIS mean-activity atomic oxygen, Gruen 1985 micrometeoroids); they place
# the SAA and the polar horns correctly and give order-of-magnitude fluxes
# and doses. For design/qualification numbers run SPENVIS or IRENE (AE9/AP9).
# ---------------------------------------------------------------------------

# Eccentric tilted dipole of the geomagnetic field: axis toward the
# geomagnetic north pole (IGRF-13 2020) and the dipole centre displaced from
# the Earth centre toward the west Pacific - the displacement is what makes
# the field anomalously weak over the South Atlantic (the SAA).
_DIPOLE_B0 = 3.12e-5  # [T] surface field at the magnetic equator (SMAD)
_DIPOLE_POLE = (radians(80.7), radians(-72.7))  # geomagnetic north pole lat, lon
_DIPOLE_OFFSET_M = 560e3  # [m] eccentric dipole centre displacement
_DIPOLE_OFFSET_DIR = (radians(22.5), radians(141.5))  # displacement direction lat, lon

# Trapped-particle flux parametrisation, omnidirectional integral flux
# [cm^-2 s^-1]: a Gaussian radial belt profile in L, a power-law decay in
# B/B_eq (off-equator attenuation along the field line) and a drift-loss
# gate - inner-belt particles survive only where their drift shell stays
# above the atmosphere at the weakest-field longitude, which the eccentric
# dipole offset turns into the South Atlantic Anomaly. Calibrated to AE8/AP8
# LEO maps: protons >10 MeV ~1e3 in the SAA core at 800 km, electrons >1 MeV
# a few 1e4 in the polar horns.
_PROTON_PEAK = 2.0e3   # >10 MeV low-x amplitude of the inner belt (L=1.45)
_ELECTRON_INNER = 2.0e5  # >1 MeV inner belt amplitude at L=1.55
_ELECTRON_OUTER = 5.0e6  # >1 MeV outer belt equatorial peak at L=4.8
_FLUX_FLOOR = 1.0  # [cm^-2 s^-1] fluxes below this are reported as zero

# Fluence-to-dose conversion behind aluminium shielding [rad cm^2 /particle]:
# protons >10 MeV penetrate several mm with slow attenuation, electrons
# >1 MeV are stopped in the first few mm (residual = bremsstrahlung floor).
_SHIELD_MM = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0])
_DOSE_PER_PROTON = lambda t_mm: 2.0e-7 * np.exp(-t_mm / 12.0)
_DOSE_PER_ELECTRON = lambda t_mm: 8.0e-8 * np.exp(-t_mm / 0.8) + 1.0e-10 * np.exp(-t_mm / 40.0)

# Atomic oxygen number density [m^-3] versus altitude [km], mean solar
# activity (MSIS-class mean profile), log-interpolated between the nodes
_AO_TABLE = np.array([
    [150, 3.0e16], [200, 5.0e15], [250, 1.2e15], [300, 3.5e14], [350, 1.1e14],
    [400, 3.5e13], [450, 1.3e13], [500, 5.0e12], [550, 2.0e12], [600, 9.0e11],
    [650, 4.0e11], [700, 1.8e11], [750, 8.0e10], [800, 3.5e10], [900, 8.0e9],
    [1000, 2.0e9]])
_AO_KAPTON_YIELD = 3.0e-24  # [cm^3/atom] kapton erosion yield (LEO standard)


def dipole_l_b(pos_ecf):
    """McIlwain L-shell [-], field strength B [T] and distance to the dipole
    centre [m] at an ECF position, from the eccentric tilted dipole (the
    offset centre is what creates the SAA)."""
    pole_lat, pole_lon = _DIPOLE_POLE
    m_hat = np.array([cos(pole_lat) * cos(pole_lon),
                      cos(pole_lat) * sin(pole_lon), sin(pole_lat)])
    off_lat, off_lon = _DIPOLE_OFFSET_DIR
    centre = _DIPOLE_OFFSET_M * np.array([cos(off_lat) * cos(off_lon),
                                          cos(off_lat) * sin(off_lon), sin(off_lat)])
    r_vec = np.asarray(pos_ecf, dtype=float) - centre
    r = np.linalg.norm(r_vec)
    sin_maglat = np.clip(np.dot(r_vec / r, m_hat), -1.0, 1.0)
    cos2_maglat = max(1.0 - sin_maglat ** 2, 1e-6)
    l_shell = (r / R_EARTH) / cos2_maglat
    b_field = _DIPOLE_B0 * (R_EARTH / r) ** 3 * np.sqrt(1.0 + 3.0 * sin_maglat ** 2)
    return l_shell, b_field, r


def trapped_flux(l_shell, b_field, r_dipole):
    """AE8/AP8-style omnidirectional integral flux estimates [cm^-2 s^-1]:
    (protons >10 MeV, electrons >1 MeV) at a point (l_shell, b_field) at
    distance r_dipole from the dipole centre. The drift-loss gate keeps
    inner-belt flux only where the drift shell clears the atmosphere at the
    weakest-field longitude: the same field point sits an offset-distance
    lower over the SAA side, so at LEO only the SAA region (and a smooth
    fringe) sees the inner belts. First-order engineering model."""
    b_eq = _DIPOLE_B0 / l_shell ** 3  # equatorial field on the drift shell
    x = max(b_field / b_eq, 1.0)  # off-equator attenuation parameter
    # Minimum geodetic altitude of this drift shell (weakest-field longitude);
    # smooth 100..300 km ramp instead of a hard atmospheric cutoff
    drift_alt_min_km = (r_dipole - _DIPOLE_OFFSET_M - R_EARTH) / 1000.0
    gate = min(max((drift_alt_min_km - 100.0) / 200.0, 0.0), 1.0)
    j_proton = _PROTON_PEAK * np.exp(-0.5 * ((l_shell - 1.45) / 0.25) ** 2) * x ** -3.0 * gate
    j_electron = (_ELECTRON_INNER * np.exp(-0.5 * ((l_shell - 1.55) / 0.30) ** 2) * x ** -3.0 * gate +
                  _ELECTRON_OUTER * np.exp(-0.5 * ((l_shell - 4.80) / 0.80) ** 2) * x ** -0.9)
    if j_proton < _FLUX_FLOOR:
        j_proton = 0.0
    if j_electron < _FLUX_FLOOR:
        j_electron = 0.0
    return j_proton, j_electron


def atomic_oxygen_density(altitude_m):
    """Atomic oxygen number density [m^-3] at altitude, mean solar activity
    (log-interpolation of an MSIS-class mean profile)."""
    h_km = altitude_m / 1000.0
    if h_km <= _AO_TABLE[0, 0]:
        return _AO_TABLE[0, 1]
    if h_km >= _AO_TABLE[-1, 0]:  # extrapolate with the last scale height
        h0, h1 = _AO_TABLE[-2, 0], _AO_TABLE[-1, 0]
        scale = (h1 - h0) / np.log(_AO_TABLE[-2, 1] / _AO_TABLE[-1, 1])
        return _AO_TABLE[-1, 1] * np.exp(-(h_km - h1) / scale)
    return float(np.exp(np.interp(h_km, _AO_TABLE[:, 0], np.log(_AO_TABLE[:, 1]))))


def gruen_flux(mass_g):
    """Gruen et al. (1985) interplanetary micrometeoroid flux: cumulative
    number of particles with mass > mass_g [g] per m^2 per second at 1 AU
    (random tumbling plate)."""
    m = np.asarray(mass_g, dtype=float)
    return ((2.2e3 * m ** 0.306 + 15.0) ** -4.38 +
            1.3e-9 * (m + 1.0e11 * m ** 2 + 1.0e27 * m ** 4) ** -0.36 +
            1.3e-16 * (m + 1.0e6 * m ** 2) ** -0.85)


class AnalysisOrbEnvironment(AnalysisBase):
    """Space environment along the orbit, SPENVIS-style summary sheet:
    trapped radiation (eccentric-dipole L/B with AE8/AP8-style flux
    estimates, SAA and polar horn crossings, annual fluences and a total
    ionizing dose versus aluminium shielding curve), atomic oxygen (density
    at the orbit altitude, ram fluence and kapton erosion over the mission)
    and Gruen micrometeoroid flux with the expected impact count on the
    spacecraft area over the mission. All models are first-order engineering
    estimates - use SPENVIS/IRENE for design values."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection
        self.mission_years = 5.0  # Fluence/erosion/impact accumulation period
        self.surface_area = 10.0  # [m^2] area exposed to micrometeoroids
        self.idx_found_satellite = 0
        self.metric = None  # (num_epoch, 8): lat, lon, alt, L, B, j_p, j_e, n_AO

    def read_config(self, node):
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('SatelliteID') is not None:
            self.satellite_id = int(node.find('SatelliteID').text)
        if node.find('MissionYears') is not None:
            self.mission_years = float(node.find('MissionYears').text)
        if node.find('SurfaceArea') is not None:
            self.surface_area = float(node.find('SurfaceArea').text)

    def before_loop(self, sm):
        self.idx_found_satellite = 0
        for idx_sat, satellite in enumerate(sm.satellites):
            if self.constellation_id > 0 and \
                    satellite.constellation_id != self.constellation_id:
                continue
            if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:
                continue
            self.idx_found_satellite = idx_sat
            break
        self.metric = np.zeros((sm.num_epoch, 8))

    def in_loop(self, sm):
        satellite = sm.satellites[self.idx_found_satellite]
        satellite.det_lla()
        pos_ecf = np.asarray(satellite.pos_ecf, dtype=float)
        altitude = np.linalg.norm(pos_ecf) - misc_fn.earth_radius_lat(satellite.lla[0])
        l_shell, b_field, r_dipole = dipole_l_b(pos_ecf)
        j_proton, j_electron = trapped_flux(l_shell, b_field, r_dipole)
        self.metric[sm.cnt_epoch] = [degrees(satellite.lla[0]), degrees(satellite.lla[1]),
                                     altitude / 1000.0, l_shell, b_field,
                                     j_proton, j_electron,
                                     atomic_oxygen_density(altitude)]

    def after_loop(self, sm):
        times = np.asarray(self.times_f_doy)
        lat, lon, alt_km = self.metric[:, 0], self.metric[:, 1], self.metric[:, 2]
        l_shell, b_gauss = self.metric[:, 3], self.metric[:, 4] * 1e4
        j_p, j_e = self.metric[:, 5], self.metric[:, 6]
        n_ao = self.metric[:, 7]
        seconds_per_year = 365.25 * 86400.0

        # Annual fluences from the orbit-average fluxes, dose vs shielding
        fluence_p_yr = j_p.mean() * seconds_per_year
        fluence_e_yr = j_e.mean() * seconds_per_year
        dose_p = fluence_p_yr * _DOSE_PER_PROTON(_SHIELD_MM)
        dose_e = fluence_e_yr * _DOSE_PER_ELECTRON(_SHIELD_MM)
        dose_total = dose_p + dose_e

        # Atomic oxygen: ram flux [cm^-2 s^-1] at orbital velocity, mission
        # fluence and kapton erosion depth
        v_orbit = np.sqrt(GM_EARTH / (R_EARTH + alt_km.mean() * 1000.0))
        ao_flux = n_ao.mean() * v_orbit * 1e-4  # [cm^-2 s^-1]
        ao_fluence_mission = ao_flux * self.mission_years * seconds_per_year
        ao_erosion_um = ao_fluence_mission * _AO_KAPTON_YIELD * 1e4

        # Micrometeoroids: Gruen flux with Earth shielding and gravitational
        # focusing at the orbit altitude, expected impacts over the mission
        r_orbit = R_EARTH + alt_km.mean() * 1000.0
        sin_eta = min((R_EARTH + 100e3) / r_orbit, 1.0)  # Earth + atmosphere shield
        shielding = (1.0 + np.sqrt(1.0 - sin_eta ** 2)) / 2.0
        focusing = 1.0 + R_EARTH / r_orbit * 0.76  # v_esc^2/v_inf^2 ~ 0.76 R/r (20 km/s)
        masses = np.logspace(-9, 0, 46)
        mm_flux_m2yr = gruen_flux(masses) * shielding * focusing * seconds_per_year
        mm_impacts = mm_flux_m2yr * self.surface_area * self.mission_years

        # In the SAA when the inner-belt proton flux is substantial (an order
        # of magnitude below the core value, i.e. not the drift-loss fringe)
        in_saa = (j_p > 100.0) & (np.abs(lat) < 45.0)
        saa_fraction = in_saa.mean() * 100.0
        idx_4mm = int(np.argmin(np.abs(_SHIELD_MM - 4.0)))
        impacts_1ug = float(np.interp(1e-6, masses, mm_impacts))
        ls.logger.info(f'{self.type}: SAA crossings {saa_fraction:.1f}% of the time, '
                       f'orbit-average flux protons >10 MeV {j_p.mean():.1f}, '
                       f'electrons >1 MeV {j_e.mean():.1f} /cm2/s')
        ls.logger.info(f'{self.type}: annual fluence protons {fluence_p_yr:.2e}, '
                       f'electrons {fluence_e_yr:.2e} /cm2/yr, TID behind 4 mm Al '
                       f'{dose_total[idx_4mm]:.0f} rad/yr '
                       f'({dose_total[idx_4mm] * self.mission_years / 1000.0:.1f} krad '
                       f'over {self.mission_years:.0f} years)')
        ls.logger.info(f'{self.type}: atomic oxygen {n_ao.mean():.2e} /m3, ram fluence '
                       f'{ao_fluence_mission:.2e} /cm2 over the mission, kapton erosion '
                       f'{ao_erosion_um:.1f} um')
        ls.logger.info(f'{self.type}: micrometeoroid impacts >1 ug over the mission on '
                       f'{self.surface_area:.0f} m2: {impacts_1ug:.1f}')

        # Six separate plots (<type>_trapped_flux.png, ...) instead of a
        # single six-panel environment sheet. All panels are first-order
        # engineering models - use SPENVIS/IRENE for design values.
        fig, ax1 = make_map_cyl()
        self.decorate_map2d(sm, ax1)
        total_flux = j_p + j_e
        quiet = total_flux <= 0
        ax1.scatter(lon[quiet], lat[quiet], s=1, c='lightgrey',
                    transform=ccrs.PlateCarree())
        sc = ax1.scatter(lon[~quiet], lat[~quiet], s=4, c=np.log10(total_flux[~quiet]),
                         cmap=plt.cm.jet, transform=ccrs.PlateCarree())
        plt.colorbar(sc, ax=ax1, shrink=0.85, label='log10 trapped flux [/cm2/s]')
        ax1.set_title('Trapped radiation: SAA and polar horn crossings '
                      '(first-order model)')
        plt.savefig(sm.output_path(self.type + '_trapped_flux.png'))
        plt.show()

        fig, ax2 = plt.subplots(figsize=(10, 6))
        ax2.semilogy(times, np.where(j_p > 0, j_p, np.nan), 'r.', markersize=2,
                     label='Protons >10 MeV')
        ax2.semilogy(times, np.where(j_e > 0, j_e, np.nan), 'b.', markersize=2,
                     label='Electrons >1 MeV')
        ax2.set_xlabel('Day of Year (DOY)')
        ax2.set_ylabel('Flux [/cm2/s]')
        ax2.set_title('Trapped particle flux (first-order model)')
        ax2.grid(True, alpha=0.4)
        ax2.legend(fontsize=8)
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '_flux_timeseries.png'))
        plt.show()

        fig, ax3 = plt.subplots(figsize=(10, 6))
        ax3.plot(times, l_shell, 'g-', linewidth=0.6)
        ax3.set_xlabel('Day of Year (DOY)')
        ax3.set_ylabel('L-shell [-]', color='g')
        ax3.set_ylim(1, min(l_shell.max() * 1.1, 12))
        ax3.grid(True, alpha=0.4)
        ax3b = ax3.twinx()
        ax3b.plot(times, b_gauss, 'm-', linewidth=0.6)
        ax3b.set_ylabel('B [gauss]', color='m')
        ax3.set_title('Magnetic drift-shell coordinates (eccentric dipole)')
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '_drift_shell.png'))
        plt.show()

        fig, ax4 = plt.subplots(figsize=(10, 6))
        ax4.semilogy(_SHIELD_MM, dose_total, 'k-o', markersize=4, label='Total')
        ax4.semilogy(_SHIELD_MM, dose_p, 'r--', label='Trapped protons')
        ax4.semilogy(_SHIELD_MM, dose_e, 'b--', label='Trapped electrons')
        ax4.set_xlabel('Aluminium shielding [mm]')
        ax4.set_ylabel('Total ionizing dose [rad/yr]')
        ax4.set_title('Dose-depth curve (first-order model - use SPENVIS/IRENE '
                      'for design values)')
        ax4.grid(True, which='both', alpha=0.4)
        ax4.legend(fontsize=8)
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '_dose_depth.png'))
        plt.show()

        fig, ax5 = plt.subplots(figsize=(10, 6))
        years = np.linspace(0.0, self.mission_years, 50)
        ax5.plot(years, years * ao_erosion_um / self.mission_years, 'g-')
        ax5.set_xlabel('Mission time [years]')
        ax5.set_ylabel('Kapton erosion depth [um]')
        ax5.set_title(f'Atomic oxygen erosion (ram flux {ao_flux:.2e} /cm2/s)')
        ax5.grid(True, alpha=0.4)
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '_atomic_oxygen.png'))
        plt.show()

        fig, ax6 = plt.subplots(figsize=(10, 6))
        ax6.loglog(masses, mm_flux_m2yr, 'b-', label='Flux [/m2/yr]')
        ax6.loglog(masses, mm_impacts, 'r--',
                   label=f'Impacts, {self.surface_area:.0f} m2, '
                         f'{self.mission_years:.0f} yr')
        ax6.set_xlabel('Particle mass [g]')
        ax6.set_ylabel('Cumulative flux / impact count')
        ax6.set_title('Micrometeoroids (Gruen 1985)')
        ax6.grid(True, which='both', alpha=0.4)
        ax6.legend(fontsize=8)
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '_micrometeoroids.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'lat_deg', 'lon_deg', 'alt_km', 'l_shell', 'b_gauss',
                            'proton_flux_cm2s', 'electron_flux_cm2s', 'ao_density_m3'],
                       np.column_stack([times, lat, lon, alt_km, l_shell, b_gauss,
                                        j_p, j_e, n_ao]))
        self.write_csv(sm, ['shielding_mm', 'dose_protons_rad_yr', 'dose_electrons_rad_yr',
                            'dose_total_rad_yr'],
                       np.column_stack([_SHIELD_MM, dose_p, dose_e, dose_total]),
                       suffix='dose')
        self.write_csv(sm, ['mass_g', 'flux_m2yr', 'impacts_mission'],
                       np.column_stack([masses, mm_flux_m2yr, mm_impacts]),
                       suffix='meteoroid')
