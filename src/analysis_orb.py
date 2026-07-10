from math import degrees, radians

import numpy as np
import matplotlib.pyplot as plt

# Project modules
from constants import GM_EARTH, OMEGA_EARTH, R_EARTH
from analysis import AnalysisBase
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
        fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
        panels = [  # (key, scale factor, label)
            ('sma', 1e-3, 'Semi-major axis [km]'),
            ('ecc', 1.0, 'Eccentricity [-]'),
            ('incl', np.degrees(1.0), 'Inclination [deg]'),
            ('raan', np.degrees(1.0), 'RAAN [deg]'),
            ('arg_perigee', np.degrees(1.0), 'Argument of perigee [deg]'),
            ('mean_anomaly', np.degrees(1.0), 'Mean anomaly [deg]'),
        ]
        times = np.asarray(self.times_f_doy)
        plotted = False
        for idx_sat, satellite in enumerate(sm.satellites):
            if not self._selected(satellite):
                continue
            used = ~np.isnan(self.sat_metric[idx_sat, :, 0])
            if not used.any():
                continue
            elements = rv2kepler(self.sat_metric[idx_sat, used, 0:3], self.sat_metric[idx_sat, used, 3:6])
            sma_km = elements['sma'] / 1000.0
            # Secular change: difference of the mean over the first and last
            # ~100 epochs, which averages out the J2 short-period oscillation
            n_avg = min(100, len(sma_km))
            first, last = np.mean(sma_km[:n_avg]), np.mean(sma_km[-n_avg:])
            ls.logger.info(f'Satellite {satellite.sat_id}: mean SMA first epochs {first:.3f} km, '
                           f'last epochs {last:.3f} km, change {(last-first)*1000:.1f} m')
            for ax, (key, scale, label) in zip(axes.flat, panels):
                ax.plot(times[used], elements[key] * scale, '-', linewidth=0.9,
                        label=f'Sat {satellite.sat_id}')
            plotted = True
        if not plotted:
            ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                            f'SatelliteID {self.satellite_id}. No plot produced.')
            return
        for ax, (key, scale, label) in zip(axes.flat, panels):
            ax.set_ylabel(label)
            ax.grid(True)
        for ax in axes[-1, :]:
            ax.set_xlabel('DOY [-]')
        axes.flat[0].legend(fontsize=8)
        fig.suptitle('Osculating Kepler elements')
        fig.tight_layout()
        plt.savefig('../output/' + self.type + '.png')
        plt.show()


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
        plt.savefig('../output/' + self.type + '.png')
        plt.show()


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
        plt.savefig('../output/' + self.type + '.png')
        plt.show()


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
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        ax1.plot(times, self.metric[:, 0], 'r-', label='xp')
        ax1.plot(times, self.metric[:, 1], 'b-', label='yp')
        ax1.set_xlabel('DOY [-]')
        ax1.set_ylabel('Polar motion [arcsec]')
        ax1.grid(True)
        ax1.legend()
        ax2.plot(self.metric[:, 0], self.metric[:, 1], 'g-')
        ax2.plot(self.metric[0, 0], self.metric[0, 1], 'go', label='start')
        ax2.plot(self.metric[-1, 0], self.metric[-1, 1], 'rs', label='end')
        ax2.set_xlabel('xp [arcsec]')
        ax2.set_ylabel('yp [arcsec]')
        ax2.set_aspect('equal', adjustable='datalim')
        ax2.grid(True)
        ax2.legend()
        fig.suptitle('Earth pole wobble (IERS polar motion)')
        fig.tight_layout()
        plt.savefig('../output/' + self.type + '.png')
        plt.show()


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
        plt.savefig('../output/' + self.type + '.png')
        plt.show()
