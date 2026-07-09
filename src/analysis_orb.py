import numpy as np
import matplotlib.pyplot as plt

# Project modules
from constants import GM_EARTH, OMEGA_EARTH
from analysis import AnalysisBase
import logging_svs as ls


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
