import numpy as np
from math import degrees
import matplotlib.pyplot as plt

# Project modules
from constants import GM_EARTH, R_EARTH
from analysis import AnalysisBase
import logging_svs as ls


class AnalysisOrbSemiMajorAxis(AnalysisBase):
    """Osculating semi-major axis over the simulation time, computed from the
    ECI state vector with the vis-viva equation a = 1/(2/r - v^2/mu). With the
    HPOP propagator and drag enabled this shows the orbital decay; the analysis
    also works with the other propagators (constant a for Keplerian, mean-element
    variations for SGP4)."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection

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
        for satellite in sm.satellites:
            satellite.metric = np.full(sm.num_epoch, np.nan)

    def in_loop(self, sm):
        for satellite in sm.satellites:
            if not self._selected(satellite):
                continue
            r = np.linalg.norm(satellite.pos_eci)
            v2 = float(np.dot(satellite.vel_eci, satellite.vel_eci))
            satellite.metric[sm.cnt_epoch] = 1.0 / (2.0 / r - v2 / GM_EARTH)

    def after_loop(self, sm):
        fig = plt.figure(figsize=(10, 6))
        plotted = False
        for satellite in sm.satellites:
            if not self._selected(satellite):
                continue
            sma_km = satellite.metric / 1000.0
            # Secular change: difference of the mean SMA over the first and last
            # ~100 epochs, which averages out the J2 short-period oscillation
            n_avg = min(100, len(sma_km))
            first, last = np.nanmean(sma_km[:n_avg]), np.nanmean(sma_km[-n_avg:])
            plt.plot(self.times_f_doy, sma_km, '-',
                     label=f'Sat {satellite.sat_id} (secular change {(last-first)*1000:.0f} m)')
            ls.logger.info(f'Satellite {satellite.sat_id}: mean SMA first epochs {first:.3f} km, '
                           f'last epochs {last:.3f} km, change {(last-first)*1000:.1f} m')
            plotted = True
        if not plotted:
            ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                            f'SatelliteID {self.satellite_id}. No plot produced.')
            return
        plt.xlabel('DOY [-]')
        plt.ylabel('Osculating semi-major axis [km]')
        plt.grid()
        plt.legend()
        plt.savefig('../output/' + self.type + '.png')
        plt.show()
