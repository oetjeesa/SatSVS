import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import sin, cos, asin, degrees, radians, log10
import itur
import astropy.units as u
from astropy.time import Time
import ast
# Project modules
from constants import K_BOLTZMANN, C_LIGHT
from analysis import AnalysisBase
import antenna
import misc_fn
import logging_svs as ls


def _off_nadir_angle(satellite, station):
    """Angle at the satellite between nadir and the line of sight to the
    station: the off-boresight angle of a nadir-pointed satellite antenna."""
    los = station.pos_ecf - satellite.pos_ecf
    return misc_fn.angle_two_vectors(los, -satellite.pos_ecf,
                                     norm(los), norm(satellite.pos_ecf))


class AnalysisComGr2SpBudget(AnalysisBase):

    def __init__(self):
        super().__init__()
        self.station_id = None
        self.transmitter_object = None
        self.carrier_frequency = None
        self.transmit_power = None
        self.transmit_losses = None
        self.transmit_gain = None
        self.p_exceed = None
        self.include_rain = None
        self.include_gas = None
        self.include_scintillation = None
        self.include_clouds = None
        self.receive_gain = None
        self.receive_losses = None
        self.receive_temp = None
        self.modulation_type = None
        self.ber = None
        self.data_rate = None
        self.tx_pattern_file = None  # Optional GRASP .cut/.grd antenna patterns
        self.rx_pattern_file = None
        self.tx_pattern = None
        self.rx_pattern = None
        self.sat_pointing = 'nadir'  # Satellite antenna: nadir-fixed or tracking

        self.idx_found_station = None

        self.metric = None
        self.cn0_required = 0

    def read_config(self, node):
        if node.find('GroundStationID') is not None:
            self.station_id = int(node.find('GroundStationID').text)
        if node.find('TransmitterObject') is not None:
            self.transmitter_object = node.find('TransmitterObject').text.lower()
        if node.find('CarrierFrequency') is not None:
            self.carrier_frequency = float(node.find('CarrierFrequency').text)

        if node.find('TransmitPowerW') is not None:
            self.transmit_power = float(node.find('TransmitPowerW').text)
        if node.find('TransmitLossesdB') is not None:
            self.transmit_losses = float(node.find('TransmitLossesdB').text)
        if node.find('TransmitGaindB') is not None:
            self.transmit_gain = float(node.find('TransmitGaindB').text)
        if node.find('TransmitAntennaPatternFile') is not None:
            self.tx_pattern_file = node.find('TransmitAntennaPatternFile').text
        if node.find('ReceiveAntennaPatternFile') is not None:
            self.rx_pattern_file = node.find('ReceiveAntennaPatternFile').text
        if node.find('SatelliteAntennaPointing') is not None:
            self.sat_pointing = node.find('SatelliteAntennaPointing').text.lower()

        if node.find('PExceedPerc') is not None:
            self.p_exceed = float(node.find('PExceedPerc').text)
        if node.find('IncludeRain') is not None:
            self.include_rain = misc_fn.str2bool(node.find('IncludeRain').text)
        if node.find('IncludeGas') is not None:
            self.include_gas = misc_fn.str2bool(node.find('IncludeGas').text)
        if node.find('IncludeScintillation') is not None:
            self.include_scintillation = misc_fn.str2bool(node.find('IncludeScintillation').text)
        if node.find('IncludeClouds') is not None:
            self.include_clouds = misc_fn.str2bool(node.find('IncludeClouds').text)

        if node.find('ReceiveGaindB') is not None:
            self.receive_gain = float(node.find('ReceiveGaindB').text)
        if node.find('ReceiveLossesdB') is not None:
            self.receive_losses = float(node.find('ReceiveLossesdB').text)
        if node.find('ReceiveTempK') is not None:
            self.receive_temp = float(node.find('ReceiveTempK').text)

        if node.find('ModulationType') is not None:
            self.modulation_type = node.find('ModulationType').text
        if node.find('BitErrorRate') is not None:
            self.ber = float(node.find('BitErrorRate').text)
        if node.find('DataRateBitPerSec') is not None:
            self.data_rate = float(node.find('DataRateBitPerSec').text)

    def before_loop(self, sm):
        for idx_station, station in enumerate(sm.stations):
            if station.station_id == self.station_id:
                self.idx_found_station = idx_station
                break
        if self.tx_pattern_file is not None:
            self.tx_pattern = antenna.AntennaPattern.from_file(self.tx_pattern_file)
        if self.rx_pattern_file is not None:
            self.rx_pattern = antenna.AntennaPattern.from_file(self.rx_pattern_file)
        if self.tx_pattern is not None or self.rx_pattern is not None:
            antenna.plot_patterns(
                [(label, p) for label, p in (('Tx', self.tx_pattern),
                                             ('Rx', self.rx_pattern)) if p is not None],
                sm.output_path(self.type + '_antenna.png'))
        self.metric = np.zeros((sm.num_epoch, 10))
        if self.modulation_type is not None:
            self.cn0_required = misc_fn.comp_cn0_required(self.modulation_type, self.ber, self.data_rate)

    def in_loop(self, sm):
        idx_station = self.idx_found_station
        for idx_sat in sm.stations[idx_station].idx_sat_in_view:
            elevation = sm.gr2sp[idx_station][idx_sat].elevation
            distance = sm.gr2sp[idx_station][idx_sat].distance
            # Antenna gains: with a pattern file the ground antenna tracks the
            # satellite (peak gain); the satellite antenna is nadir-pointed
            # (gain at the epoch's off-nadir angle, right for isoflux/horn
            # antennas) or, with SatelliteAntennaPointing Tracking, steered at
            # the station (peak gain). The scalar config gains are used
            # otherwise
            tx_gain, rx_gain = self.transmit_gain, self.receive_gain
            if self.tx_pattern is not None or self.rx_pattern is not None:
                off_sat = 0.0 if self.sat_pointing == 'tracking' else \
                    _off_nadir_angle(sm.satellites[idx_sat], sm.stations[idx_station])
                if self.transmitter_object == 'satellite':
                    if self.tx_pattern is not None:
                        tx_gain = self.tx_pattern.gain(off_sat)
                    if self.rx_pattern is not None:
                        rx_gain = self.rx_pattern.peak
                else:  # Station transmits (uplink)
                    if self.tx_pattern is not None:
                        tx_gain = self.tx_pattern.peak
                    if self.rx_pattern is not None:
                        rx_gain = self.rx_pattern.gain(off_sat)
            eirp = 10*log10(self.transmit_power) + tx_gain - self.transmit_losses
            fsl = 20*log10(distance/1000) + 20*log10(self.carrier_frequency/1e9) + 92.45
            # a_g = misc_fn.comp_gas_attenuation(self.carrier_frequency, elevation)  # Fast method but unaccurate <5 deg
            # fast method of the rain model
            # a_r = misc_fn.comp_rain_attenuation(self.carrier_frequency, elevation,
            #                                          sm.stations[idx_station].lla[0], sm.stations[idx_station].lla[2],
            #                                          self.rain_p_exceed, self.rain_rate, self.rain_height)
            a_g, a_c, a_r, a_s, a_t = itur.atmospheric_attenuation_slant_path(
                degrees(sm.stations[idx_station].lla[0]),
                degrees(sm.stations[idx_station].lla[1]),
                self.carrier_frequency / 1e9 * itur.u.GHz,
                degrees(elevation), self.p_exceed,
                D=1.0 * itur.u.m,
                include_gas=self.include_gas,
                include_rain=self.include_rain,
                include_scintillation=self.include_scintillation,
                include_clouds=self.include_clouds,
                return_contributions=True)
            if self.transmitter_object == 'satellite':
                ant_temp = 10 + misc_fn.temp_brightness(self.carrier_frequency, elevation)
            else:  # station is the transmitter
                ant_temp = 290
            temp_sys = self.receive_temp + ant_temp
            cn0 = eirp - fsl \
                  - a_g.value - a_r.value - a_c.value - a_s.value \
                  + rx_gain - self.receive_losses - \
                  K_BOLTZMANN - 10*np.log10(temp_sys)
            self.metric[sm.cnt_epoch,:] = [self.times_f_doy[sm.cnt_epoch], degrees(elevation),
                                           cn0, a_g.value, a_r.value, a_c.value, a_s.value, a_t.value,
                                           fsl, self.cn0_required]

    def after_loop(self, sm):
        self.metric = self.metric[~np.all(self.metric == 0, axis=1)]  # Clean up empty rows
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.subplots_adjust(left=.1, right=.95, top=0.95, bottom=0.07)
        time_list = self.metric[:, 0]
        plt.plot(time_list, self.metric[:, 1], 'k.', label='Elevation')
        plt.plot(time_list, self.metric[:, 2], 'b.', label='CN0 Computed')
        if self.include_gas:
            plt.plot(time_list, self.metric[:, 3], 'g+', label='Gas Attenuation')
        if self.include_rain:
            plt.plot(time_list, self.metric[:, 4], 'm+', label='Rain Attenuation')
        if self.include_clouds:
            plt.plot(time_list, self.metric[:, 5], 'r+', label='Cloud Attenuation')
        if self.include_scintillation:
            plt.plot(time_list, self.metric[:, 6], 'y+', label='Scintillation Attenuation')
        plt.plot(time_list, self.metric[:, 7], 'r.', label='Total Atmospheric Attenuation')
        plt.plot(time_list, self.metric[:, 8]-100, 'c+', label='Free space loss - 100dB')
        if self.modulation_type is not None:
            plt.plot(time_list, self.metric[:, 9], 'b+', label='CN0 Required')
        plt.xlabel('Day Of Year DOY [-]'); plt.ylabel('Elevation [deg], Power values [dB]')
        plt.legend(); plt.grid()
        plt.savefig(sm.output_path(self.type + '.png'))
        ax.ticklabel_format(useOffset=False, style='plain')
        plt.show()

        self.write_csv(sm, ['doy', 'elevation_deg', 'cn0_dbhz', 'gas_att_db', 'rain_att_db',
                            'cloud_att_db', 'scint_att_db', 'total_att_db', 'fsl_db',
                            'cn0_required_dbhz'], self.metric)


class AnalysisComSp2SpBudget(AnalysisBase):

    def __init__(self):
        super().__init__()
        self.sat_id1 = None
        self.sat_id2 = None

        self.carrier_frequency = None
        self.transmit_power = None
        self.transmit_losses = None
        self.transmit_gain = None
        self.receive_losses = None
        self.receive_temp = None
        self.modulation_type = None
        self.ber = None
        self.data_rate = None
        self.tx_pattern_file = None  # Optional GRASP .cut/.grd antenna patterns
        self.rx_pattern_file = None

        self.idx_found_sat1 = None
        self.idx_found_sat2 = None

        self.eirp = None
        self.metric = None
        self.cn0_required = 0

    def read_config(self, node):
        if node.find('SatelliteID1') is not None:
            self.sat_id1 = int(node.find('SatelliteID1').text)
        if node.find('SatelliteID2') is not None:
            self.sat_id2 = int(node.find('SatelliteID2').text)

        if node.find('CarrierFrequency') is not None:
            self.carrier_frequency = float(node.find('CarrierFrequency').text)
        if node.find('TransmitPowerW') is not None:
            self.transmit_power = float(node.find('TransmitPowerW').text)
        if node.find('TransmitLossesdB') is not None:
            self.transmit_losses = float(node.find('TransmitLossesdB').text)
        if node.find('TransmitGaindB') is not None:
            self.transmit_gain = float(node.find('TransmitGaindB').text)

        if node.find('ReceiveGaindB') is not None:
            self.receive_gain = float(node.find('ReceiveGaindB').text)
        if node.find('ReceiveLossesdB') is not None:
            self.receive_losses = float(node.find('ReceiveLossesdB').text)
        if node.find('ReceiveTempK') is not None:
            self.receive_temp = float(node.find('ReceiveTempK').text)

        if node.find('ModulationType') is not None:
            self.modulation_type = node.find('ModulationType').text
        if node.find('BitErrorRate') is not None:
            self.ber = float(node.find('BitErrorRate').text)
        if node.find('DataRateBitPerSec') is not None:
            self.data_rate = float(node.find('DataRateBitPerSec').text)
        if node.find('TransmitAntennaPatternFile') is not None:
            self.tx_pattern_file = node.find('TransmitAntennaPatternFile').text
        if node.find('ReceiveAntennaPatternFile') is not None:
            self.rx_pattern_file = node.find('ReceiveAntennaPatternFile').text

    def before_loop(self, sm):
        for idx_sat, satellite in enumerate(sm.satellites):
            if satellite.sat_id == self.sat_id1:
                self.idx_found_sat1 = idx_sat
            if satellite.sat_id == self.sat_id2:
                self.idx_found_sat2 = idx_sat
        # ISL antennas track each other: pattern files contribute their peak gain
        patterns = []
        if self.tx_pattern_file is not None:
            pattern = antenna.AntennaPattern.from_file(self.tx_pattern_file)
            self.transmit_gain = pattern.peak
            patterns.append(('Tx', pattern))
        if self.rx_pattern_file is not None:
            pattern = antenna.AntennaPattern.from_file(self.rx_pattern_file)
            self.receive_gain = pattern.peak
            patterns.append(('Rx', pattern))
        if patterns:
            antenna.plot_patterns(patterns, sm.output_path(self.type + '_antenna.png'))
        self.eirp = 10*log10(self.transmit_power) + self.transmit_gain - self.transmit_losses
        self.metric = np.zeros((sm.num_epoch, 5))
        if self.modulation_type is not None:
            self.cn0_required = misc_fn.comp_cn0_required(self.modulation_type, self.ber, self.data_rate)

    def in_loop(self, sm):
        if self.idx_found_sat2 in sm.satellites[self.idx_found_sat1].idx_sat_in_view:
            elevation = sm.sp2sp[self.idx_found_sat1][self.idx_found_sat2].elevation
            distance = sm.sp2sp[self.idx_found_sat1][self.idx_found_sat2].distance
            fsl = 20*log10(distance/1000) + 20*log10(self.carrier_frequency/1e9) + 92.45
            temp_sys = self.receive_temp + 20
            cn0 = self.eirp - fsl \
                  +self.receive_gain - self.receive_losses - \
                  K_BOLTZMANN - 10*np.log10(temp_sys)
            self.metric[sm.cnt_epoch,:] = [self.times_f_doy[sm.cnt_epoch], degrees(elevation),
                                           cn0, fsl, self.cn0_required]

    def after_loop(self, sm):
        self.metric = self.metric[~np.all(self.metric == 0, axis=1)]  # Clean up empty rows
        fig = plt.figure(figsize=(10, 6))
        plt.subplots_adjust(left=.1, right=.95, top=0.95, bottom=0.07)
        plt.plot(self.metric[:, 0], self.metric[:, 1], 'k.', label='Elevation')
        plt.plot(self.metric[:, 0], self.metric[:, 2], 'b.', label='CN0 Computed')
        plt.plot(self.metric[:, 0], self.metric[:, 3]-100, 'y.', label='Free space loss - 100dB')
        if self.modulation_type is not None:
            plt.plot(self.metric[:, 0], self.metric[:, 4], 'r.', label='CN0 Required')
        plt.xlabel('DOY[-]'); plt.ylabel('Elevation [deg], Power values [dB]')
        plt.legend(); plt.grid()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'elevation_deg', 'cn0_dbhz', 'fsl_db',
                            'cn0_required_dbhz'], self.metric)


class AnalysisComDoppler(AnalysisBase):

    def __init__(self):
        super().__init__()
        self.station_id = None
        self.carrier_frequency = None
        self.metric = None

    def read_config(self, node):
        if node.find('StationID') is not None:
            self.station_id = int(node.find('StationID').text)
        if node.find('CarrierFrequency') is not None:
            self.carrier_frequency = float(node.find('CarrierFrequency').text)

    def before_loop(self, sm):
        for idx_station, station in enumerate(sm.stations):
            if station.station_id == self.station_id:
                self.idx_found_station = idx_station
                break
        self.metric = np.zeros((sm.num_epoch, 3))

    def in_loop(self, sm):
        idx_station = self.idx_found_station
        for idx_sat in sm.stations[idx_station].idx_sat_in_view:
            velocity = sm.satellites[idx_sat].vel_ecf
            range = sm.gr2sp[idx_station][idx_sat].gr2sp_ecf
            range_rate = np.dot(velocity,range)/norm(range)
            doppler = self.carrier_frequency*range_rate/C_LIGHT
            elevation = sm.gr2sp[idx_station][idx_sat].elevation
            self.metric[sm.cnt_epoch,:] = [self.times_f_doy[sm.cnt_epoch], degrees(elevation), doppler]

    def after_loop(self, sm):
        self.metric = self.metric[~np.all(self.metric == 0, axis=1)]  # Clean up empty rows
        fig, ax1 = plt.subplots(figsize=(10, 6))
        plt.grid()
        ax1.set_ylabel('Doppler [kHz]')
        ax1.yaxis.label.set_color('red')
        ax1.tick_params(axis='y', colors='red')
        ax1.plot(self.metric[:, 0], self.metric[:, 2]/1000, 'r.', label='Doppler [kHz]')
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        ax2.set_ylabel('Elevation [deg]')
        ax2.yaxis.label.set_color('blue')
        ax2.tick_params(axis='y', colors='blue')
        ax2.plot(self.metric[:, 0], self.metric[:, 1], 'b.', label='Elevation [deg]')
        plt.xlabel('DOY[-]');
        plt.legend();
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'elevation_deg', 'doppler_hz'], self.metric)


class AnalysisComGr2SpBudgetInterference(AnalysisBase):

    def __init__(self):
        super().__init__()
        self.station_id = None
        self.transmitter_object = None
        self.carrier_frequency = None
        self.bandwith = None
        self.transmit_power = None
        self.transmit_losses = None
        self.transmit_gain = None
        self.transmit_gain_manual = None
        self.p_exceed = None
        self.include_rain = None
        self.include_gas = None
        self.include_scintillation = None
        self.include_clouds = None
        self.receive_gain = None
        self.receive_losses = None
        self.receive_temp = None
        self.modulation_type = None
        self.ber = None
        self.data_rate = None
        self.tx_pattern_file = None  # Optional GRASP .cut/.grd antenna patterns
        self.rx_pattern_file = None
        self.tx_pattern = None
        self.rx_pattern = None

        self.idx_found_station = None
        self.eirp = None

        self.metric = None
        self.cn0_required = 0

    def read_config(self, node):
        if node.find('GroundStationID') is not None:
            self.station_id = int(node.find('GroundStationID').text)
        if node.find('TransmitterObject') is not None:
            self.transmitter_object = node.find('TransmitterObject').text.lower()
        if node.find('CarrierFrequency') is not None:
            self.carrier_frequency = float(node.find('CarrierFrequency').text)
        if node.find('BandWidth') is not None:
            self.bandwidth = float(node.find('BandWidth').text)
        if node.find('TransmitAntennaPatternFile') is not None:
            self.tx_pattern_file = node.find('TransmitAntennaPatternFile').text
        if node.find('ReceiveAntennaPatternFile') is not None:
            self.rx_pattern_file = node.find('ReceiveAntennaPatternFile').text

        if node.find('TransmitPowerW') is not None:
            self.transmit_power = float(node.find('TransmitPowerW').text)
        if node.find('TransmitLossesdB') is not None:
            self.transmit_losses = float(node.find('TransmitLossesdB').text)
        if node.find('TransmitGaindB') is not None:
            self.transmit_gain = float(node.find('TransmitGaindB').text)
        if node.find('TransmitGainManualdB') is not None:
            self.transmit_gain_manual = np.array(list(ast.literal_eval(node.find('TransmitGainManualdB').text)))
        if node.find('TransmitAntennaDiameter') is not None:
            self.transmit_ant_dia = float(node.find('TransmitAntennaDiameter').text)

        if node.find('PExceedPerc') is not None:
            self.p_exceed = float(node.find('PExceedPerc').text)
        if node.find('IncludeRain') is not None:
            self.include_rain = misc_fn.str2bool(node.find('IncludeRain').text)
        if node.find('IncludeGas') is not None:
            self.include_gas = misc_fn.str2bool(node.find('IncludeGas').text)
        if node.find('IncludeScintillation') is not None:
            self.include_scintillation = misc_fn.str2bool(node.find('IncludeScintillation').text)
        if node.find('IncludeClouds') is not None:
            self.include_clouds = misc_fn.str2bool(node.find('IncludeClouds').text)

        if node.find('ReceiveGaindB') is not None:
            self.receive_gain = float(node.find('ReceiveGaindB').text)
        if node.find('ReceiveAntennaDiameter') is not None:
            self.receive_ant_dia = float(node.find('ReceiveAntennaDiameter').text)
        if node.find('ReceiveLossesdB') is not None:
            self.receive_losses = float(node.find('ReceiveLossesdB').text)
        if node.find('ReceiveTempK') is not None:
            self.receive_temp = float(node.find('ReceiveTempK').text)

        if node.find('ModulationType') is not None:
            self.modulation_type = node.find('ModulationType').text
        if node.find('BitErrorRate') is not None:
            self.ber = float(node.find('BitErrorRate').text)
        if node.find('DataRateBitPerSec') is not None:
            self.data_rate = float(node.find('DataRateBitPerSec').text)

    def before_loop(self, sm):
        for idx_station, station in enumerate(sm.stations):
            if station.station_id == self.station_id:
                self.idx_found_station = idx_station
                break
        if self.tx_pattern_file is not None:
            self.tx_pattern = antenna.AntennaPattern.from_file(self.tx_pattern_file)
        if self.rx_pattern_file is not None:
            self.rx_pattern = antenna.AntennaPattern.from_file(self.rx_pattern_file)
        if self.tx_pattern is not None or self.rx_pattern is not None:
            antenna.plot_patterns(
                [(label, p) for label, p in (('Tx', self.tx_pattern),
                                             ('Rx', self.rx_pattern)) if p is not None],
                sm.output_path(self.type + '_antenna.png'))
        self.metric = np.zeros((sm.num_epoch, 14))
        if self.modulation_type is not None:
            self.cn0_required = misc_fn.comp_cn0_required(self.modulation_type, self.ber, self.data_rate)

    def in_loop(self, sm):
        idx_station = self.idx_found_station
        if len(sm.stations[idx_station].idx_sat_in_view) > 1:

            elevation = sm.gr2sp[idx_station][0].elevation
            distance = sm.gr2sp[idx_station][0].distance

            # Nominal gains: both dishes point at each other, so pattern files
            # contribute their peak (boresight) gain on the nominal link
            tx_gain_nom, rx_gain_nom = self.transmit_gain, self.receive_gain
            if self.transmit_gain_manual is not None:
                tx_gain_nom = misc_fn.dish_pattern_manual(self.transmit_gain_manual, 0)
            if self.tx_pattern is not None:
                tx_gain_nom = self.tx_pattern.peak
            if self.rx_pattern is not None:
                rx_gain_nom = self.rx_pattern.peak

            # Nominal satellite C/N0
            self.eirp = 10 * log10(self.transmit_power) + tx_gain_nom - self.transmit_losses
            fsl = 20 * log10(distance / 1000) + 20 * log10(self.carrier_frequency / 1e9) + 92.45
            a_g, a_c, a_r, a_s, a_t = itur.atmospheric_attenuation_slant_path(
                degrees(sm.stations[idx_station].lla[0]),
                degrees(sm.stations[idx_station].lla[1]),
                self.carrier_frequency / 1e9 * itur.u.GHz,
                degrees(elevation), self.p_exceed,
                D=1.0 * itur.u.m,
                include_gas=self.include_gas,
                include_rain=self.include_rain,
                include_scintillation=self.include_scintillation,
                include_clouds=self.include_clouds,
                return_contributions=True)
            if self.transmitter_object == 'satellite':
                ant_temp = 10 + misc_fn.temp_brightness(self.carrier_frequency, elevation)
            else:  # station is the transmitter
                ant_temp = 290
            temp_sys = self.receive_temp + ant_temp
            cn0 = self.eirp - fsl \
                  - a_g.value - a_r.value - a_c.value - a_s.value \
                  + rx_gain_nom - self.receive_losses - \
                  K_BOLTZMANN - 10*np.log10(temp_sys)

            # Interference computation from second satellite
            u = sm.satellites[0].pos_ecf - sm.stations[idx_station].pos_ecf
            v = sm.satellites[1].pos_ecf - sm.stations[idx_station].pos_ecf
            off_boresight_angle = misc_fn.angle_two_vectors(u,v,np.linalg.norm(u),np.linalg.norm(v))
            C = self.eirp - fsl \
                  - a_g.value - a_r.value - a_c.value - a_s.value \
                  + rx_gain_nom - self.receive_losses
            C_fact = np.power(10, C/10) # all in factors since C0/(N0+I0)
            N0 = K_BOLTZMANN + 10*np.log10(temp_sys)
            N0_fact = np.power(10, N0/10)
            # Interferer gains: both antennas discriminated by the leader/
            # interferer separation angle, with the pattern files as drop-in
            # replacement of the analytic dish patterns
            if self.tx_pattern is not None:
                tx_gain = self.tx_pattern.gain(off_boresight_angle)
            elif self.transmit_gain_manual is not None:
                tx_gain = misc_fn.dish_pattern_manual(self.transmit_gain_manual,off_boresight_angle)
            else:
                tx_gain = misc_fn.dish_pattern(self.carrier_frequency,
                                           self.transmit_ant_dia,self.transmit_gain,off_boresight_angle)
            if self.rx_pattern is not None:
                rx_gain = self.rx_pattern.gain(off_boresight_angle)
            else:
                rx_gain = misc_fn.dish_pattern(self.carrier_frequency,
                                               self.receive_ant_dia,self.receive_gain,off_boresight_angle)
            I = 10 * log10(self.transmit_power) + tx_gain - self.transmit_losses + \
                - fsl - a_g.value - a_r.value - a_c.value - a_s.value + \
                + rx_gain - self.receive_losses
            I0 = I - 10*np.log10(self.bandwidth)
            I0_fact = np.power(10, I0/10)
            cni0 = 10*np.log10(C_fact / (N0_fact+I0_fact))

            self.metric[sm.cnt_epoch,:] = [self.times_f_doy[sm.cnt_epoch], degrees(elevation),
                                           cn0, cni0, a_g.value, a_r.value, a_c.value, a_s.value, a_t.value,
                                           degrees(off_boresight_angle), tx_gain, rx_gain,
                                           tx_gain_nom, rx_gain_nom]

    def after_loop(self, sm):
        self.metric = self.metric[~np.all(self.metric == 0, axis=1)]  # Clean up empty rows
        if self.metric.size == 0:
            ls.logger.error(f'No epoch had more than one satellite in view of GroundStationID '
                            f'{self.station_id}: the interference analysis needs the leading satellite '
                            f'(first in the space segment) and an interferer (second) visible '
                            f'simultaneously. Configure at least 2 satellites. No plot produced.')
            return
        fig = plt.figure(figsize=(10, 6))
        plt.subplots_adjust(left=.1, right=.95, top=0.95, bottom=0.07)
        time_list = self.metric[:, 0]
        plt.plot(time_list, self.metric[:, 1], label='Elevation leading satellite')
        plt.plot(time_list, self.metric[:, 2]-100, label='CN0 Nominal - 100')
        plt.plot(time_list, self.metric[:, 3]-100, label='CN0 with Interference - 100')
        plt.plot(time_list, self.metric[:, 8], label='Total Atmospheric Attenuation')
        plt.plot(time_list, self.metric[:, 9], label='Off boresight angle interferer')
        plt.plot(time_list, self.metric[:, 10], label='Tx Gain Interferer')
        plt.plot(time_list, self.metric[:, 11], label='Rx Gain Interferer')
        plt.plot(time_list, self.metric[:, 12], label='Tx Gain Nom')
        plt.plot(time_list, self.metric[:, 13], label='Rx Gain Nom')
        plt.gca().ticklabel_format(axis='both', style='plain', useOffset=False)
        plt.xlabel('Day Of Year DOY [-]'); plt.ylabel('Angles [deg], Power values [dB]')
        plt.legend(); plt.grid()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        fig = plt.figure(figsize=(10, 6))
        plt.subplots_adjust(left=.1, right=.95, top=0.95, bottom=0.07)
        plt.plot(time_list, self.metric[:, 2]-self.metric[:, 3], label='CN0 drop [dB]')
        plt.gca().ticklabel_format(axis='x', style='plain', useOffset=False)
        plt.gca().ticklabel_format(axis='y', style='sci', useOffset=False)
        plt.xlabel('Day Of Year DOY [-]'); plt.ylabel('Power values [dB]')
        plt.legend(); plt.grid()
        plt.savefig(sm.output_path(self.type + '_drop.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'elevation_deg', 'cn0_dbhz', 'cni0_dbhz', 'gas_att_db',
                            'rain_att_db', 'cloud_att_db', 'scint_att_db', 'total_att_db',
                            'off_boresight_deg', 'tx_gain_interferer_db', 'rx_gain_interferer_db',
                            'tx_gain_nominal_db', 'rx_gain_nominal_db'], self.metric)


class AnalysisComContactPlan(AnalysisBase):
    """Ground station contact plan: every station-satellite pass over the
    simulation window as an AOS/LOS table (CSV + human-readable text file),
    with the pass duration, maximum elevation, downlinkable data volume
    (optional DownlinkRateMbps) and a station-conflict flag when two passes
    overlap at the same station (one antenna cannot track two satellites).
    The plot is a pass timeline (Gantt) per station, overlapping passes edged
    in red. The log summarises contact time and volume per station and day."""

    def __init__(self):
        super().__init__()
        self.station_id = 0  # Optional: 0 = all stations
        self.constellation_id = 0  # Optional: 0 = all constellations
        self.min_duration = 0.0  # Passes shorter than this [s] are dropped
        self.downlink_rate = None  # Optional [Mbps] for the volume column
        self.elev_hist = None  # (num_station, num_sat, num_epoch) elevation [deg]

    def read_config(self, node):
        if node.find('GroundStationID') is not None:
            self.station_id = int(node.find('GroundStationID').text)
        if node.find('ConstellationID') is not None:
            self.constellation_id = int(node.find('ConstellationID').text)
        if node.find('MinDuration') is not None:
            self.min_duration = float(node.find('MinDuration').text)
        if node.find('DownlinkRateMbps') is not None:
            self.downlink_rate = float(node.find('DownlinkRateMbps').text)

    def before_loop(self, sm):
        self.elev_hist = np.full((sm.num_station, sm.num_sat, sm.num_epoch),
                                 np.nan, dtype=np.float32)

    def in_loop(self, sm):
        for idx_station, station in enumerate(sm.stations):
            if self.station_id > 0 and station.station_id != self.station_id:
                continue
            for idx_sat in station.idx_sat_in_view:
                if self.constellation_id > 0 and \
                        sm.satellites[idx_sat].constellation_id != self.constellation_id:
                    continue
                self.elev_hist[idx_station, idx_sat, sm.cnt_epoch] = \
                    degrees(sm.gr2sp[idx_station][idx_sat].elevation)

    def after_loop(self, sm):
        times_mjd = np.asarray(self.times_mjd)
        times_doy = np.asarray(self.times_f_doy)
        step = sm.time_step
        passes = []  # One dict per pass
        for idx_station, station in enumerate(sm.stations):
            for idx_sat, satellite in enumerate(sm.satellites):
                elevations = self.elev_hist[idx_station, idx_sat]
                idx_vis = np.flatnonzero(~np.isnan(elevations))
                if idx_vis.size == 0:
                    continue
                for run in np.split(idx_vis, np.flatnonzero(np.diff(idx_vis) > 1) + 1):
                    i0, i1 = run[0], run[-1]
                    duration = (i1 - i0 + 1) * step
                    if duration < self.min_duration:
                        continue
                    passes.append({
                        'idx_station': idx_station, 'station': station,
                        'satellite': satellite, 'i0': i0, 'i1': i1,
                        'aos_mjd': times_mjd[i0], 'los_mjd': times_mjd[i1],
                        'aos_doy': times_doy[i0], 'duration': duration,
                        'max_el': float(np.max(elevations[run])),
                        'volume': duration * self.downlink_rate / 1000.0
                                  if self.downlink_rate else 0.0,
                        'overlap': 0})
        passes.sort(key=lambda p: (p['idx_station'], p['aos_mjd']))
        # Station conflicts: passes at the same station that overlap in time
        for i, p in enumerate(passes):
            for q in passes[i + 1:]:
                if q['idx_station'] != p['idx_station']:
                    break  # Sorted by station: no further pairs at this station
                if q['i0'] <= p['i1']:
                    p['overlap'] += 1
                    q['overlap'] += 1
        if not passes:
            ls.logger.warning(f'{self.type}: no passes found in the simulation window. '
                              f'No plot produced.')
            self.write_csv(sm, ['station_id', 'sat_id', 'aos_mjd', 'los_mjd', 'aos_doy',
                                'duration_s', 'max_elevation_deg', 'overlap',
                                'volume_gbit'], [])
            return

        window_days = sm.num_epoch * step / 86400.0
        for idx_station, station in enumerate(sm.stations):
            own = [p for p in passes if p['idx_station'] == idx_station]
            if not own:
                continue
            total = sum(p['duration'] for p in own)
            ls.logger.info(f"{self.type}: {station.station_name}: {len(own)} passes, "
                           f"{total / 60.0 / window_days:.1f} min/day contact "
                           f"(mean {np.mean([p['duration'] for p in own]):.0f} s, max "
                           f"{max(p['duration'] for p in own):.0f} s), "
                           f"{sum(p['volume'] for p in own) / window_days:.1f} Gbit/day, "
                           f"{sum(1 for p in own if p['overlap'])} overlapping")

        # Human-readable pass table next to the CSV
        with open(sm.output_path(self.type + '.txt'), 'w') as f:
            f.write(f"{'Station':<14}{'Satellite':>10}  {'AOS (UTC)':<20}"
                    f"{'LOS (UTC)':<20}{'Dur [s]':>8}{'MaxEl':>7}{'Gbit':>8}{'Ovl':>5}\n")
            for p in sorted(passes, key=lambda p: p['aos_mjd']):
                f.write(f"{p['station'].station_name:<14}{p['satellite'].sat_id:>10}  "
                        f"{Time(p['aos_mjd'], format='mjd').iso[:19]:<20}"
                        f"{Time(p['los_mjd'], format='mjd').iso[:19]:<20}"
                        f"{p['duration']:>8.0f}{p['max_el']:>7.1f}"
                        f"{p['volume']:>8.2f}{p['overlap']:>5d}\n")
        ls.logger.info(f'Written pass table {self.type}.txt ({len(passes)} passes)')

        # Pass timeline: one lane per station with a sub-lane per satellite,
        # so simultaneous passes (station conflicts) are visible side by side
        sat_ids = sorted({p['satellite'].sat_id for p in passes})
        colors = plt.cm.tab20(np.linspace(0, 0.95, max(len(sat_ids), 2)))
        color_of = {sat_id: colors[i] for i, sat_id in enumerate(sat_ids)}
        lane_of = {sat_id: i for i, sat_id in enumerate(sat_ids)}
        lane_height = 0.7 / len(sat_ids)
        fig, ax = plt.subplots(figsize=(12, 1.2 + 1.1 * sm.num_station))
        for p in passes:
            y0 = p['idx_station'] - 0.35 + lane_of[p['satellite'].sat_id] * lane_height
            ax.broken_barh([(p['aos_doy'], p['duration'] / 86400.0)],
                           (y0, lane_height * 0.9),
                           facecolors=[color_of[p['satellite'].sat_id]],
                           edgecolors='red' if p['overlap'] else 'black',
                           linewidth=1.2 if p['overlap'] else 0.3)
        ax.set_yticks(range(sm.num_station))
        ax.set_yticklabels([station.station_name for station in sm.stations])
        ax.set_ylim(-0.6, sm.num_station - 0.4)
        ax.set_xlabel('DOY [-]')
        ax.grid(True, axis='x')
        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=color_of[sat_id])
                   for sat_id in sat_ids[:12]]
        ax.legend(handles, [f'Sat {sat_id}' for sat_id in sat_ids[:12]],
                  fontsize=8, loc='upper right')
        ax.set_title(f'Contact plan: {len(passes)} passes, red edge = station conflict')
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['station_id', 'sat_id', 'aos_mjd', 'los_mjd', 'aos_doy',
                            'duration_s', 'max_elevation_deg', 'overlap', 'volume_gbit'],
                       [[p['station'].station_id, p['satellite'].sat_id, p['aos_mjd'],
                         p['los_mjd'], p['aos_doy'], p['duration'], p['max_el'],
                         p['overlap'], p['volume']] for p in passes])


class AnalysisComPfd(AnalysisBase):
    """Power flux density produced at the ground by the satellite emissions,
    versus elevation, against an ITU-R Article 21 style limit mask. Per epoch
    the PFD at the ground station is EIRP spectral density (transmit power
    spread over BandWidth, taken in the ReferenceBandwidth) plus the satellite
    antenna gain at the epoch's off-nadir angle (nadir-pointed pattern file or
    a fixed gain) minus the spreading loss 10log10(4 pi d^2). The limit mask
    follows the standard Article 21 shape: PfdLimit up to 5 deg elevation,
    +0.5 dB/deg between 5 and 25 deg, PfdLimit+10 above 25 deg. Atmospheric
    attenuation is not credited (conservative, as in the Radio Regulations)."""

    def __init__(self):
        super().__init__()
        self.station_id = None  # Optional: default first station
        self.carrier_frequency = None
        self.transmit_power = None
        self.transmit_losses = 0.0
        self.bandwidth = None  # Occupied bandwidth [Hz] the power is spread over
        self.ref_bandwidth = 4000.0  # ITU reference bandwidth [Hz] (4 kHz < 15 GHz)
        self.transmit_gain = None
        self.tx_pattern_file = None  # Optional GRASP .cut/.grd pattern (nadir pointed)
        self.tx_pattern = None
        self.pfd_limit = -150.0  # Limit at <= 5 deg elevation [dB(W/m2) in ref BW]
        self.idx_found_station = 0
        self.metric = None  # (num_epoch, 5): doy, el, pfd, limit, margin

    def read_config(self, node):
        if node.find('GroundStationID') is not None:
            self.station_id = int(node.find('GroundStationID').text)
        if node.find('CarrierFrequency') is not None:
            self.carrier_frequency = float(node.find('CarrierFrequency').text)
        if node.find('TransmitPowerW') is not None:
            self.transmit_power = float(node.find('TransmitPowerW').text)
        if node.find('TransmitLossesdB') is not None:
            self.transmit_losses = float(node.find('TransmitLossesdB').text)
        if node.find('BandWidth') is not None:
            self.bandwidth = float(node.find('BandWidth').text)
        if node.find('ReferenceBandwidth') is not None:
            self.ref_bandwidth = float(node.find('ReferenceBandwidth').text)
        if node.find('TransmitGaindB') is not None:
            self.transmit_gain = float(node.find('TransmitGaindB').text)
        if node.find('TransmitAntennaPatternFile') is not None:
            self.tx_pattern_file = node.find('TransmitAntennaPatternFile').text
        if node.find('PfdLimit') is not None:
            self.pfd_limit = float(node.find('PfdLimit').text)

    def _limit(self, elevation_deg):
        """Article 21 style mask: flat to 5 deg, +0.5 dB/deg to 25 deg, then flat."""
        return self.pfd_limit + 0.5 * np.clip(elevation_deg - 5.0, 0.0, 20.0)

    def before_loop(self, sm):
        if self.station_id is not None:
            for idx_station, station in enumerate(sm.stations):
                if station.station_id == self.station_id:
                    self.idx_found_station = idx_station
                    break
        if self.tx_pattern_file is not None:
            self.tx_pattern = antenna.AntennaPattern.from_file(self.tx_pattern_file)
            antenna.plot_patterns([('Tx', self.tx_pattern)],
                                  sm.output_path(self.type + '_antenna.png'))
        self.metric = np.zeros((sm.num_epoch, 5))

    def in_loop(self, sm):
        idx_station = self.idx_found_station
        station = sm.stations[idx_station]
        worst = None  # Keep the in-view satellite with the smallest margin
        for idx_sat in station.idx_sat_in_view:
            link = sm.gr2sp[idx_station][idx_sat]
            elevation_deg = degrees(link.elevation)
            off_nadir = _off_nadir_angle(sm.satellites[idx_sat], station)
            gain = self.tx_pattern.gain(off_nadir) if self.tx_pattern is not None \
                else self.transmit_gain
            # Power in the reference bandwidth, spread over the ground distance sphere
            pfd = 10 * log10(self.transmit_power * self.ref_bandwidth / self.bandwidth) \
                - self.transmit_losses + gain - 10 * log10(4 * np.pi * link.distance ** 2)
            margin = self._limit(elevation_deg) - pfd
            if worst is None or margin < worst[4]:
                worst = [self.times_f_doy[sm.cnt_epoch], elevation_deg, pfd,
                         self._limit(elevation_deg), margin]
        if worst is not None:
            self.metric[sm.cnt_epoch, :] = worst

    def after_loop(self, sm):
        self.metric = self.metric[~np.all(self.metric == 0, axis=1)]  # Clean up empty rows
        self.write_csv(sm, ['doy', 'elevation_deg', 'pfd_dbwm2', 'limit_dbwm2',
                            'margin_db'], self.metric)
        if self.metric.size == 0:
            ls.logger.warning(f'{self.type}: no satellite was in view of the station. '
                              f'No plot produced.')
            return
        idx_worst = np.argmin(self.metric[:, 4])
        ls.logger.info(f'{self.type}: worst margin {self.metric[idx_worst, 4]:.1f} dB at '
                       f'{self.metric[idx_worst, 1]:.1f} deg elevation (PFD '
                       f'{self.metric[idx_worst, 2]:.1f} vs limit '
                       f'{self.metric[idx_worst, 3]:.1f} dB(W/m2) in '
                       f'{self.ref_bandwidth / 1000:.0f} kHz)')

        fig, ax = plt.subplots(figsize=(10, 6))
        sc = ax.scatter(self.metric[:, 1], self.metric[:, 2], s=8, c='b',
                        label='PFD at the station')
        el_line = np.linspace(0.0, 90.0, 181)
        ax.plot(el_line, self._limit(el_line), 'r-', linewidth=1.2,
                label='PFD limit mask')
        ax.set_xlabel('Elevation [deg]')
        ax.set_ylabel(f'PFD [dB(W/m$^2$) in {self.ref_bandwidth / 1000:.0f} kHz]')
        ax.grid(True)
        ax.legend()
        ax.set_title(f'Power flux density at {sm.stations[self.idx_found_station].station_name}: '
                     f'worst margin {self.metric[idx_worst, 4]:.1f} dB')
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()
