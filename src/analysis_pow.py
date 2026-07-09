import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import get_sun
from astropy.time import Time
import astropy.units as u
from analysis import AnalysisBase
import logging_svs as ls
from constants import R_EARTH
import misc_fn

class AnalysisPowDepthDischarge(AnalysisBase):
    def __init__(self):
        super().__init__()
        self.battery_capacity_wh = 0
        self.initial_soc = 1.0
        self.panel_area = 0
        self.efficiency = 0
        self.p_bus_w = 0
        self.p_payload_w = 0
        self.metric = None # Stores [Time, SoC, P_gen, Eclipse_Status]

    def read_config(self, node):
        self.battery_capacity_wh = float(node.find('BatteryCapacityWh').text)
        self.initial_soc = float(node.find('InitialSoC').text)
        self.panel_area = float(node.find('SolarPanelArea').text)
        self.efficiency = float(node.find('PanelEfficiency').text)
        self.p_bus_w = float(node.find('BasePowerDrawW').text)
        self.p_payload_w = float(node.find('InstrumentPowerDrawW').text)
        # Read the latitude limit from XML
        if node.find('PayloadLatitudeLimit') is not None:
            self.lat_limit = float(node.find('PayloadLatitudeLimit').text)
        else:
            self.lat_limit = 90.0 # Default to "always on" if not specified

    def before_loop(self, sm):
        # Initialize metric array: [Time, SoC, Generation, PowerDraw]
        self.metric = np.zeros((sm.num_epoch, 4))
        self.current_soc = self.initial_soc
        ls.logger.info("Power Analysis Initialized")

    def in_loop(self, sm):
        # 1. Determine Sun Visibility (Eclipse Check)
        t = Time(sm.time_mjd, format='mjd')
        sm.satellites[0].det_lla() # update LLA
        # Calculate the local Earth radius based on current latitude
        local_r = misc_fn.earth_radius_lat(sm.satellites[0].lla[0])
        # 1. Normalize the direction to the Sun
        sun_pos_eci = get_sun(t).cartesian.xyz.to(u.m).value
        sun_dist = np.linalg.norm(sun_pos_eci)
        sun_dir = sun_pos_eci / sun_dist

        # 2. Project satellite position onto the Sun-direction vector
        # This determines how far 'forward' or 'backward' the satellite is relative to Earth
        sat_pos = sm.satellites[0].pos_eci
        projection = np.dot(sat_pos, sun_dir)

        # 3. Calculate minimum distance from Earth Center (0,0,0) to the Sat-Sun line
        # Use the Pythagorean theorem: dist^2 + proj^2 = sat_pos_norm^2
        dist_to_sun_line = np.sqrt(np.linalg.norm(sat_pos)**2 - projection**2)

        # 4. Refined Eclipse Condition
        # A satellite is in eclipse ONLY if:
        # A) It is on the night side of Earth (projection < 0)
        # B) The line of sight to the Sun is blocked by Earth's radius
        in_eclipse = (projection < 0) and (dist_to_sun_line < local_r)

        # 2. Power Generation
        solar_constant = 1361.0 # W/m^2
        p_gen = 0 if in_eclipse else (solar_constant * self.panel_area * self.efficiency)

        # 3. Power Consumption
        # Determine current latitude in degrees

        is_active = abs(np.degrees(sm.satellites[0].lla[0])) <= self.lat_limit
        
        p_draw = self.p_bus_w + (self.p_payload_w if is_active else 0)

        # 4. Update SoC (Wh)
        delta_hours = sm.time_step / 3600.0
        energy_balance_wh = (p_gen - p_draw) * delta_hours
        
        # Convert capacity to Wh and update
        self.current_soc += energy_balance_wh / self.battery_capacity_wh
        self.current_soc = np.clip(self.current_soc, 0, 1) # Keep between 0% and 100%

        self.metric[sm.cnt_epoch, :] = [self.times_f_doy[-1], self.current_soc, p_gen, p_draw]

    def after_loop(self, sm):
        dod = (1.0 - self.metric[:, 1]) * 100
        
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.set_xlabel('Day of Year (DOY)')
        ax1.set_ylabel('Depth of Discharge DOD (%)', color='tab:red')
        ax1.plot(self.metric[:, 0], dod, color='tab:red', label='DoD %')
        ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.set_ylabel('Power Generated-Blue Draw-Green (W)', color='tab:blue')
        ax2.step(self.metric[:, 0], self.metric[:, 2], color='tab:blue', label='P_Gen', where='post')
        ax2.step(self.metric[:, 0], self.metric[:, 3], color='tab:green', label='P_Draw', where='post')
        
        fig.tight_layout()
        plt.savefig('../output/pow_depth_discharge.png')
        plt.show()

class AnalysisPowEclipseDuration(AnalysisBase):
    def __init__(self):
        super().__init__()
        self.eclipse_durations = [] # List of (Time, Duration_minutes)
        self.in_eclipse_prev = False
        self.eclipse_start_time = 0
        self.local_r = 0

    def read_config(self, node):
        # No specific config needed for this, but could add altitude-specific masks if desired
        pass

    def before_loop(self, sm):
        self.eclipse_durations = []
        self.in_eclipse_prev = False
        ls.logger.info("Eclipse Duration Analysis Initialized")

    def in_loop(self, sm):
        t = Time(sm.time_mjd, format='mjd')
        sat = sm.satellites[0]
        sat.det_lla()
        local_r = misc_fn.earth_radius_lat(sat.lla[0])

        # Vector Projection Eclipse Logic
        sun_pos_eci = get_sun(t).cartesian.xyz.to(u.m).value
        sun_dir = sun_pos_eci / np.linalg.norm(sun_pos_eci)
        sat_pos = sat.pos_eci
        projection = np.dot(sat_pos, sun_dir)
        dist_to_sun_line = np.sqrt(np.linalg.norm(sat_pos)**2 - projection**2)

        # Condition for eclipse
        in_eclipse_now = (projection < 0) and (dist_to_sun_line < local_r)

        # Detect transition from sunlight to eclipse (Entry)
        if in_eclipse_now and not self.in_eclipse_prev:
            self.eclipse_start_time = sm.time_mjd
        
        # Detect transition from eclipse to sunlight (Exit)
        elif not in_eclipse_now and self.in_eclipse_prev:
            duration_days = sm.time_mjd - self.eclipse_start_time
            duration_minutes = duration_days * 24 * 60
            self.eclipse_durations.append([self.times_f_doy[-1], duration_minutes])

        self.in_eclipse_prev = in_eclipse_now

    def after_loop(self, sm):
        if not self.eclipse_durations:
            ls.logger.warning("No eclipse events detected during the simulation. No plot produced.")
            return

        data = np.array(self.eclipse_durations)
        
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(data[:, 0], data[:, 1], 'k.', markersize=2, label='Eclipse Duration')
        ax.set_xlabel('Day of Year (DOY)')
        ax.set_ylabel('Eclipse Duration [minutes]')
        ax.set_title(f'Eclipse Duration per Orbit')
        ax.grid(True, which='both', linestyle='--', alpha=0.5)
        
        plt.savefig('../output/pow_eclipse_duration.png')
        plt.show()