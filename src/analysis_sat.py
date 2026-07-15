"""
Satellite platform (subsystem) analyses:
- sat_thermal: single-node thermal balance over the orbit
- sat_aocs: AOCS disturbance torques and momentum buildup
- sat_battery_depth_discharge / sat_eclipse_duration: electrical power
- sat_data_storage / sat_data_latency: data handling (SSR fill state, latency)

All work with any orbit propagator; sat_aocs samples the HPOP atmosphere
model when available and falls back to a built-in exponential atmosphere.
"""
import numpy as np
import matplotlib.pyplot as plt
from math import degrees, radians
from astropy.coordinates import get_sun
from astropy.time import Time
import astropy.units as u

# Project modules
from constants import R_EARTH, GM_EARTH, SOLAR_FLUX, SOLAR_PRESSURE, ALPHA
from analysis import AnalysisBase
import misc_fn
import logging_svs as ls

STEFAN_BOLTZMANN = 5.670374419e-8  # [W/m2/K4]
EARTH_IR_FLUX = 237.0  # Mean Earth infrared emission at the surface [W/m2]

# Centred dipole model of the geomagnetic field (SMAD): field constant at the
# surface on the magnetic equator and the geomagnetic north pole (IGRF-13 2020)
B0_EQUATOR = 3.12e-5  # [T]
DIPOLE_LAT = radians(80.7)  # Geomagnetic north pole latitude
DIPOLE_LON = radians(-72.7)  # Geomagnetic north pole longitude

# Piecewise-exponential atmosphere (Vallado, Fundamentals of Astrodynamics,
# Table 8-4): base altitude [km], base density [kg/m3], scale height [km]
_ATMOSPHERE_TABLE = np.array([
    [0, 1.225, 7.249], [25, 3.899e-2, 6.349], [30, 1.774e-2, 6.682],
    [40, 3.972e-3, 7.554], [50, 1.057e-3, 8.382], [60, 3.206e-4, 7.714],
    [70, 8.770e-5, 6.549], [80, 1.905e-5, 5.799], [90, 3.396e-6, 5.382],
    [100, 5.297e-7, 5.877], [110, 9.661e-8, 7.263], [120, 2.438e-8, 9.473],
    [130, 8.484e-9, 12.636], [140, 3.845e-9, 16.149], [150, 2.070e-9, 22.523],
    [180, 5.464e-10, 29.740], [200, 2.789e-10, 37.105], [250, 7.248e-11, 45.546],
    [300, 2.418e-11, 53.628], [350, 9.518e-12, 53.298], [400, 3.725e-12, 58.515],
    [450, 1.585e-12, 60.828], [500, 6.967e-13, 63.822], [600, 1.454e-13, 71.835],
    [700, 3.614e-14, 88.667], [800, 1.170e-14, 124.64], [900, 5.245e-15, 181.05],
    [1000, 3.019e-15, 268.00]])


def air_density_exponential(altitude_m):
    """Atmospheric density [kg/m3] from the piecewise-exponential model
    (Vallado Table 8-4); crude compared to NRLMSISE00 but propagator-free."""
    h_km = max(altitude_m / 1000.0, 0.0)
    idx = int(np.searchsorted(_ATMOSPHERE_TABLE[:, 0], h_km, side='right')) - 1
    h0, rho0, scale_h = _ATMOSPHERE_TABLE[idx]
    return rho0 * np.exp(-(h_km - h0) / scale_h)


def sun_direction_and_eclipse(satellite, time_mjd):
    """Unit vector to the Sun in ECI and the cylindrical-shadow eclipse flag
    for the satellite (same shadow test as the power analyses)."""
    t = Time(time_mjd, format='mjd')
    sun_pos_eci = get_sun(t).cartesian.xyz.to(u.m).value
    sun_dir = sun_pos_eci / np.linalg.norm(sun_pos_eci)
    satellite.det_lla()
    local_r = misc_fn.earth_radius_lat(satellite.lla[0])
    sat_pos = satellite.pos_eci
    projection = np.dot(sat_pos, sun_dir)
    dist_to_sun_line = np.sqrt(max(np.linalg.norm(sat_pos) ** 2 - projection ** 2, 0.0))
    in_eclipse = (projection < 0) and (dist_to_sun_line < local_r)
    return sun_dir, in_eclipse


class _SelectSatellite:
    """Optional ConstellationID/SatelliteID selection shared by the sat_
    analyses: the first matching satellite is analysed."""

    def read_selection(self, node):
        self.constellation_id = int(node.find('ConstellationID').text) \
            if node.find('ConstellationID') is not None else 0
        self.satellite_id = int(node.find('SatelliteID').text) \
            if node.find('SatelliteID') is not None else 0

    def find_satellite(self, sm):
        for idx_sat, satellite in enumerate(sm.satellites):
            if self.constellation_id > 0 and satellite.constellation_id != self.constellation_id:
                continue
            if self.satellite_id > 0 and satellite.sat_id != self.satellite_id:
                continue
            return idx_sat
        ls.logger.error(f'No satellite matched ConstellationID {self.constellation_id} / '
                        f'SatelliteID {self.satellite_id}. Analysis skipped.')
        return None


class AnalysisSatThermal(AnalysisBase, _SelectSatellite):
    """Single-node spacecraft thermal balance over the orbit: direct solar
    flux (with eclipses), Earth albedo, Earth infrared and internal
    dissipation against the radiated heat, integrated per time step to a
    temperature history. The classic hot-case/cold-case equilibrium
    temperatures are logged alongside."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection (first match is analysed)
        self.idx_found_satellite = None
        self.surface_area = None  # Total radiating surface area [m2]
        self.cross_section_sun = None  # Projected area towards the Sun [m2]
        self.cross_section_earth = None  # Projected area towards the Earth [m2]
        self.absorptivity = 0.3  # Solar absorptivity alpha
        self.emissivity = 0.8  # Infrared emissivity epsilon
        self.internal_power = 0.0  # Dissipated electrical power [W]
        self.heat_capacity = None  # Thermal capacitance [J/K]
        self.initial_temp = None  # [K]; default: equilibrium at the first epoch
        self.temperature = None  # Current node temperature [K]
        self.metric = None  # (num_epoch, 7): T [K], heat flows [W], eclipse flag
        self.enabled = True

    def read_config(self, node):
        self.read_selection(node)
        self.surface_area = float(node.find('SurfaceArea').text)
        self.heat_capacity = float(node.find('HeatCapacity').text)
        # A sphere presents a quarter of its surface to any direction: the
        # default cross sections, override for other shapes
        self.cross_section_sun = float(node.find('CrossSectionSun').text) \
            if node.find('CrossSectionSun') is not None else self.surface_area / 4.0
        self.cross_section_earth = float(node.find('CrossSectionEarth').text) \
            if node.find('CrossSectionEarth') is not None else self.surface_area / 4.0
        if node.find('Absorptivity') is not None:
            self.absorptivity = float(node.find('Absorptivity').text)
        if node.find('Emissivity') is not None:
            self.emissivity = float(node.find('Emissivity').text)
        if node.find('InternalPowerW') is not None:
            self.internal_power = float(node.find('InternalPowerW').text)
        if node.find('InitialTemperature') is not None:  # degC in the config
            self.initial_temp = float(node.find('InitialTemperature').text) + 273.15

    def before_loop(self, sm):
        self.idx_found_satellite = self.find_satellite(sm)
        self.enabled = self.idx_found_satellite is not None
        self.metric = np.full((sm.num_epoch, 7), np.nan)
        self.temperature = self.initial_temp  # None: set at the first epoch
        ls.logger.info(f'Thermal analysis: area {self.surface_area} m2, alpha '
                       f'{self.absorptivity}, epsilon {self.emissivity}, C '
                       f'{self.heat_capacity} J/K, Q_int {self.internal_power} W')

    def _equilibrium_temp(self, q_in):
        return (q_in / (self.emissivity * STEFAN_BOLTZMANN * self.surface_area)) ** 0.25

    def in_loop(self, sm):
        if not self.enabled:
            return
        satellite = sm.satellites[self.idx_found_satellite]
        sun_dir, in_eclipse = sun_direction_and_eclipse(satellite, sm.time_mjd)

        r = np.linalg.norm(satellite.pos_eci)
        view_factor = (R_EARTH / r) ** 2  # Nadir-facing plate view factor of Earth
        # Direct solar flux, zero in the Earth shadow
        q_sun = 0.0 if in_eclipse else \
            SOLAR_FLUX * self.cross_section_sun * self.absorptivity
        # Albedo: reflected sunlight from the sunlit Earth below the satellite
        cos_sun_zenith = max(np.dot(satellite.pos_eci / r, sun_dir), 0.0)
        q_albedo = SOLAR_FLUX * ALPHA * self.cross_section_earth * \
            self.absorptivity * view_factor * cos_sun_zenith
        # Earth infrared, always on
        q_ir = EARTH_IR_FLUX * self.cross_section_earth * self.emissivity * view_factor
        q_in = q_sun + q_albedo + q_ir + self.internal_power

        if self.temperature is None:  # Start in equilibrium with the first epoch
            self.temperature = self._equilibrium_temp(q_in)
        # Explicit Euler with bounded substeps (the radiative time constant is
        # hours for realistic C/A, so 60 s substeps are comfortably stable)
        num_sub = max(1, int(np.ceil(sm.time_step / 60.0)))
        dt = sm.time_step / num_sub
        for _ in range(num_sub):
            q_rad = self.emissivity * STEFAN_BOLTZMANN * self.surface_area * \
                self.temperature ** 4
            self.temperature += (q_in - q_rad) * dt / self.heat_capacity
        q_rad = self.emissivity * STEFAN_BOLTZMANN * self.surface_area * self.temperature ** 4
        self.metric[sm.cnt_epoch] = [self.temperature, q_sun, q_albedo, q_ir,
                                     self.internal_power, q_rad, float(in_eclipse)]

    def after_loop(self, sm):
        if not self.enabled:
            return
        times = np.asarray(self.times_f_doy)
        temp_c = self.metric[:, 0] - 273.15
        q_in = self.metric[:, 1:5].sum(axis=1)
        # Hot case: full sun plus albedo/IR; cold case: eclipse (IR + internal only)
        t_eq_hot = self._equilibrium_temp(np.nanmax(q_in)) - 273.15
        t_eq_cold = self._equilibrium_temp(np.nanmin(q_in)) - 273.15
        ls.logger.info(f'{self.type}: temperature {np.nanmin(temp_c):.1f} .. '
                       f'{np.nanmax(temp_c):.1f} degC, equilibrium hot case '
                       f'{t_eq_hot:.1f} degC, cold case {t_eq_cold:.1f} degC')

        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(times, temp_c, 'r-', linewidth=1.2, label='Temperature')
        ax1.set_xlabel('DOY [-]')
        ax1.set_ylabel('Temperature [degC]', color='red')
        ax1.tick_params(axis='y', colors='red')
        ax1.grid(True)
        ax2 = ax1.twinx()
        ax2.plot(times, self.metric[:, 1], 'y-', linewidth=0.8, label='Solar')
        ax2.plot(times, self.metric[:, 2], 'c-', linewidth=0.8, label='Albedo')
        ax2.plot(times, self.metric[:, 3], 'g-', linewidth=0.8, label='Earth IR')
        ax2.plot(times, self.metric[:, 5], 'k--', linewidth=0.8, label='Radiated')
        ax2.set_ylabel('Heat flow [W]')
        ax1.legend(loc='upper left', fontsize=8)
        ax2.legend(loc='upper right', fontsize=8)
        sat_id = sm.satellites[self.idx_found_satellite].sat_id
        fig.suptitle(f'Single-node thermal balance, satellite {sat_id} '
                     f'(hot eq. {t_eq_hot:.0f} degC, cold eq. {t_eq_cold:.0f} degC)')
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'temperature_k', 'q_sun_w', 'q_albedo_w',
                            'q_earth_ir_w', 'q_internal_w', 'q_radiated_w', 'eclipse'],
                       np.column_stack([times, self.metric]))


class AnalysisSatAocs(AnalysisBase, _SelectSatellite):
    """AOCS disturbance torques over the orbit with the standard worst-case
    models (SMAD): gravity gradient, aerodynamic, solar radiation pressure
    (zero in eclipse) and magnetic residual-dipole torque, plus the momentum
    buildup (integral of the total torque) for reaction wheel sizing. The
    atmospheric density comes from the HPOP DragModel when that propagator is
    active, and from a built-in exponential atmosphere otherwise."""

    def __init__(self):
        super().__init__()
        self.constellation_id = 0  # Optional selection
        self.satellite_id = 0  # Optional selection (first match is analysed)
        self.idx_found_satellite = None
        self.inertia = None  # Principal moments of inertia [kg m2]
        self.pointing_offset = radians(1.0)  # Attitude deviation for gravity gradient
        self.residual_dipole = 1.0  # Residual magnetic dipole [A m2]
        self.drag_area = None  # [m2]; default satellite FrontalArea
        self.drag_cd = 2.2
        self.srp_area = None  # [m2]; default DragArea
        self.reflectivity = 0.6  # SRP reflectance factor q
        self.cop_offset = 0.1  # Centre of pressure - centre of mass offset [m]
        self.wheel_momentum = None  # Optional wheel capacity [N m s] for the plot
        self.dipole_axis_ecf = np.array([  # Unit vector to the geomagnetic north pole
            np.cos(DIPOLE_LAT) * np.cos(DIPOLE_LON),
            np.cos(DIPOLE_LAT) * np.sin(DIPOLE_LON),
            np.sin(DIPOLE_LAT)])
        self.momentum = 0.0  # Accumulated angular momentum [N m s]
        self.metric = None  # (num_epoch, 6): torques [N m] + cumulative momentum
        self.enabled = True

    def read_config(self, node):
        self.read_selection(node)
        self.inertia = np.array([float(node.find('InertiaXX').text),
                                 float(node.find('InertiaYY').text),
                                 float(node.find('InertiaZZ').text)])
        if node.find('MaxPointingOffset') is not None:
            self.pointing_offset = radians(float(node.find('MaxPointingOffset').text))
        if node.find('ResidualDipole') is not None:
            self.residual_dipole = float(node.find('ResidualDipole').text)
        if node.find('DragArea') is not None:
            self.drag_area = float(node.find('DragArea').text)
        if node.find('DragCoefficient') is not None:
            self.drag_cd = float(node.find('DragCoefficient').text)
        if node.find('SrpArea') is not None:
            self.srp_area = float(node.find('SrpArea').text)
        if node.find('Reflectivity') is not None:
            self.reflectivity = float(node.find('Reflectivity').text)
        if node.find('CopOffset') is not None:
            self.cop_offset = float(node.find('CopOffset').text)
        if node.find('WheelMomentum') is not None:
            self.wheel_momentum = float(node.find('WheelMomentum').text)

    def before_loop(self, sm):
        self.idx_found_satellite = self.find_satellite(sm)
        self.enabled = self.idx_found_satellite is not None
        if self.enabled and self.drag_area is None:
            frontal = sm.satellites[self.idx_found_satellite].frontal_area
            self.drag_area = frontal if frontal is not None else 1.0
            ls.logger.info(f'{self.type}: no <DragArea> given, using '
                           f'{self.drag_area} m2 (constellation FrontalArea or 1)')
        if self.srp_area is None and self.drag_area is not None:
            self.srp_area = self.drag_area
        self.metric = np.full((sm.num_epoch, 6), np.nan)
        self.momentum = 0.0
        ls.logger.info(f'AOCS analysis: inertia {self.inertia} kg m2, dipole '
                       f'{self.residual_dipole} A m2, drag area {self.drag_area} m2, '
                       f'CoP offset {self.cop_offset} m')

    def in_loop(self, sm):
        if not self.enabled:
            return
        satellite = sm.satellites[self.idx_found_satellite]
        sun_dir, in_eclipse = sun_direction_and_eclipse(satellite, sm.time_mjd)
        r = np.linalg.norm(satellite.pos_eci)
        altitude = np.linalg.norm(satellite.pos_ecf) - \
            misc_fn.earth_radius_lat(satellite.lla[0])

        # Gravity gradient torque at the worst-case attitude deviation
        t_gg = 1.5 * GM_EARTH / r ** 3 * (np.max(self.inertia) - np.min(self.inertia)) \
            * abs(np.sin(2.0 * self.pointing_offset))
        # Aerodynamic torque: dynamic pressure on the drag area with the CoP lever
        if sm.hpop is not None:
            density = sm.hpop.air_density(satellite.pos_ecf, sm.time_mjd)
        else:
            density = air_density_exponential(altitude)
        velocity = np.linalg.norm(satellite.vel_eci)
        t_aero = 0.5 * density * velocity ** 2 * self.drag_cd * self.drag_area * \
            self.cop_offset
        # Solar radiation pressure torque, zero in eclipse
        t_srp = 0.0 if in_eclipse else \
            SOLAR_PRESSURE * self.srp_area * (1.0 + self.reflectivity) * self.cop_offset
        # Magnetic torque: residual dipole in the (tilted) dipole field
        sin_mag_lat = np.dot(satellite.pos_ecf / np.linalg.norm(satellite.pos_ecf),
                             self.dipole_axis_ecf)
        b_field = B0_EQUATOR * (R_EARTH / r) ** 3 * np.sqrt(1.0 + 3.0 * sin_mag_lat ** 2)
        t_mag = self.residual_dipole * b_field

        t_total = t_gg + t_aero + t_srp + t_mag  # Conservative: magnitudes add up
        self.momentum += t_total * sm.time_step
        self.metric[sm.cnt_epoch] = [t_gg, t_aero, t_srp, t_mag, t_total, self.momentum]

    def after_loop(self, sm):
        if not self.enabled:
            return
        times = np.asarray(self.times_f_doy)
        satellite = sm.satellites[self.idx_found_satellite]
        period = 2.0 * np.pi * np.sqrt((np.linalg.norm(satellite.pos_eci)) ** 3 / GM_EARTH)
        days = sm.num_epoch * sm.time_step / 86400.0
        per_orbit = self.momentum * period / (sm.num_epoch * sm.time_step)
        ls.logger.info(f'{self.type}: worst-case torques [N m] gravity gradient '
                       f'{np.nanmax(self.metric[:, 0]):.2e}, aero '
                       f'{np.nanmax(self.metric[:, 1]):.2e}, SRP '
                       f'{np.nanmax(self.metric[:, 2]):.2e}, magnetic '
                       f'{np.nanmax(self.metric[:, 3]):.2e}; momentum buildup '
                       f'{per_orbit:.3f} N m s/orbit, {self.momentum / days:.2f} N m s/day')

        fig, ax1 = plt.subplots(figsize=(10, 6))
        labels = ['Gravity gradient', 'Aerodynamic', 'Solar radiation pressure', 'Magnetic']
        colors = ['b', 'g', 'y', 'm']
        for i, (label, color) in enumerate(zip(labels, colors)):
            ax1.semilogy(times, self.metric[:, i], color + '-', linewidth=0.8, label=label)
        ax1.semilogy(times, self.metric[:, 4], 'k-', linewidth=1.2, label='Total')
        ax1.set_xlabel('DOY [-]')
        ax1.set_ylabel('Disturbance torque [N m]')
        ax1.grid(True)
        ax1.legend(loc='upper left', fontsize=8)
        ax2 = ax1.twinx()
        ax2.plot(times, self.metric[:, 5], 'r-', linewidth=1.2, label='Momentum buildup')
        if self.wheel_momentum is not None:
            ax2.axhline(self.wheel_momentum, color='red', linestyle=':',
                        linewidth=1.0, label='Wheel capacity')
        ax2.set_ylabel('Angular momentum [N m s]', color='red')
        ax2.tick_params(axis='y', colors='red')
        ax2.set_ylim(bottom=0)
        ax2.legend(loc='upper right', fontsize=8)
        fig.suptitle(f'AOCS disturbance torques, satellite {satellite.sat_id} '
                     f'({per_orbit:.3f} N m s/orbit momentum buildup)')
        fig.tight_layout()
        plt.savefig(sm.output_path(self.type + '.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'torque_gravity_gradient_nm', 'torque_aero_nm',
                            'torque_srp_nm', 'torque_magnetic_nm', 'torque_total_nm',
                            'momentum_nms'],
                       np.column_stack([times, self.metric]))


# ---------------------------------------------------------------------------
# Power subsystem analyses (historically the pow_* types in analysis_pow.py)
# ---------------------------------------------------------------------------

class AnalysisSatBatteryDepthDischarge(AnalysisBase):
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
        # 1. Determine Sun Visibility (shared cylindrical-shadow eclipse test)
        satellite = sm.satellites[0]
        sun_dir, in_eclipse = sun_direction_and_eclipse(satellite, sm.time_mjd)

        # 2. Power Generation
        solar_constant = 1361.0 # W/m^2
        p_gen = 0 if in_eclipse else (solar_constant * self.panel_area * self.efficiency)

        # 3. Power Consumption: instrument only active below the latitude limit
        is_active = abs(np.degrees(satellite.lla[0])) <= self.lat_limit
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
        plt.savefig(sm.output_path('sat_battery_depth_discharge.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'state_of_charge', 'power_generated_w', 'power_draw_w'],
                       self.metric)


class AnalysisSatEclipseDuration(AnalysisBase):
    def __init__(self):
        super().__init__()
        self.eclipse_durations = [] # List of (Time, Duration_minutes)
        self.in_eclipse_prev = False
        self.eclipse_start_time = 0

    def read_config(self, node):
        # No specific config needed for this, but could add altitude-specific masks if desired
        pass

    def before_loop(self, sm):
        self.eclipse_durations = []
        self.in_eclipse_prev = False
        ls.logger.info("Eclipse Duration Analysis Initialized")

    def in_loop(self, sm):
        satellite = sm.satellites[0]
        sun_dir, in_eclipse_now = sun_direction_and_eclipse(satellite, sm.time_mjd)

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
        self.write_csv(sm, ['doy', 'eclipse_duration_min'], self.eclipse_durations)
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

        plt.savefig(sm.output_path('sat_eclipse_duration.png'))
        plt.show()


# ---------------------------------------------------------------------------
# Data handling subsystem analyses (historically the dat_* types in analysis_dat.py)
# ---------------------------------------------------------------------------

class AnalysisSatDataStorage(AnalysisBase):
    def __init__(self):
        super().__init__()
        self.ssr_capacity_gbits = 0
        self.initial_fill_gbits = 0
        self.instrument_rate_mbps = 0
        self.downlink_rate_mbps = 0
        self.lat_limit = 90.0
        self.metric = None # Stores [Time, SSR_Level_Gbits, Is_Downlinking, Is_Recording]

    def read_config(self, node):
        self.ssr_capacity_gbits = float(node.find('SSRCapacityGbits').text)
        self.initial_fill_gbits = float(node.find('InitialFillGbits').text)
        self.instrument_rate_mbps = float(node.find('InstrumentRateMbps').text)
        self.downlink_rate_mbps = float(node.find('DownlinkRateMbps').text)
        if node.find('PayloadLatitudeLimit') is not None:
            self.lat_limit = float(node.find('PayloadLatitudeLimit').text)

    def before_loop(self, sm):
        self.metric = np.zeros((sm.num_epoch, 4))
        self.current_fill = self.initial_fill_gbits
        ls.logger.info("Data Storage Analysis Initialized")

    def in_loop(self, sm):
        sat = sm.satellites[0]
        sat.det_lla()

        # 1. Recording Logic (matching the Power module latitude trigger)
        is_recording = abs(np.degrees(sat.lla[0])) <= self.lat_limit

        # 2. Downlinking Logic (Ground Station in view)
        is_downlinking = len(sat.idx_stat_in_view) > 0

        # 3. Data Budget Calculation (Mbps * sec / 1000 = Gbits)
        inflow = (self.instrument_rate_mbps / 1000.0) * sm.time_step if is_recording else 0
        outflow = (self.downlink_rate_mbps / 1000.0) * sm.time_step if is_downlinking else 0

        self.current_fill += (inflow - outflow)

        # Constraints: Cannot be negative, cannot exceed capacity
        self.current_fill = np.clip(self.current_fill, 0, self.ssr_capacity_gbits)

        # Store [DOY, Fill Level, Downlink Status, Record Status]
        self.metric[sm.cnt_epoch, :] = [self.times_f_doy[-1], self.current_fill,
                                        float(is_downlinking), float(is_recording)]

    def after_loop(self, sm):
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(self.metric[:, 0], self.metric[:, 1], 'g-', label='SSR Fill Level')
        ax1.set_ylabel('Data Stored (Gbits)')
        ax1.set_xlabel('Day of Year (DOY)')

        # Shade regions to show activity for visual debugging
        ax1.fill_between(self.metric[:, 0], 0, self.ssr_capacity_gbits,
                        where=self.metric[:, 2] > 0, color='blue', alpha=0.1, label='Downlink Active')

        plt.grid(True)
        plt.legend()
        plt.savefig(sm.output_path('sat_data_storage.png'))
        plt.show()

        self.write_csv(sm, ['doy', 'ssr_fill_gbits', 'is_downlinking', 'is_recording'],
                       self.metric)


class AnalysisSatDataLatency(AnalysisBase):
    def __init__(self):
        super().__init__()
        self.ssr_capacity_gbits = 0
        self.initial_fill_gbits = 0
        self.instrument_rate_mbps = 0
        self.downlink_rate_mbps = 0
        self.lat_limit = 90.0
        self.ground_proc_min = 0 # x minutes for ground processing

        self.data_queue = []
        self.latency_metrics = []
        self.metric = None

    def read_config(self, node):
        self.ssr_capacity_gbits = float(node.find('SSRCapacityGbits').text)
        self.initial_fill_gbits = float(node.find('InitialFillGbits').text)
        self.instrument_rate_mbps = float(node.find('InstrumentRateMbps').text)
        self.downlink_rate_mbps = float(node.find('DownlinkRateMbps').text)
        if node.find('PayloadLatitudeLimit') is not None:
            self.lat_limit = float(node.find('PayloadLatitudeLimit').text)
        # Load the ground processing delay (x) from XML
        if node.find('GroundProcessingMin') is not None:
            self.ground_proc_min = float(node.find('GroundProcessingMin').text)

    def before_loop(self, sm):
        self.metric = np.zeros((sm.num_epoch, 4))
        self.current_fill = self.initial_fill_gbits
        if self.initial_fill_gbits > 0:
            self.data_queue.append([sm.time_mjd, self.initial_fill_gbits])
        ls.logger.info(f"Data Analysis Initialized with {self.ground_proc_min}min ground delay")

    def in_loop(self, sm):
        sat = sm.satellites[0]
        sat.det_lla()

        is_recording = abs(np.degrees(sat.lla[0])) <= self.lat_limit
        if is_recording:
            generated_gbits = (self.instrument_rate_mbps / 1000.0) * sm.time_step
            self.data_queue.append([sm.time_mjd, generated_gbits])
            self.current_fill += generated_gbits

        is_downlinking = len(sat.idx_stat_in_view) > 0
        if is_downlinking and self.current_fill > 0:
            downlink_capacity = (self.downlink_rate_mbps / 1000.0) * sm.time_step

            while downlink_capacity > 0 and len(self.data_queue) > 0:
                packet_time, packet_size = self.data_queue[0]

                if packet_size <= downlink_capacity:
                    # Latency = (Time on Orbit) + (Ground Processing x)
                    latency_h = ((sm.time_mjd - packet_time) * 24.0) + (self.ground_proc_min / 60.0)
                    self.latency_metrics.append([self.times_f_doy[-1], latency_h])
                    downlink_capacity -= packet_size
                    self.current_fill -= packet_size
                    self.data_queue.pop(0)
                else:
                    self.data_queue[0][1] -= downlink_capacity
                    self.current_fill -= downlink_capacity
                    latency_h = ((sm.time_mjd - packet_time) * 24.0) + (self.ground_proc_min / 60.0)
                    self.latency_metrics.append([self.times_f_doy[-1], latency_h])
                    downlink_capacity = 0

        self.current_fill = np.clip(self.current_fill, 0, self.ssr_capacity_gbits)
        self.metric[sm.cnt_epoch, :] = [self.times_f_doy[-1], self.current_fill,
                                        float(is_downlinking), float(is_recording)]

    def after_loop(self, sm):
        self.write_csv(sm, ['doy', 'latency_hours'], self.latency_metrics)
        if not self.latency_metrics:
            return

        lat_data = np.array(self.latency_metrics)
        latencies = lat_data[:, 1]

        # Calculate statistics
        mean_lat = np.mean(latencies)
        p95_lat = np.percentile(latencies, 95)
        max_lat = np.max(latencies)

        # Calculate percentage < 2 hours
        pct_under_2h = (np.sum(latencies < 2.0) / len(latencies)) * 100

        # Two separate plots (<type>_timeseries.png and <type>_histogram.png)
        # instead of a single two-panel figure
        fig, ax2 = plt.subplots(figsize=(12, 6))
        ax2.scatter(lat_data[:, 0], latencies, c='red', s=5, alpha=0.3)
        ax2.axhline(2.0, color='black', linestyle=':', label='2h Threshold')
        ax2.set_xlabel('Day of Year (DOY)')
        ax2.set_ylabel('Latency (Hours)')
        ax2.set_title(f'{pct_under_2h:.1f}% of data received in < 2 hours')
        ax2.legend()
        ax2.grid(True)
        plt.tight_layout()
        plt.savefig(sm.output_path('sat_data_latency_timeseries.png'))
        plt.show()

        fig, ax3 = plt.subplots(figsize=(12, 6))
        ax3.hist(latencies, bins=50, color='skyblue', edgecolor='black', alpha=0.7)
        ax3.axvline(mean_lat, color='blue', linestyle='--', label=f'Mean: {mean_lat:.2f}h')
        ax3.axvline(p95_lat, color='orange', linestyle='--', label=f'95%: {p95_lat:.2f}h')
        ax3.set_xlabel('Latency (Hours)')
        ax3.set_ylabel('Number of data packets [-]')
        ax3.set_title('Data latency histogram')
        ax3.legend()
        ax3.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(sm.output_path('sat_data_latency_histogram.png'))
        plt.show()
