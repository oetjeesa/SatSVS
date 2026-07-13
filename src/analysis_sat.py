"""
Satellite platform analyses: single-node thermal balance (sat_thermal) and
AOCS disturbance torques / momentum buildup (sat_aocs). Both work with any
orbit propagator; sat_aocs samples the HPOP atmosphere model when available
and falls back to a built-in exponential atmosphere otherwise.
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
    for the satellite (same shadow test as the pow_ analyses)."""
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
