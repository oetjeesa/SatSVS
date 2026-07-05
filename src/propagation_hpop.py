"""
High Precision Orbit Propagation (HPOP) for SatSVS, based on the Orekit
astrodynamics library (java, accessed through the orekit_jpype bridge).

Selected with <OrbitPropagator>HPOP</OrbitPropagator> and configured through an
optional <HPOP> block in the SimulationManager section of Config.xml (see
readme.md for the full schema). Each perturbation is individually switchable:
geopotential (degree/order), atmospheric drag (area/Cd/atmosphere model), solar
radiation pressure, third-body attraction (Sun/Moon/planets), solid and ocean
tides, Earth pole rotation (full IERS EOP vs simplified transforms) and
relativity.

Requirements (only when HPOP is selected):
- python packages orekit_jpype and jdk4py (bundled JVM), and
- the Orekit physical data archive at input/orekit-data.zip
  (https://gitlab.orekit.org/orekit/orekit-data).

Frames: the propagation itself is done in GCRF with force models evaluated in
ITRF. The rest of SatSVS relates its ECI and ECF frames by a plain GMST
rotation, so to keep the ground geometry exact the ITRF position/velocity is
spun *forward* by GMST into the tool's pseudo-ECI frame before it is handed to
the Satellite object; Satellite.det_posvel_ecf then recovers the exact ITRF
coordinates. Satellite velocity therefore is the Earth-fixed (ITRF-relative)
velocity expressed in the tool frames, which is also the quantity the swath
and Doppler analyses actually need.

Each satellite is integrated once over the whole simulation window when the
propagator is built (dense-output ephemeris); the per-epoch calls in the time
loop only interpolate that ephemeris, which keeps the time loop fast.
"""
import os
from math import radians, degrees

import numpy as np

import misc_fn
import logging_svs as ls

OREKIT_DATA_ZIP = '../input/orekit-data.zip'

_vm_started = False


def _init_orekit():
    """Start the JVM (once) and load the Orekit physical data archive."""
    global _vm_started
    if _vm_started:
        return
    if 'JAVA_HOME' not in os.environ:
        try:
            import jdk4py
            os.environ['JAVA_HOME'] = str(jdk4py.JAVA_HOME)
        except ImportError:
            pass  # rely on a system JVM
    import orekit_jpype
    orekit_jpype.initVM()
    if not os.path.isfile(OREKIT_DATA_ZIP):
        ls.logger.error(f'HPOP needs the Orekit data archive at {OREKIT_DATA_ZIP} '
                        f'(download from https://gitlab.orekit.org/orekit/orekit-data)')
        exit()
    from orekit_jpype.pyhelpers import setup_orekit_data
    setup_orekit_data(filenames=OREKIT_DATA_ZIP, from_pip_library=False)
    _vm_started = True
    ls.logger.info('HPOP: Orekit JVM started and physical data loaded')


class HpopConfig:
    """Parsed <HPOP> configuration block (defaults give a full force model)."""

    def __init__(self):
        # Integrator (variable step Dormand-Prince 8(5,3))
        self.integrator_min_step = 0.001  # s
        self.integrator_max_step = 300.0  # s
        self.integrator_position_tolerance = 1.0  # m
        self.mass = 1000.0  # kg

        # Geopotential
        self.geopotential = True
        self.geopotential_degree = 21
        self.geopotential_order = 21

        # Earth orientation: True = full IERS EOP incl. pole rotation,
        # False = simplified (no polar motion / tidal EOP corrections)
        self.earth_pole_rotation = True

        # Atmospheric drag
        self.drag = True
        self.drag_area = 1.0  # m^2
        self.drag_cd = 2.2
        self.drag_model = 'NRLMSISE00'  # NRLMSISE00 | DTM2000 | HarrisPriester

        # Solar radiation pressure
        self.solar_radiation_pressure = True
        self.srp_area = 1.0  # m^2
        self.srp_cr = 1.5

        # Third bodies
        self.third_body_sun = True
        self.third_body_moon = True
        self.third_body_planets = False  # Venus, Mars, Jupiter point masses

        # Tides
        self.solid_tides = True
        self.ocean_tides = False
        self.ocean_tides_degree = 4
        self.ocean_tides_order = 4

        # Relativistic correction (Schwarzschild)
        self.relativity = False

    def read_config(self, node):
        def get_float(tag, default):
            return float(node.find(tag).text) if node.find(tag) is not None else default

        def get_int(tag, default):
            return int(node.find(tag).text) if node.find(tag) is not None else default

        def get_bool(tag, default):
            return misc_fn.str2bool(node.find(tag).text) if node.find(tag) is not None else default

        self.integrator_min_step = get_float('IntegratorMinStep', self.integrator_min_step)
        self.integrator_max_step = get_float('IntegratorMaxStep', self.integrator_max_step)
        self.integrator_position_tolerance = get_float('IntegratorPositionTolerance',
                                                       self.integrator_position_tolerance)
        self.mass = get_float('Mass', self.mass)

        self.geopotential = get_bool('Geopotential', self.geopotential)
        self.geopotential_degree = get_int('GeopotentialDegree', self.geopotential_degree)
        self.geopotential_order = get_int('GeopotentialOrder', self.geopotential_order)

        self.earth_pole_rotation = get_bool('EarthPoleRotation', self.earth_pole_rotation)

        self.drag = get_bool('Drag', self.drag)
        self.drag_area = get_float('DragArea', self.drag_area)
        self.drag_cd = get_float('DragCd', self.drag_cd)
        if node.find('DragModel') is not None:
            self.drag_model = node.find('DragModel').text

        self.solar_radiation_pressure = get_bool('SolarRadiationPressure',
                                                 self.solar_radiation_pressure)
        self.srp_area = get_float('SRPArea', self.srp_area)
        self.srp_cr = get_float('SRPCr', self.srp_cr)

        self.third_body_sun = get_bool('ThirdBodySun', self.third_body_sun)
        self.third_body_moon = get_bool('ThirdBodyMoon', self.third_body_moon)
        self.third_body_planets = get_bool('ThirdBodyPlanets', self.third_body_planets)

        self.solid_tides = get_bool('SolidTides', self.solid_tides)
        self.ocean_tides = get_bool('OceanTides', self.ocean_tides)
        self.ocean_tides_degree = get_int('OceanTidesDegree', self.ocean_tides_degree)
        self.ocean_tides_order = get_int('OceanTidesOrder', self.ocean_tides_order)

        self.relativity = get_bool('Relativity', self.relativity)

        ls.logger.info(f'HPOP configuration: {self.__dict__}')


class HpopPropagation:
    """Builds one Orekit numerical propagator per satellite and serves the
    per-epoch position/velocity requests of the SatSVS time loop."""

    def __init__(self, sm):
        _init_orekit()
        cfg = sm.hpop_config if sm.hpop_config is not None else HpopConfig()
        self.cfg = cfg

        from org.orekit.frames import FramesFactory
        from org.orekit.utils import IERSConventions, Constants
        from org.orekit.time import TimeScalesFactory
        from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
        from org.orekit.forces.gravity.potential import GravityFieldFactory

        self.utc = TimeScalesFactory.getUTC()
        self.gcrf = FramesFactory.getGCRF()
        # simpleEOP=True skips the tidal corrections of the EOP interpolation;
        # EarthPoleRotation False additionally drops polar motion by using the
        # equinox-based TOD-like transforms without EOP corrections
        if cfg.earth_pole_rotation:
            self.itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)
        else:
            self.itrf = FramesFactory.getITRFEquinox(IERSConventions.IERS_1996, True)
        self.mu = Constants.EIGEN5C_EARTH_MU
        self.earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                      Constants.WGS84_EARTH_FLATTENING, self.itrf)
        self.sun = CelestialBodyFactory.getSun()
        self.moon = CelestialBodyFactory.getMoon()
        self.gravity_provider = GravityFieldFactory.getNormalizedProvider(
            cfg.geopotential_degree, cfg.geopotential_order)

        # Dense-output ephemeris per satellite over the whole simulation window
        pad = 2 * sm.time_step
        self.date_start = self._mjd2date(sm.start_time)
        date_end = self._mjd2date(sm.stop_time + pad / 86400.0)
        self.ephemerides = []
        for satellite in sm.satellites:
            self.ephemerides.append(self._build_ephemeris(satellite, date_end))
        ls.logger.info(f'HPOP: generated dense ephemerides for {len(self.ephemerides)} satellite(s) '
                       f'from MJD {sm.start_time} to {sm.stop_time}')

    # ------------------------------------------------------------------ time
    def _mjd2date(self, mjd):
        from org.orekit.time import AbsoluteDate
        return AbsoluteDate.createMJDDate(int(mjd), (mjd - int(mjd)) * 86400.0, self.utc)

    # -------------------------------------------------------- initial states
    def _initial_orbit(self, satellite):
        """Initial osculating orbit at the simulation start: from the TLE if the
        satellite was defined by one, else from its configured Kepler elements."""
        from org.orekit.orbits import KeplerianOrbit, CartesianOrbit, PositionAngleType

        if satellite.tle_line1:  # TLE definition: SGP4 state at simulation start
            from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
            tle = TLE(satellite.tle_line1.strip(), satellite.tle_line2.strip())
            state = TLEPropagator.selectExtrapolator(tle).propagate(self.date_start)
            pv_gcrf = state.getPVCoordinates(self.gcrf)
            return CartesianOrbit(pv_gcrf, self.gcrf, self.date_start, self.mu)

        epoch = self._mjd2date(satellite.kepler.epoch_mjd)
        return KeplerianOrbit(satellite.kepler.semi_major_axis,
                              satellite.kepler.eccentricity,
                              satellite.kepler.inclination,
                              satellite.kepler.arg_perigee,
                              satellite.kepler.right_ascension,
                              satellite.kepler.mean_anomaly,
                              PositionAngleType.MEAN, self.gcrf, epoch, self.mu)

    def _build_propagator(self, satellite):
        from org.orekit.orbits import OrbitType
        from org.orekit.propagation import SpacecraftState, ToleranceProvider
        from org.orekit.propagation.numerical import NumericalPropagator
        from org.hipparchus.ode.nonstiff import DormandPrince853Integrator

        cfg = self.cfg
        orbit0 = self._initial_orbit(satellite)
        if not satellite.tle_line1:
            gap = abs(orbit0.getDate().durationFrom(self.date_start)) / 86400.0
            if gap > 30:
                ls.logger.warning(f'HPOP satellite {satellite.sat_id}: Kepler epoch is '
                                  f'{gap:.0f} days away from the simulation start; the whole '
                                  f'gap is integrated numerically, which is slow and drifts. '
                                  f'Set EpochMJD close to StartDate for HPOP runs.')

        tol = ToleranceProvider.getDefaultToleranceProvider(
            cfg.integrator_position_tolerance).getTolerances(orbit0, OrbitType.CARTESIAN)
        integrator = DormandPrince853Integrator(cfg.integrator_min_step,
                                                cfg.integrator_max_step, tol[0], tol[1])
        propagator = NumericalPropagator(integrator)
        propagator.setOrbitType(OrbitType.CARTESIAN)
        propagator.setInitialState(SpacecraftState(orbit0).withMass(cfg.mass))
        for force in self._force_models():
            propagator.addForceModel(force)
        return propagator

    def _force_models(self):
        """The perturbation set selected in the <HPOP> block. Note the Newtonian
        central attraction is always added by NumericalPropagator itself."""
        from org.orekit.forces.gravity import (HolmesFeatherstoneAttractionModel,
                                               ThirdBodyAttraction, SolidTides, OceanTides,
                                               Relativity)
        from org.orekit.forces.drag import DragForce, IsotropicDrag
        from org.orekit.forces.radiation import (SolarRadiationPressure,
                                                 IsotropicRadiationSingleCoefficient)
        from org.orekit.models.earth.atmosphere import NRLMSISE00, DTM2000, HarrisPriester
        from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
        from org.orekit.bodies import CelestialBodyFactory
        from org.orekit.utils import IERSConventions
        from org.orekit.time import TimeScalesFactory

        cfg = self.cfg
        forces = []

        if cfg.geopotential:
            forces.append(HolmesFeatherstoneAttractionModel(self.itrf, self.gravity_provider))

        if cfg.drag:
            model = cfg.drag_model.lower()
            if model == 'harrispriester':
                atmosphere = HarrisPriester(self.sun, self.earth)
            else:
                space_weather = CssiSpaceWeatherData(CssiSpaceWeatherData.DEFAULT_SUPPORTED_NAMES)
                if model == 'dtm2000':
                    atmosphere = DTM2000(space_weather, self.sun, self.earth)
                elif model == 'nrlmsise00':
                    atmosphere = NRLMSISE00(space_weather, self.sun, self.earth)
                else:
                    ls.logger.error(f'HPOP: unknown DragModel {cfg.drag_model} '
                                    f'(use NRLMSISE00, DTM2000 or HarrisPriester)')
                    exit()
            forces.append(DragForce(atmosphere, IsotropicDrag(cfg.drag_area, cfg.drag_cd)))

        if cfg.solar_radiation_pressure:
            forces.append(SolarRadiationPressure(
                self.sun, self.earth,
                IsotropicRadiationSingleCoefficient(cfg.srp_area, cfg.srp_cr)))

        if cfg.third_body_sun:
            forces.append(ThirdBodyAttraction(self.sun))
        if cfg.third_body_moon:
            forces.append(ThirdBodyAttraction(self.moon))
        if cfg.third_body_planets:
            for body in ('VENUS', 'MARS', 'JUPITER'):
                forces.append(ThirdBodyAttraction(getattr(CelestialBodyFactory, 'get' + body.capitalize())()))

        if cfg.solid_tides or cfg.ocean_tides:
            conventions = IERSConventions.IERS_2010
            ut1 = TimeScalesFactory.getUT1(conventions, not cfg.earth_pole_rotation)
            gp = self.gravity_provider
            if cfg.solid_tides:
                forces.append(SolidTides(self.itrf, gp.getAe(), gp.getMu(), gp.getTideSystem(),
                                         conventions, ut1, [self.sun, self.moon]))
            if cfg.ocean_tides:
                forces.append(OceanTides(self.itrf, gp.getAe(), gp.getMu(),
                                         cfg.ocean_tides_degree, cfg.ocean_tides_order,
                                         conventions, ut1))

        if cfg.relativity:
            forces.append(Relativity(self.mu))

        return forces

    def _build_ephemeris(self, satellite, date_end):
        propagator = self._build_propagator(satellite)
        generator = propagator.getEphemerisGenerator()
        propagator.propagate(self.date_start, date_end)
        return generator.getGeneratedEphemeris()

    # ------------------------------------------------------------ time loop
    def update_satellite(self, satellite, idx_sat, mjd, gmst):
        """Set satellite.pos_eci/vel_eci for the requested epoch. The values are
        the exact ITRF coordinates spun forward by GMST (tool pseudo-ECI), so
        that Satellite.det_posvel_ecf recovers the exact Earth-fixed state."""
        date = self._mjd2date(mjd)
        pv = self.ephemerides[idx_sat].getPVCoordinates(date, self.itrf)
        pos = pv.getPosition()
        vel = pv.getVelocity()
        pos_itrf = np.array([pos.getX(), pos.getY(), pos.getZ()])
        vel_itrf = np.array([vel.getX(), vel.getY(), vel.getZ()])
        satellite.pos_eci, satellite.vel_eci = misc_fn.spin_vector(gmst, pos_itrf, vel_itrf)

    def sample_gcrf(self, idx_sat, mjd):
        """Inertial GCRF position/velocity [m, m/s] at the requested epoch —
        used by the benchmark scripts to validate the propagation itself
        without the tool frame conventions."""
        date = self._mjd2date(mjd)
        pv = self.ephemerides[idx_sat].getPVCoordinates(date, self.gcrf)
        pos = pv.getPosition()
        vel = pv.getVelocity()
        return (np.array([pos.getX(), pos.getY(), pos.getZ()]),
                np.array([vel.getX(), vel.getY(), vel.getZ()]))
