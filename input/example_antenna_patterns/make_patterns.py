"""
Generates the example antenna pattern files in this folder in the TICRA GRASP
file formats (.cut polar cuts and .grd uv grid), for the com_* link-budget
analyses (<TransmitAntennaPatternFile>/<ReceiveAntennaPatternFile>).

The patterns are realistic analytic models: a uniform circular aperture
(Airy pattern) with the ITU-R sidelobe envelope for the parabolic dishes, and
a shaped isoflux pattern typical for LEO X-band downlink antennas (the gain
increase towards the edge of coverage compensates the slant-range increase).
The far fields are stored with the standard GRASP power normalisation
(gain dBi = 10*log10(re^2 + im^2) summed over the components).

Run once:  py input/example_antenna_patterns/make_patterns.py
"""
import os

import numpy as np
from scipy.special import j1

HERE = os.path.dirname(os.path.abspath(__file__))
C_LIGHT = 299792458.0


def dish_gain(theta_deg, freq_hz, diameter, efficiency=0.6):
    """Uniform circular aperture (Airy) pattern with an ITU-R 32-25log(theta)
    style sidelobe envelope, in dBi."""
    lam = C_LIGHT / freq_hz
    theta = np.abs(np.asarray(theta_deg, dtype=float))
    x = np.pi * diameter / lam * np.sin(np.radians(theta))
    with np.errstate(divide='ignore', invalid='ignore'):
        airy = np.where(x > 1e-9, 2.0 * j1(x) / x, 1.0)
    peak = 10.0 * np.log10(efficiency * (np.pi * diameter / lam) ** 2)
    gain = peak + 20.0 * np.log10(np.maximum(np.abs(airy), 1e-12))
    envelope = np.maximum(32.0 - 25.0 * np.log10(np.maximum(theta, 0.1)), -10.0)
    return np.maximum(gain, np.minimum(envelope, peak - 20.0))


def isoflux_gain(theta_deg, nadir_dbi=3.0, peak_dbi=6.5, theta_peak=58.0,
                 theta_edge=64.0, floor_dbi=-15.0):
    """Shaped isoflux pattern of a LEO downlink antenna: gain rising from
    nadir to the edge-of-coverage peak, then rolling off steeply, in dBi."""
    theta = np.abs(np.asarray(theta_deg, dtype=float))
    rising = nadir_dbi + (peak_dbi - nadir_dbi) * (theta / theta_peak) ** 2
    falling = peak_dbi - 25.0 * ((theta - theta_peak) / (theta_edge - theta_peak)) ** 2
    gain = np.where(theta <= theta_peak, rising, falling)
    return np.maximum(gain, floor_dbi)


def write_cut(file_name, theta_deg, gain_dbi, text, phi_cuts=(0.0, 90.0)):
    """GRASP tabulated cut file: co-polar amplitude 10^(gain/20), tiny
    cross-polar component (ICOMP=3 linear co/cx, ICUT=1 polar cuts)."""
    amp = 10.0 ** (np.asarray(gain_dbi, dtype=float) / 20.0)
    v_ini, v_inc, v_num = theta_deg[0], theta_deg[1] - theta_deg[0], len(theta_deg)
    with open(os.path.join(HERE, file_name), 'w') as f:
        for phi in phi_cuts:
            f.write(f'{text}, phi = {phi:.0f} deg cut\n')
            f.write(f'{v_ini:.4f} {v_inc:.6f} {v_num} {phi:.1f} 3 1 2\n')
            for a in amp:
                f.write(f'{a:.6e} 0.000000e+00 1.000000e-08 0.000000e+00\n')
    print(f'wrote {file_name}: peak {np.max(gain_dbi):.1f} dBi')


def write_grd(file_name, freq_hz, diameter, text, u_max=0.2, n=201):
    """GRASP ASCII uv-grid file (KTYPE=1, IGRID=1, NCOMP=2, KLIMIT=0) of a
    circular-aperture dish pattern."""
    u = np.linspace(-u_max, u_max, n)
    uu, vv = np.meshgrid(u, u, indexing='xy')
    sin_theta = np.clip(np.hypot(uu, vv), 0.0, 1.0)
    gain = dish_gain(np.degrees(np.arcsin(sin_theta)), freq_hz, diameter)
    amp = 10.0 ** (gain / 20.0)
    with open(os.path.join(HERE, file_name), 'w') as f:
        f.write(f'{text}\n')
        f.write(f'Frequency [GHz]: {freq_hz / 1e9:.4f}\n')
        f.write('++++\n')
        f.write(' 1\n')  # KTYPE
        f.write(' 1 3 2 1\n')  # NSET ICOMP NCOMP IGRID
        f.write(' 0 0\n')  # IX IY
        f.write(f' {-u_max:.6f} {-u_max:.6f} {u_max:.6f} {u_max:.6f}\n')  # XS YS XE YE
        f.write(f' {n} {n} 0\n')  # NX NY KLIMIT
        for j in range(n):  # x index varying fastest
            for i in range(n):
                f.write(f'{amp[j, i]:.6e} 0.000000e+00 1.000000e-08 0.000000e+00\n')
    print(f'wrote {file_name}: peak {np.max(gain):.1f} dBi')


# LEO X-band downlink isoflux antenna (e.g. TerraSAR-X class downlink)
theta = np.arange(-90.0, 90.25, 0.25)
write_cut('isoflux_x_band.cut', theta, isoflux_gain(theta),
          'Shaped isoflux LEO X-band downlink antenna, 8.025 GHz')

# X-band 3 m ground station dish as uv grid (~46 dBi at 8.025 GHz)
write_grd('dish_x_band_3m.grd', 8.025e9, 3.0,
          'X-band 3 m ground station dish, 8.025 GHz, uniform aperture model')

# Ka-band 25 cm satellite downlink dish (~34 dBi at 26.25 GHz)
theta = np.arange(0.0, 90.05, 0.05)
write_cut('dish_ka_band_25cm.cut', theta, dish_gain(theta, 26.25e9, 0.25),
          'Ka-band 0.25 m satellite downlink dish, 26.25 GHz, aperture model')

# Ka-band 6.8 m ground station dish (~63 dBi at 26.25 GHz)
theta = np.arange(0.0, 90.01, 0.01)
write_cut('dish_ka_band_6_8m.cut', theta, dish_gain(theta, 26.25e9, 6.8),
          'Ka-band 6.8 m ground station dish, 26.25 GHz, aperture model')

print('done')
