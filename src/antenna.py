"""
Antenna gain patterns from the common antenna tool file formats, for the
com_* link-budget analyses (<TransmitAntennaPatternFile> /
<ReceiveAntennaPatternFile>).

Supported formats:
- TICRA GRASP tabulated cut files (.cut): polar cuts (ICUT=1), any number of
  phi cuts and field components.
- TICRA GRASP grid files (.grd, ASCII): uv far-field grids (IGRID=1),
  KLIMIT 0 or 1.

The far field is assumed stored with the standard GRASP power normalisation,
so the gain in dBi of a point is 10*log10(sum of re^2+im^2 over the field
components). The 2D pattern is reduced to gain versus off-boresight angle
theta by power-averaging over all phi cuts / grid directions — the
link-budget analyses treat antennas as rotationally symmetric around their
boresight. Example files in the real formats are provided in
input/example_antenna_patterns (see make_patterns.py there).
"""
import os

import numpy as np
import matplotlib.pyplot as plt

import logging_svs as ls
import misc_fn


class AntennaPattern:
    """Gain [dBi] versus off-boresight angle [deg], loaded from a GRASP
    .cut or .grd file. Use gain(off_boresight_rad) in the link budgets."""

    def __init__(self, theta_deg, gain_dbi, name):
        self.theta_deg = np.asarray(theta_deg, dtype=float)
        self.gain_dbi = np.asarray(gain_dbi, dtype=float)
        self.name = name

    @property
    def peak(self):
        """Peak gain [dBi] (the tracking-antenna boresight gain)."""
        return float(np.max(self.gain_dbi))

    def gain(self, off_boresight_rad):
        """Gain [dBi] at the off-boresight angle [rad]; clamped to the first/
        last tabulated value outside the tabulated range."""
        return float(np.interp(np.degrees(abs(off_boresight_rad)),
                               self.theta_deg, self.gain_dbi))

    @classmethod
    def from_file(cls, file_name):
        """Load a pattern; the format is chosen by the file extension."""
        file_name = misc_fn.resolve_path(file_name)
        ext = os.path.splitext(file_name)[1].lower()
        if ext == '.cut':
            theta, gain = _read_cut(file_name)
        elif ext == '.grd':
            theta, gain = _read_grd(file_name)
        else:
            ls.logger.error(f'Unknown antenna pattern format {ext} of {file_name} '
                            f'(use a GRASP .cut or .grd file)')
            exit()
        pattern = cls(theta, gain, os.path.basename(file_name))
        ls.logger.info(f'Loaded antenna pattern {pattern.name}: peak '
                       f'{pattern.peak:.1f} dBi, theta {pattern.theta_deg[0]:.1f}..'
                       f'{pattern.theta_deg[-1]:.1f} deg ({len(pattern.theta_deg)} points)')
        return pattern


def _fold_to_theta_curve(theta_deg, power):
    """(theta, power) samples of any sign/cut folded into a single gain [dBi]
    versus |theta| curve, power-averaging coincident angles."""
    theta_deg = np.abs(np.asarray(theta_deg, dtype=float))
    power = np.asarray(power, dtype=float)
    uniq, inverse = np.unique(np.round(theta_deg, 4), return_inverse=True)
    mean_power = np.bincount(inverse, weights=power) / np.bincount(inverse)
    return uniq, 10.0 * np.log10(np.maximum(mean_power, 1e-30))


def _read_cut(file_name):
    """GRASP tabulated cut file: per cut a text line, a header line
    'V_INI V_INC V_NUM C ICOMP ICUT NCOMP' and V_NUM field lines."""
    with open(file_name) as f:
        lines = [line for line in f.read().splitlines() if line.strip()]
    thetas, powers = [], []
    i = 0
    while i + 1 < len(lines):
        header = lines[i + 1].split()
        v_ini, v_inc, v_num = float(header[0]), float(header[1]), int(header[2])
        icut = int(header[5]) if len(header) > 5 else 1
        ncomp = int(header[6]) if len(header) > 6 else 2
        if icut != 1:
            ls.logger.warning(f'{file_name}: only polar cuts (ICUT=1) are '
                              f'supported, cut at line {i + 2} skipped')
            i += 2 + v_num
            continue
        for j in range(v_num):
            fields = [float(x) for x in lines[i + 2 + j].split()]
            thetas.append(v_ini + v_inc * j)
            powers.append(sum(f * f for f in fields[:2 * ncomp]))
        i += 2 + v_num
    if not thetas:
        ls.logger.error(f'No polar cut data found in {file_name}')
        exit()
    return _fold_to_theta_curve(thetas, powers)


def _read_grd(file_name):
    """GRASP ASCII grid file: header up to the ++++ line, then KTYPE,
    'NSET ICOMP NCOMP IGRID', IX/IY per set, 'XS YS XE YE', 'NX NY KLIMIT'
    and the field values (x index varying fastest). Only uv grids (IGRID=1)
    of the first field set are read."""
    with open(file_name) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines) and not lines[i].startswith('++++'):
        i += 1
    if i >= len(lines):
        ls.logger.error(f'No ++++ header terminator found in {file_name}')
        exit()
    tokens = ' '.join(lines[i + 1:]).split()
    pos = 0

    def take(n):
        nonlocal pos
        vals = [float(tokens[pos + k]) for k in range(n)]
        pos += n
        return vals

    take(1)  # KTYPE
    nset, icomp, ncomp, igrid = [int(v) for v in take(4)]
    if igrid != 1:
        ls.logger.error(f'{file_name}: only uv grids (IGRID=1) are supported, '
                        f'found IGRID={igrid}')
        exit()
    take(2 * nset)  # IX IY centre indices per set
    xs, ys, xe, ye = take(4)
    nx, ny, klimit = [int(v) for v in take(3)]
    thetas, powers = [], []
    for j in range(ny):
        if klimit == 1:
            i_s, i_n = [int(v) for v in take(2)]
        else:
            i_s, i_n = 1, nx
        v = ys + (ye - ys) * j / max(ny - 1, 1)
        for k in range(i_n):
            fields = take(2 * ncomp)
            u = xs + (xe - xs) * (i_s - 1 + k) / max(nx - 1, 1)
            sin_theta = np.hypot(u, v)
            if sin_theta > 1.0:  # Outside visible space
                continue
            thetas.append(np.degrees(np.arcsin(sin_theta)))
            powers.append(sum(f * f for f in fields))
    if not thetas:
        ls.logger.error(f'No visible-space grid data found in {file_name}')
        exit()
    # Bin the scattered grid directions into a theta curve
    thetas = np.asarray(thetas)
    bins = np.round(thetas / 0.1) * 0.1  # 0.1 deg bins
    return _fold_to_theta_curve(bins, powers)


def plot_patterns(patterns, file_name):
    """Verification plot of the loaded pattern(s): gain versus off-boresight
    angle. patterns: list of (label, AntennaPattern)."""
    fig = plt.figure(figsize=(10, 6))
    for label, pattern in patterns:
        plt.plot(pattern.theta_deg, pattern.gain_dbi, '-', linewidth=0.9,
                 label=f'{label}: {pattern.name} (peak {pattern.peak:.1f} dBi)')
    plt.xlabel('Off-boresight angle [deg]')
    plt.ylabel('Gain [dBi]')
    plt.grid(True)
    plt.legend(fontsize=8)
    plt.title('Antenna pattern(s) used in the link budget')
    plt.savefig(file_name)
    plt.close(fig)
    ls.logger.info(f'Saved antenna pattern plot to {file_name}')
