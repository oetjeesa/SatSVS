"""
Benchmark of the HPOP (Orekit) propagator, two parts:

A) Correctness of the propagation plumbing: HPOP with ALL perturbations
   disabled must reproduce the analytical two-body propagation of the same
   Kepler elements (misc_fn.kep2xyz). Two-body motion is frame independent,
   so the GCRF coordinates from Orekit and the tool-frame coordinates from
   kep2xyz are directly comparable. Expected agreement: metres over one day
   (integrator tolerance); GM values are identical (398600.4415e9).

B) Realism against a known orbit: the full-force HPOP trajectory of
   TerraSAR-X (orbits.txt, produced by running this test's Config.xml through
   the simulator) is compared against the SGP4 trajectory of the same TLE
   (reference_orbit_sgp4.txt, the tool's established propagation path).
   Expected agreement: a few km over one day, dominated by the TLE mean-element
   fit, drag model differences and the TEME vs ITRF-based frame conventions.

Run AFTER the test itself (which produces orbits.txt in this folder):
    py tests/run_test.py hpop_benchmark
    py tests/hpop_benchmark/benchmark_hpop.py

Outputs benchmark_report.txt and hpop_benchmark_diff.png in this folder.
"""
import os
import sys
from math import radians

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(ROOT, 'src')
os.chdir(SRC)  # project modules use ../input and ../output relative paths
sys.path.insert(0, SRC)
os.environ.setdefault('MPLBACKEND', 'Agg')

import matplotlib.pyplot as plt
import misc_fn
from segments import Satellite
import propagation_hpop

# Thresholds
TWO_BODY_MAX_M = 10.0        # part A: max |diff| over one day [m]
VS_SGP4_MAX_KM = 50.0        # part B: max |diff| over one day [km]
VS_SGP4_RMS_KM = 20.0        # part B: rms |diff| over one day [km]

MJD0, MJD1, STEP = 61072.0, 61073.0, 60.0  # 2026-02-01, one day, 60 s


def part_a_two_body():
    """HPOP with every perturbation off vs analytical Keplerian propagation."""
    sat = Satellite()
    sat.sat_id = 1
    sat.kepler.epoch_mjd = MJD0
    sat.kepler.semi_major_axis = 7078137.0
    sat.kepler.eccentricity = 0.001
    sat.kepler.inclination = radians(97.4)
    sat.kepler.right_ascension = radians(50.0)
    sat.kepler.arg_perigee = radians(20.0)
    sat.kepler.mean_anomaly = radians(10.0)

    cfg = propagation_hpop.HpopConfig()
    cfg.geopotential = False
    cfg.drag = False
    cfg.solar_radiation_pressure = False
    cfg.third_body_sun = False
    cfg.third_body_moon = False
    cfg.third_body_planets = False
    cfg.solid_tides = False
    cfg.ocean_tides = False
    cfg.relativity = False
    cfg.integrator_position_tolerance = 0.001

    class Stub:
        pass
    sm = Stub()
    sm.hpop_config = cfg
    sm.satellites = [sat]
    sm.start_time = MJD0
    sm.stop_time = MJD1
    sm.time_step = STEP

    hpop = propagation_hpop.HpopPropagation(sm)

    mjds = np.arange(MJD0, MJD1 + 1e-9, STEP / 86400.0)
    diffs = np.zeros(len(mjds))
    for i, mjd in enumerate(mjds):
        pos_hpop, _ = hpop.sample_gcrf(0, float(mjd))
        pos_kep, _ = misc_fn.kep2xyz(float(mjd), sat.kepler.epoch_mjd,
                                     sat.kepler.semi_major_axis, sat.kepler.eccentricity,
                                     sat.kepler.inclination, sat.kepler.right_ascension,
                                     sat.kepler.arg_perigee, sat.kepler.mean_anomaly)
        diffs[i] = np.linalg.norm(pos_hpop - np.asarray(pos_kep))
    hours = (mjds - MJD0) * 24.0
    return hours, diffs


def part_b_vs_sgp4():
    """Full-force HPOP trajectory vs the SGP4 reference orbit file."""
    hpop_file = os.path.join(HERE, 'orbits.txt')
    ref_file = os.path.join(HERE, 'reference_orbit_sgp4.txt')
    if not os.path.isfile(hpop_file):
        raise SystemExit('orbits.txt missing - run "py tests/run_test.py hpop_benchmark" first')
    hpop = np.genfromtxt(hpop_file, delimiter=',')
    ref = np.genfromtxt(ref_file, delimiter=',')
    n = min(len(hpop), len(ref))
    hpop, ref = hpop[:n], ref[:n]
    if not np.allclose(hpop[:, 0], ref[:, 0]):
        raise SystemExit('epoch columns of orbits.txt and reference file do not match')
    diffs = np.linalg.norm(hpop[:, 2:5] - ref[:, 2:5], axis=1)
    hours = (hpop[:, 0] - hpop[0, 0]) * 24.0
    return hours, diffs


def main():
    hours_a, diff_a = part_a_two_body()
    hours_b, diff_b = part_b_vs_sgp4()

    a_max = float(np.max(diff_a))
    a_rms = float(np.sqrt(np.mean(diff_a ** 2)))
    b_max_km = float(np.max(diff_b)) / 1000.0
    b_rms_km = float(np.sqrt(np.mean(diff_b ** 2))) / 1000.0
    a_pass = a_max < TWO_BODY_MAX_M
    b_pass = b_max_km < VS_SGP4_MAX_KM and b_rms_km < VS_SGP4_RMS_KM

    lines = [
        'HPOP (Orekit) propagation benchmark - 2026-02-01, 1 day, 60 s sampling',
        '',
        'A) Two-body HPOP vs analytical Keplerian propagation (frame independent)',
        f'   max  position difference: {a_max:10.3f} m   (threshold {TWO_BODY_MAX_M} m)',
        f'   rms  position difference: {a_rms:10.3f} m',
        f'   -> {"PASS" if a_pass else "FAIL"}',
        '',
        'B) Full-force HPOP vs SGP4 reference trajectory (TerraSAR-X TLE, known orbit file)',
        f'   max  position difference: {b_max_km:10.3f} km  (threshold {VS_SGP4_MAX_KM} km)',
        f'   rms  position difference: {b_rms_km:10.3f} km  (threshold {VS_SGP4_RMS_KM} km)',
        f'   -> {"PASS" if b_pass else "FAIL"}',
        '',
        f'OVERALL: {"PASS" if (a_pass and b_pass) else "FAIL"}',
    ]
    report = '\n'.join(lines)
    print(report)
    with open(os.path.join(HERE, 'benchmark_report.txt'), 'w') as f:
        f.write(report + '\n')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    ax1.plot(hours_a, diff_a, 'b-')
    ax1.set_ylabel('|Δr| two-body vs Kepler [m]')
    ax1.set_title('A) HPOP two-body vs analytical Keplerian propagation')
    ax1.grid(True)
    ax2.plot(hours_b, diff_b / 1000.0, 'r-')
    ax2.set_ylabel('|Δr| HPOP vs SGP4 [km]')
    ax2.set_xlabel('Time since 2026-02-01 00:00 [h]')
    ax2.set_title('B) HPOP full force model vs SGP4 reference (TerraSAR-X)')
    ax2.grid(True)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'hpop_benchmark_diff.png'))

    sys.exit(0 if (a_pass and b_pass) else 1)


if __name__ == '__main__':
    main()
