"""
Regression tests: run every tests/<name>/Config.xml through the simulator
(into a scratch output directory) and compare the numeric CSV data dumps
against the golden copies stored in the test folder. Numbers are compared
with tolerances, so the tests are robust against matplotlib/cartopy version
changes that would break pixel comparison of the PNGs.

    py -m pytest tests                  # full suite (takes a while: MP4/HPOP tests)
    py -m pytest tests -m "not hpop"    # skip the Orekit-based tests
    py -m pytest tests -k sat_thermal   # one test
    py -m pytest tests --update-golden  # refresh the golden CSVs from this run

A test without golden CSVs is skipped with a hint (generate them with
--update-golden or tests/run_test.py, which also refreshes the plots).
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from conftest import TESTS_DIR, discover_tests

SRC_DIR = os.path.join(os.path.dirname(TESTS_DIR), 'src')
TIMEOUT_S = 2400
RTOL, ATOL = 1e-6, 1e-8


def run_simulator(config_file, output_dir):
    env = dict(os.environ, MPLBACKEND='Agg')
    result = subprocess.run(
        [sys.executable, 'main.py', config_file, '--output-dir', output_dir],
        cwd=SRC_DIR, env=env, capture_output=True, text=True, timeout=TIMEOUT_S)
    if result.returncode != 0:
        tail = '\n'.join(result.stderr.strip().splitlines()[-15:])
        pytest.fail(f'simulator exited with rc={result.returncode}:\n{tail}',
                    pytrace=False)


def compare_csv(golden_file, produced_file):
    golden = np.genfromtxt(golden_file, delimiter=',', names=True)
    produced = np.genfromtxt(produced_file, delimiter=',', names=True)
    assert golden.dtype.names == produced.dtype.names, \
        f'column names differ: {golden.dtype.names} vs {produced.dtype.names}'
    assert golden.shape == produced.shape, \
        f'row count differs: {golden.shape} vs {produced.shape}'
    for column in golden.dtype.names or []:
        assert np.allclose(golden[column], produced[column],
                           rtol=RTOL, atol=ATOL, equal_nan=True), \
            f'column {column} differs beyond tolerance'


@pytest.mark.parametrize('name', discover_tests())
def test_analysis(name, tmp_path, update_golden):
    test_dir = os.path.join(TESTS_DIR, name)
    output_dir = str(tmp_path / 'output')
    run_simulator(os.path.join(test_dir, 'Config.xml'), output_dir)

    if update_golden:
        copied = 0
        for produced in sorted(glob.glob(os.path.join(output_dir, '*.csv'))):
            shutil.copy2(produced, test_dir)
            copied += 1
        assert copied > 0, 'run produced no CSV data dumps to store as golden'
        return

    golden_files = sorted(glob.glob(os.path.join(test_dir, '*.csv')))
    if not golden_files:
        pytest.skip('no golden CSVs in the test folder yet - generate with '
                    '"py -m pytest tests --update-golden" or tests/run_test.py')
    for golden_file in golden_files:
        produced_file = os.path.join(output_dir, os.path.basename(golden_file))
        assert os.path.isfile(produced_file), \
            f'run did not produce {os.path.basename(golden_file)}'
        compare_csv(golden_file, produced_file)
