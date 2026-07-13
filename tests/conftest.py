"""
Pytest configuration of the SatSVS regression suite (see test_analyses.py).
"""
import os

import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))


def pytest_addoption(parser):
    parser.addoption(
        '--update-golden', action='store_true', default=False,
        help='Refresh the golden CSV data dumps in the test folders from this run '
             'instead of comparing against them')


@pytest.fixture
def update_golden(request):
    return request.config.getoption('--update-golden')


def discover_tests():
    """One test case per tests/<name>/Config.xml, marked hpop when the
    scenario uses the Orekit propagator (deselect with -m 'not hpop')."""
    cases = []
    for name in sorted(os.listdir(TESTS_DIR)):
        config = os.path.join(TESTS_DIR, name, 'Config.xml')
        if not os.path.isfile(config):
            continue
        marks = []
        with open(config) as f:
            if '<OrbitPropagator>HPOP</OrbitPropagator>' in f.read():
                marks.append(pytest.mark.hpop)
        cases.append(pytest.param(name, marks=marks))
    return cases
