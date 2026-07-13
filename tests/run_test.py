"""
Test runner for SatSVS analysis tests: runs a test scenario and refreshes
every output file (plots, CSV data dumps, movies, main.log) stored in its
test folder. The pytest suite (test_analyses.py) compares the CSV dumps of a
fresh run against these files; use this runner to (re)generate them together
with the reference plots.

Each test lives in tests/<name>/ and contains a Config.xml (plus any input
files it needs, e.g. TLE files referenced as ../tests/<name>/xxx.txt since
the simulator runs from src/).

Usage (from anywhere):
    py tests/run_test.py <test_name> [more_test_names...]
    py tests/run_test.py --all

Each test runs with `python main.py <config> --output-dir <scratch>` into a
temporary directory (outside the repo, so sync clients cannot lock files),
after which every produced file is copied into the test folder.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
SRC = os.path.join(ROOT, 'src')
TIMEOUT_S = 2400


def run_test(name):
    test_dir = os.path.join(ROOT, 'tests', name)
    cfg = os.path.join(test_dir, 'Config.xml')
    if not os.path.isfile(cfg):
        print(f'[{name}] FAIL: no Config.xml in {test_dir}')
        return False

    out_dir = tempfile.mkdtemp(prefix=f'satsvs_{name}_')
    env = dict(os.environ, MPLBACKEND='Agg')
    t_start = time.time()
    try:
        res = subprocess.run([sys.executable, 'main.py', cfg, '--output-dir', out_dir],
                             cwd=SRC, env=env, capture_output=True, text=True,
                             timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        print(f'[{name}] FAIL: timeout after {TIMEOUT_S}s')
        return False
    elapsed = time.time() - t_start

    # Copy every output file of this run into the test folder
    copied = []
    for f in sorted(glob.glob(os.path.join(out_dir, '*'))):
        if os.path.isfile(f) and os.path.basename(f) != 'orbits_internal.txt':
            shutil.copy2(f, test_dir)
            copied.append(os.path.basename(f))
    shutil.rmtree(out_dir, ignore_errors=True)

    # Keep real errors, drop the benign headless-backend warning plt.show() emits under Agg
    stderr_clean = '\n'.join(
        ln for ln in res.stderr.splitlines()
        if 'FigureCanvasAgg is non-interactive' not in ln and ln.strip() != 'plt.show()')
    if stderr_clean.strip():
        with open(os.path.join(test_dir, 'run_stderr.txt'), 'w') as fh:
            fh.write(stderr_clean)

    status = 'OK' if res.returncode == 0 else f'FAIL rc={res.returncode}'
    print(f'[{name}] {status} in {elapsed:.1f}s, outputs: {", ".join(copied) or "none"}')
    if res.returncode != 0:
        tail = res.stderr.strip().splitlines()[-15:]
        print('--- stderr tail ---')
        print('\n'.join(tail))
    return res.returncode == 0


def all_test_names():
    return sorted(name for name in os.listdir(os.path.join(ROOT, 'tests'))
                  if os.path.isfile(os.path.join(ROOT, 'tests', name, 'Config.xml')))


if __name__ == '__main__':
    names = sys.argv[1:]
    if names == ['--all']:
        names = all_test_names()
    if not names:
        print('usage: py tests/run_test.py <test_name> [...] | --all')
        sys.exit(2)
    ok = all([run_test(n) for n in names])
    sys.exit(0 if ok else 1)
