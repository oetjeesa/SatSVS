"""
Test runner for SatSVS analysis tests.

Each test lives in tests/<name>/ and contains a Config.xml (plus any input
files it needs, e.g. TLE files referenced as ../tests/<name>/xxx.txt since
the simulator runs from src/).

Usage (from anywhere):
    py tests/run_test.py <test_name> [more_test_names...]

For each test the runner:
 1. copies tests/<name>/Config.xml to input/Config.xml
 2. runs `python main.py` from src/ with MPLBACKEND=Agg (headless plots)
 3. copies every file in output/ that was created/updated by the run
    (including main.log) back into tests/<name>/
"""
import glob
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
SRC = os.path.join(ROOT, 'src')
OUTPUT = os.path.join(ROOT, 'output')
TIMEOUT_S = 2400


def run_test(name):
    test_dir = os.path.join(ROOT, 'tests', name)
    cfg = os.path.join(test_dir, 'Config.xml')
    if not os.path.isfile(cfg):
        print(f'[{name}] FAIL: no Config.xml in {test_dir}')
        return False
    shutil.copyfile(cfg, os.path.join(ROOT, 'input', 'Config.xml'))

    # Snapshot output/ so only files created or rewritten by this run are collected
    before = {f: os.path.getmtime(f) for f in glob.glob(os.path.join(OUTPUT, '*'))
              if os.path.isfile(f)}

    env = dict(os.environ, MPLBACKEND='Agg')
    t_start = time.time()
    try:
        res = subprocess.run([sys.executable, 'main.py'], cwd=SRC, env=env,
                             capture_output=True, text=True, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        print(f'[{name}] FAIL: timeout after {TIMEOUT_S}s')
        return False
    elapsed = time.time() - t_start

    # Copy every output file touched by this run into the test folder
    copied = []
    for f in sorted(glob.glob(os.path.join(OUTPUT, '*'))):
        if os.path.isfile(f) and os.path.getmtime(f) > before.get(f, -1):
            shutil.copy2(f, test_dir)
            copied.append(os.path.basename(f))

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


if __name__ == '__main__':
    names = sys.argv[1:]
    if not names:
        print('usage: py tests/run_test.py <test_name> [...]')
        sys.exit(2)
    ok = all([run_test(n) for n in names])
    sys.exit(0 if ok else 1)
