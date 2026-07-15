"""
Mission report generator: writes a single self-contained report.html into a
results directory, assembling the run metadata, the per-analysis summary
lines from main.log and every produced plot, with links to the CSV data
dumps, pass tables and movies. Two ways to use it:

- In-run: <Report>True</Report> in the <SimulationManager> block writes the
  report at the end of the run (wired in main.py).
- Standalone on any existing results directory (e.g. a projects/ mission
  folder), no re-run needed:

      py report.py <results_dir> [--title "Mission name"]
"""
import argparse
import datetime
import glob
import html
import os
import re

from config_checks import ANALYSIS_PARAMS

# Files that are caches/inputs rather than results
_EXCLUDE = {'report.html', 'orbits_internal.txt', 'Config.xml'}

_CSS = """
body { font-family: 'Segoe UI', Arial, sans-serif; margin: 0; color: #222;
       background: #f5f6f8; }
.page { max-width: 1100px; margin: 0 auto; padding: 24px 32px 60px; }
h1 { margin: 8px 0 2px; font-size: 26px; }
.generated { color: #777; font-size: 13px; margin-bottom: 18px; }
table.meta { border-collapse: collapse; margin: 12px 0 24px; background: #fff;
             box-shadow: 0 1px 3px rgba(0,0,0,.12); }
table.meta td { border: 1px solid #ddd; padding: 6px 12px; font-size: 14px; }
table.meta td:first-child { font-weight: 600; background: #fafafa;
                            white-space: nowrap; }
nav { margin-bottom: 18px; font-size: 14px; }
nav a { margin-right: 14px; color: #0b62a4; text-decoration: none; }
nav a:hover { text-decoration: underline; }
section { background: #fff; border-radius: 6px; padding: 16px 20px 20px;
          margin-bottom: 22px; box-shadow: 0 1px 3px rgba(0,0,0,.12); }
h2 { margin: 2px 0 10px; font-size: 20px; border-bottom: 2px solid #0b62a4;
     padding-bottom: 4px; }
pre.log { background: #f4f7fa; border-left: 3px solid #0b62a4; padding: 8px 12px;
          font-size: 12.5px; overflow-x: auto; white-space: pre-wrap; }
pre.warn { background: #fdf3e7; border-left: 3px solid #d9822b; }
figure { margin: 14px 0; }
figure img { max-width: 100%; border: 1px solid #e0e0e0; }
figcaption { color: #777; font-size: 12px; margin-top: 2px; }
ul.files { font-size: 13.5px; }
footer { color: #999; font-size: 12px; margin-top: 30px; }
"""


def _known_types():
    """Registered analysis type names, longest first so e.g.
    obs_swath_push_broom matches before a shorter prefix would."""
    return sorted(ANALYSIS_PARAMS.keys(), key=len, reverse=True)


def _group_of(file_base, types):
    """Analysis type owning an output file name (without extension),
    including numbered repeated analyses (type_2)."""
    for type_name in types:
        if file_base == type_name or file_base.startswith(type_name + '_'):
            rest = file_base[len(type_name):]
            m = re.match(r'^_(\d+)(?:_|$)', rest)
            return type_name + (f'_{m.group(1)}' if m else '')
    return None


def _read_log(results_dir):
    """(all message, level) tuples from main.log, or [] when there is none."""
    log_file = os.path.join(results_dir, 'main.log')
    if not os.path.isfile(log_file):
        return []
    messages = []
    with open(log_file, errors='replace') as f:
        for line in f:
            m = re.search(r' - main - (\w+) - (.*)$', line.rstrip('\n'))
            if m:
                messages.append((m.group(2), m.group(1)))
    return messages


def write_report(results_dir, title=None, meta=None):
    """Write report.html into results_dir from the files found there.
    meta: optional list of (label, value) rows shown in the header table
    (the in-run caller passes precise scenario data; standalone use falls
    back to what the log file provides)."""
    types = _known_types()
    files = sorted(f for f in os.listdir(results_dir)
                   if os.path.isfile(os.path.join(results_dir, f))
                   and f not in _EXCLUDE and not f.endswith('.log'))
    groups = {}  # group -> {'png': [...], 'other': [...]}
    ungrouped = []
    for f in files:
        group = _group_of(os.path.splitext(f)[0], types)
        if group is None:
            ungrouped.append(f)
            continue
        kind = 'png' if f.lower().endswith('.png') else 'other'
        groups.setdefault(group, {'png': [], 'other': []})[kind].append(f)

    log = _read_log(results_dir)
    warnings = [m for m, level in log if level in ('WARNING', 'ERROR')]
    if title is None:
        title = os.path.basename(os.path.abspath(results_dir))
    if meta is None:
        meta = []
        for message, _ in log:
            if message.startswith('Loaded simulation'):
                meta.append(('Simulation', message.replace('Loaded simulation, ', '')))
                break

    e = html.escape
    out = [f'<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
           f'<title>{e(title)} - SatSVS report</title>\n<style>{_CSS}</style>\n'
           f'</head>\n<body>\n<div class="page">\n'
           f'<h1>{e(title)}</h1>\n'
           f'<div class="generated">SatSVS simulation report, generated '
           f'{datetime.datetime.now():%Y-%m-%d %H:%M}</div>\n']
    if meta:
        out.append('<table class="meta">\n')
        for label, value in meta:
            out.append(f'<tr><td>{e(str(label))}</td><td>{e(str(value))}</td></tr>\n')
        out.append('</table>\n')

    out.append('<nav>' + ' '.join(f'<a href="#{e(g)}">{e(g)}</a>'
                                  for g in sorted(groups)) + '</nav>\n')

    if warnings:
        out.append('<section><h2>Warnings</h2>\n<pre class="log warn">'
                   + e('\n'.join(warnings)) + '</pre></section>\n')

    for group in sorted(groups):
        out.append(f'<section id="{e(group)}">\n<h2>{e(group)}</h2>\n')
        summary = [m for m, level in log
                   if level == 'INFO' and m.startswith(group + ':')]
        if summary:
            out.append('<pre class="log">' + e('\n'.join(summary)) + '</pre>\n')
        for f in groups[group]['png']:
            out.append(f'<figure><img src="{e(f)}" alt="{e(f)}" loading="lazy">'
                       f'<figcaption>{e(f)}</figcaption></figure>\n')
        if groups[group]['other']:
            out.append('<ul class="files">\n')
            for f in groups[group]['other']:
                out.append(f'<li><a href="{e(f)}">{e(f)}</a></li>\n')
            out.append('</ul>\n')
        out.append('</section>\n')

    if ungrouped:
        out.append('<section id="other"><h2>Other files</h2>\n<ul class="files">\n')
        for f in ungrouped:
            out.append(f'<li><a href="{e(f)}">{e(f)}</a></li>\n')
        out.append('</ul></section>\n')

    out.append('<footer>Generated by SatSVS (report.py); the images and links '
               'are the files in this directory.</footer>\n</div>\n</body>\n</html>\n')
    report_file = os.path.join(results_dir, 'report.html')
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(''.join(out))
    return report_file


def write_report_from_sm(sm):
    """In-run report (<Report>True</Report>): title and metadata taken from
    the loaded scenario instead of being reconstructed from the log."""
    from astropy.time import Time
    config_file = os.path.abspath(sm.file_name)
    if os.path.basename(config_file).lower() == 'config.xml':
        title = os.path.basename(os.path.dirname(config_file))
    else:
        title = os.path.splitext(os.path.basename(config_file))[0]
    meta = [
        ('Scenario', config_file),
        ('Window', f'{Time(sm.start_time, format="mjd").iso[:19]} to '
                   f'{Time(sm.stop_time, format="mjd").iso[:19]} UTC, '
                   f'step {sm.time_step:g} s ({sm.num_epoch} epochs)'),
        ('Propagator', sm.orbit_propagator),
        ('Segments', f'{sm.num_sat} satellite(s), {len(sm.stations)} ground '
                     f'station(s), {len(sm.users)} user(s)'),
        ('Analyses', ', '.join(analysis.type for analysis in sm.analyses)),
    ]
    return write_report(sm.output_dir, title=title, meta=meta)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Build report.html from an existing SatSVS results directory')
    parser.add_argument('results_dir', help='directory with the run outputs '
                                            '(plots, CSVs, main.log)')
    parser.add_argument('--title', default=None, help='report title '
                        '(default: the directory name)')
    cli = parser.parse_args()
    print('written', write_report(cli.results_dir, title=cli.title))
