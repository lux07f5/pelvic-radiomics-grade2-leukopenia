#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run the whole analysis in order.

    python run_all.py                 # full pipeline
    python run_all.py --list          # show the steps and stop
    python run_all.py --skip-icc      # no inter-observer tables available
    python run_all.py --from 03       # resume from step 03
    python run_all.py --only 03       # run a single step
    python run_all.py --repeats 20    # shorter stability analysis

Every step writes its console output to outputs/logs/. A failing step stops the
run and its log tail is printed, so the error is visible without opening the
file.
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS_DIR = os.path.join(HERE, 'steps')

sys.path.insert(0, HERE)
import config  # noqa: E402  (resolves data/ and outputs/, honours RAD_DATA / RAD_OUT)

LOG_DIR = os.path.join(config.OUT_DIR, 'logs')

PY = sys.executable or 'python3'

# (number, script, title, needs icc inputs)
STEPS = [
    ('01', '01_build_dataset.py', 'build the analysis dataset', False),
    ('02', '02_table1_univariate.py', 'Table 1 and dosimetric univariable analysis', False),
    ('03', '03_radiomics_lasso.py', 'radiomics filter and LASSO, Rad-score', False),
    ('04', '04_model_comparison.py', 'model comparison and ROC', False),
    ('05', '05_calibration_dca.py', 'calibration and decision curves', False),
    ('06', '06_icc.py', 'inter-observer agreement', True),
    ('07', '07_icc_sensitivity.py', 'ICC-filtered sensitivity analysis', True),
    ('08', '08_mainline.py', 'pelvis-volume-first main analysis', False),
    ('09', '09_stability.py', 'repeated-CV LASSO stability', False),
    ('10', '10_figures_main.py', 'Figures 1 to 3', False),
    ('11', '11_figures_supp.py', 'Supplementary Figures S1 and S2', False),
]


def icc_ready():
    return all(os.path.exists(os.path.join(config.DATA_DIR, 'icc', f))
               for f in ('observer1.csv', 'observer2.csv', 'subset.csv'))


def run_step(num, script, title, repeats):
    path = os.path.join(STEPS_DIR, script)
    if not os.path.exists(path):
        return False, 'missing script: %s' % path

    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'
    if repeats:
        env['RAD_STABILITY_REPEATS'] = str(repeats)

    os.makedirs(LOG_DIR, exist_ok=True)
    log_path = os.path.join(LOG_DIR, '%s_%s.log' % (num, script[:-3]))
    t0 = time.time()
    print('\n' + '=' * 78)
    print('[%s] %s' % (num, title))
    print('=' * 78)
    with open(log_path, 'w', encoding='utf-8') as log:
        proc = subprocess.Popen([PY, path], cwd=HERE, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, encoding='utf-8', errors='replace')
        for line in proc.stdout:
            print(line, end='')
            log.write(line)
        proc.wait()

    elapsed = time.time() - t0
    if proc.returncode != 0:
        return False, 'exit code %d after %.1f s, see %s' % (proc.returncode, elapsed, log_path)
    print('-- done in %.1f s' % elapsed)
    return True, ''


def run_stability(repeats, jobs):
    """Step 09 in parallel: N workers over disjoint seed ranges, then combine.

    Single-process is fine but takes hours; the worker script writes one JSON
    line per completed repeat, so an interrupted parallel run resumes instead of
    starting over.
    """
    repeats = repeats or 100
    env = dict(os.environ)
    env['PYTHONIOENCODING'] = 'utf-8'
    env['RAD_STABILITY_REPEATS'] = str(repeats)

    os.makedirs(LOG_DIR, exist_ok=True)
    log_path = os.path.join(LOG_DIR, '09_09_stability.log')
    t0 = time.time()
    print('\n' + '=' * 78)
    print('[09] repeated-CV LASSO stability (%d repeats on %d workers)' % (repeats, jobs))
    print('=' * 78)

    with open(log_path, 'w', encoding='utf-8') as log:
        procs = []
        for i in range(jobs):
            a, b = repeats * i // jobs, repeats * (i + 1) // jobs
            outf = os.path.join(config.OUT_DIR, 'stability_w%d.jsonl' % (i + 1))
            if a == b:
                continue
            cmd = [PY, os.path.join(STEPS_DIR, '09_stability_worker.py'), str(a), str(b), outf]
            procs.append(subprocess.Popen(cmd, cwd=HERE, env=env,
                                          stdout=log, stderr=subprocess.STDOUT))
        for p in procs:
            p.wait()
        failed = [p.returncode for p in procs if p.returncode != 0]
        if failed:
            return False, 'worker exit codes %s, see %s' % (failed, log_path)

        cmd = [PY, os.path.join(STEPS_DIR, '09_stability_combine.py')]
        rc = subprocess.call(cmd, cwd=HERE, env=env, stdout=log, stderr=subprocess.STDOUT)
        if rc != 0:
            return False, 'combine failed with code %d, see %s' % (rc, log_path)

    with open(log_path, encoding='utf-8') as fh:
        tail = fh.read().strip().split('\n')[-12:]
    print('\n'.join(tail))
    print('-- done in %.1f s' % (time.time() - t0))
    return True, ''


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--list', action='store_true', help='list the steps and exit')
    ap.add_argument('--only', metavar='NN', help='run a single step, e.g. 03')
    ap.add_argument('--from', dest='start', metavar='NN', help='start at step NN')
    ap.add_argument('--skip-icc', action='store_true', help='skip the two ICC steps')
    ap.add_argument('--repeats', type=int, metavar='N',
                    help='repeats for the stability analysis (default 100)')
    ap.add_argument('--jobs', type=int, default=1, metavar='N',
                    help='parallel workers for the stability step (default 1)')
    args = ap.parse_args()

    if args.list:
        for num, script, title, needs in STEPS:
            mark = ' (needs data/icc/)' if needs else ''
            print('%s  %-28s %s%s' % (num, script, title, mark))
        return

    has_icc = icc_ready()
    if not has_icc and not args.skip_icc:
        print('note: data/icc/ is incomplete, steps 06 and 07 will be skipped')
        print('      (this is expected when only the single-observer table is available)')

    selected = []
    started = args.start is None
    for num, script, title, needs in STEPS:
        if args.only and num != args.only:
            continue
        if args.start and not started:
            if num != args.start:
                continue
            started = True
        if needs and not has_icc:
            continue
        selected.append((num, script, title))

    if not selected:
        sys.exit('no step matched the selection')

    print('pipeline: %d step(s)' % len(selected))
    print('data    : %s' % config.DATA_DIR)
    print('outputs : %s' % config.OUT_DIR)

    for num, script, title in selected:
        if num == '09' and args.jobs > 1:
            ok, msg = run_stability(args.repeats, args.jobs)
        else:
            ok, msg = run_step(num, script, title, args.repeats)
        if not ok:
            print('\nFAILED: %s' % msg)
            sys.exit(1)

    print('\nall steps completed')


if __name__ == '__main__':
    main()
