# -*- coding: utf-8 -*-
"""One call that recomputes the central quantities of the article and compares them with the archived results.

  1. Gaucher input files: rebuilt from data/gaucher_source/gaucherdata.zip (code/gaucher_make_csv.py), byte for byte.
  2. Hit counts of Table S20, Fig. 5 and Fig. S6: code/make_hits_three_variants.py, byte for byte with
     results/raw_permutations/hits_three_variants.npy.
  3. Observed honest Q2 and AUC of Routes A1, A1u and A2 for the five strains (Table 1): code/run_nested_observed.py,
     exactly equal to results/revision/I_fallback_seeds.json.
  4. The three evaluations of the ChEMBL benchmark (hypothesis tests of Table 3, Additional files 12, 17 and 18):
     benchmark/auswertung.py, auswertung_k.py and auswertung_m.py, run on a copy of the archived data sets and
     permutation results, equal to the archived auswertung.json of ergebnisse/, ergebnisse_k/ and ergebnisse_m/.
  5. Optional, when Additional files 1 and 2 and the main text are given: every sheet of Additional file 2 against its
     source (revision2/tables/check_af2.py).
  6. Pre-specification of the second data set: size and SHA-256 of revision/lone/PRESPEC.md, equal to the values
     recorded in Additional file 5.

Steps 2 and 3 need Additional file 3 (README.md, "One input has to be placed by hand"). Nothing is written into the
archive: every script runs with its output folders in a temporary directory. Running time about three minutes.

Usage:  python check_all.py [--af1 AF1.docx --af2 AF2.xlsx --main MAIN.docx] [--keep DIR]
Exit code 0 when every check passes.
"""
import argparse, filecmp, json, os, shutil, subprocess, sys, tempfile, time

ROOT = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument('--af1'); ap.add_argument('--af2'); ap.add_argument('--main')
ap.add_argument('--keep', help='folder for the outputs of the checks (default: temporary, deleted at the end)')
A = ap.parse_args()
TMP = A.keep or tempfile.mkdtemp(prefix='check_all_'); os.makedirs(TMP, exist_ok=True)
ENV = dict(os.environ, QSAR_RECALC_DIR=os.path.join(TMP, 'recalc'), QSAR_REV_OUT=os.path.join(TMP, 'revision'),
           PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
os.makedirs(ENV['QSAR_RECALC_DIR'], exist_ok=True); os.makedirs(ENV['QSAR_REV_OUT'], exist_ok=True)
RESULTS = []

def run(args, cwd=None):
    p = subprocess.run([sys.executable] + args, cwd=cwd or TMP, env=ENV, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError('%s failed:\n%s' % (' '.join(args), (p.stdout + p.stderr)[-2000:]))
    return p.stdout

def check(name, fn):
    t0 = time.time()
    try:
        ok, detail = fn()
    except Exception as e:
        ok, detail = False, str(e)
    RESULTS.append((name, ok, detail)); print('%-4s %-62s %5.0f s  %s' % ('OK' if ok else 'FAIL', name, time.time() - t0, detail), flush=True)

def gaucher():
    out = run([os.path.join(ROOT, 'code', 'gaucher_make_csv.py')])
    n = out.count('byte-identical with data/')
    return n == 2, '%d of 2 files byte-identical' % n

def hits():
    run([os.path.join(ROOT, 'code', 'make_hits_three_variants.py')])
    same = filecmp.cmp(os.path.join(ENV['QSAR_RECALC_DIR'], 'hits_three_variants.npy'),
                       os.path.join(ROOT, 'results', 'raw_permutations', 'hits_three_variants.npy'), shallow=False)
    return same, 'hits_three_variants.npy %s' % ('byte-identical' if same else 'DIFFERS')

def observed():
    run([os.path.join(ROOT, 'code', 'run_nested_observed.py')])
    new = json.load(open(os.path.join(ENV['QSAR_RECALC_DIR'], 'nested_observed.json')))
    ref = json.load(open(os.path.join(ROOT, 'results', 'revision', 'I_fallback_seeds.json')))
    bad = ['%s %s %s' % (s, r, m) for s in ('HI', 'SA', 'SPneu', 'SPyo', 'PA') for r in ('A1', 'A1u', 'A2')
           for m in ('q2', 'auc') if new[s][r][m] != ref[s][r][m]]
    return not bad, '30 values exactly equal' if not bad else 'different: ' + ', '.join(bad)

def benchmark():
    B = os.path.join(TMP, 'benchmark'); os.makedirs(B, exist_ok=True)
    for f in ('auswertung.py', 'auswertung_k.py', 'auswertung_m.py'):
        shutil.copy2(os.path.join(ROOT, 'benchmark', f), B)
    for d in ('datasets', 'datasets_k', 'datasets_m', 'ergebnisse', 'ergebnisse_k', 'ergebnisse_m'):
        if not os.path.exists(os.path.join(B, d)):
            shutil.copytree(os.path.join(ROOT, 'benchmark', d), os.path.join(B, d))
    res = []
    for script, erg in (('auswertung.py', 'ergebnisse'), ('auswertung_k.py', 'ergebnisse_k'), ('auswertung_m.py', 'ergebnisse_m')):
        os.remove(os.path.join(B, erg, 'auswertung.json'))
        run([os.path.join(B, script)], cwd=B)
        same = json.load(open(os.path.join(B, erg, 'auswertung.json'))) == json.load(open(os.path.join(ROOT, 'benchmark', erg, 'auswertung.json')))
        res.append((erg, same))
    return all(s for _, s in res), ', '.join('%s %s' % (e, 'equal' if s else 'DIFFERS') for e, s in res)

def prespec():
    import hashlib
    b = open(os.path.join(ROOT, 'revision', 'lone', 'PRESPEC.md'), 'rb').read(); h = hashlib.sha256(b).hexdigest()
    ok = len(b) == 3266 and h == '8b563093642b911eb175088fbaa4a9048826290755391482e165707f8d9451d2'
    return ok, '%d bytes, SHA-256 %s...%s' % (len(b), h[:8], h[-4:])

def af2():
    run([os.path.join(ROOT, 'revision', 'make_base.py')])
    R2 = os.path.join(ROOT, 'revision2', 'results')
    s10 = os.path.join(R2, 'TableS10_alle_pWerte 06.11.2026 0626.json'); tb = os.path.join(R2, 'TableS5_TieBreak_0.0167 06.11.2026 0334.json')
    args = [os.path.join(ROOT, 'revision2', 'tables', 'check_af2.py'), A.af2, ROOT, A.af1, s10, tb,
            '--base', os.path.join(ENV['QSAR_RECALC_DIR'], 'base.pkl')] + (['--main', A.main] if A.main else [])
    out = run(args)
    return True, out.strip().split('\n')[-1][:150]

check('1  Gaucher input files from the original data file', gaucher)
check('2  hit-count null distributions (Table S20, Figs. 5 and S6)', hits)
check('3  observed honest Q2 and AUC, Routes A1, A1u, A2 (Table 1)', observed)
check('4  benchmark evaluations (Table 3, Additional files 12, 17, 18)', benchmark)
if A.af1 and A.af2:
    check('5  Additional file 2 against its sources', af2)
else:
    print('     5  Additional file 2: not checked (give --af1 and --af2, and --main for Main_Table_1 to 3)')
check('6  pre-specification of the second data set (PRESPEC.md)', prespec)
if not A.keep:
    shutil.rmtree(TMP, ignore_errors=True)
nbad = sum(1 for _, ok, _ in RESULTS if not ok)
print('\n%d of %d checks passed' % (len(RESULTS) - nbad, len(RESULTS)))
sys.exit(1 if nbad else 0)
