"""Results workbook of the SENSITIVITY ANALYSIS MEDIAN RULE (addendum frozen 25.10.2026 21:39,
SHA-256 eaf21b7dff370396663d4fe93943a2fc4bd02044e9a45f00ffef94ef1b48eeab). Layout as ergebnismappe.py (main analysis),
purely descriptive: no decisions on H1 to H4, no Holm correction.
Sheets: README, Class_definition, Summary, Comparison_per_series, All_results_long, DS01-DS16, Robustness,
Predictions, Exploratory, Deviations.
Usage: python3 ergebnismappe_m.py <output file>
Nothing is recomputed except descriptive key figures (scaffold share, mean Tanimoto similarity) from the
SMILES and the exploratory summaries.
"""
import json, sys
from collections import Counter
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog('rdApp.*')

OUT = sys.argv[1]
ERG, DSD = 'ergebnisse_m', 'datasets_m'
A = json.load(open(f'{ERG}/auswertung.json')); T = A['tabelle']; H = A['hypothesen']
R = json.load(open(f'{ERG}/robustheit.json'))
V = json.load(open(f'{ERG}/vorhersagen.json'))
Z = json.load(open(f'{ERG}/maxz_labelfrei.json'))
AK = json.load(open('ergebnisse_k/auswertung.json')); TK = AK['tabelle']; HK = AK['hypothesen']
RK = json.load(open('ergebnisse_k/robustheit.json'))
PROT = json.load(open(f'{DSD}/protokoll.json'))
CAND = json.load(open('cand_ranked.json')); BYASSAY = {x['assay']: x for x in CAND}
PRIM = [f'DS{i:02d}' for i in range(1, 17)]; ALL = PRIM
INFO = {d: json.load(open(f'{DSD}/{d}_info.json'))['info'] for d in ALL}
INFOK = {d: json.load(open(f'datasets_k/{d}_info.json'))['info'] for d in ALL}
ROUTES = ['A1', 'A1u', 'A2']
MODELS = ['app_A1', 'app_A1u', 'app_A2', 'A1', 'A1u', 'A2', 'ohneAuswahl', 'ohneFilter', 'RF']
MNAME = {'app_A1': 'Route A1, apparent', 'app_A1u': 'Route A1u, apparent', 'app_A2': 'Route A2, apparent',
         'A1': 'Route A1, nested', 'A1u': 'Route A1u, nested', 'A2': 'Route A2, nested',
         'ohneAuswahl': 'No selection (PLS-DA, all cleaned descriptors)', 'ohneFilter': 'No filter (PLS-DA, label-free group representatives)',
         'RF': 'Random forest (Morgan fingerprints, OOB)'}
KNAME = {'ohneFilter': 'no filter', 'ohneAuswahl': 'no selection', 'RF': 'random forest'}
ART = {'MIC': 'MIC', 'Enzym': 'Enzyme IC50', 'Rezeptor': 'Receptor Ki', 'Zytotox': 'Cytotoxicity'}
DIV = {'kongener': 'congeneric', 'divers': 'diverse'}
BAND = {'A 25-35 %': 'A (25-35 %)', 'B >35-50 %': 'B (>35-50 %)'}
UNIT = {'ug.mL-1': 'ug/mL', 'nM': 'nM'}
F = 'Arial'; HF = Font(name=F, bold=True); NF = Font(name=F); TF = Font(name=F, bold=True, size=13)
FILL = PatternFill('solid', fgColor='D9E1F2'); ORA = PatternFill('solid', fgColor='F8CBAD'); YEL = PatternFill('solid', fgColor='FFF2CC')
wb = Workbook()


def table(ws, r0, cols, rows, widths=None):
    for j, c in enumerate(cols, 1):
        x = ws.cell(r0, j, c); x.font = HF; x.fill = FILL; x.alignment = Alignment(wrap_text=True, vertical='top')
    for i, row in enumerate(rows, r0 + 1):
        for j, v in enumerate(row, 1):
            if isinstance(v, (np.floating, np.integer)): v = v.item()
            if isinstance(v, float) and not np.isfinite(v): v = str(v)
            ws.cell(i, j, v).font = NF
    if widths:
        for j, w in enumerate(widths, 1):
            if w: ws.column_dimensions[get_column_letter(j)].width = max(ws.column_dimensions[get_column_letter(j)].width or 0, w)
    return r0 + len(rows) + 2


def sheet(name, cols, rows, widths=None):
    ws = wb.create_sheet(name); table(ws, 1, cols, rows, widths); ws.freeze_panes = 'A2'
    ws.auto_filter.ref = ws.dimensions; return ws


def text(ws, r, t, f=NF):
    ws.cell(r, 1, t).font = f; return r + 1


def zelle(d):
    z = INFO[d]['zelle']
    return (ART[z[0]], DIV[z[1]], BAND[z[2]])


def g4(x):
    return float(f'{x:.6g}')


GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
STRUK = {}
for d in ALL:
    tab = pd.read_csv(f'{DSD}/{d}_y.csv'); ms = [Chem.MolFromSmiles(s) for s in tab.smiles]
    sc = [Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m))) for m in ms]
    fps = [GEN.GetFingerprint(m) for m in ms]
    sims = [s for i in range(len(fps)) for s in DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])]
    STRUK[d] = dict(top=Counter(sc).most_common(1)[0][1] / len(sc), mts=float(np.mean(sims)))

# ---------------------------------------------------------------- README
ws = wb.active; ws.title = 'README'
L = [('ChEMBL benchmark – results workbook – SENSITIVITY ANALYSIS: MEDIAN RULE (descriptive only)', TF), ('', NF),
     ('Manuscript: Werle J, Rondevaldova J, Werle E. Honest performance estimation for small-sample QSAR (Journal of Cheminformatics, revision 2).', NF),
     ('Pre-specification frozen 22.10.2026 18:10, SHA-256 a13bb2c0b513289b160de744bd8d6bc47662d46553ab5a2fe3db6d336a0ae419. It names the median rule (active = better than the median of the data set) as sensitivity analysis, computed for all 16 data sets and reported descriptively only, seeds analysis code + 100.', NF),
     ('Supplementary specification (main analysis) frozen 24.10.2026 17:20, SHA-256 36668775daafd821c1871a15c6a697ea5e3c25afebbb3d6b6a4db818c902df87. The median rule is applied to its 16 series.', NF),
     ('Addendum for the median rule frozen 25.10.2026 21:39, before any model with median-rule classes was fitted, SHA-256 eaf21b7dff370396663d4fe93943a2fc4bd02044e9a45f00ffef94ef1b48eeab. It fixes what the pre-specification left open: treatment of values equal to the median and of censored values, and the seeds (analysis code + 300, because + 100 and + 200 are already used).', NF),
     ('Nothing in this workbook enters the decisions on H1 to H4. Tests are reported with their p-values, without decision and without Holm correction.', HF), ('', NF),
     ('Rule', HF),
     ('Same compounds and descriptor matrices as in the main analysis. Only the classes change. Boundary: median of the recorded values of the series (per compound the median of the log values, for censored-only compounds the bound). Active = value strictly below the median (lower = more potent). Equal to the median (log comparison, tolerance 1e-9) = inactive. ">" with bound at or above the median = inactive, ">" with bound below the median = not decidable, compound excluded. Exclusion occurs only in DS13 (6 compounds), which was rebuilt with the pipeline of the main analysis (cleaning, correlation groups, label-free representatives).', NF),
     ('DS02 and DS06 have the same classes as in the main analysis. Their permutation results were taken over (byte-identical), as specified. The bootstrap uses seed 308 for all 16 series.', NF), ('', NF),
     ('Models, permutations, inflation factors and bootstrap exactly as in the main analysis: nested models 2,000 permutations, apparent models 1,000, hit counts 1,000, random forest 500, permutation p = (1 + #{null >= observed - 1e-9}) / (1 + N), paired stratified bootstrap with 20,000 draws, 95 % percentile interval. Seeds numpy default_rng([93, data set number, analysis code + 300, permutation index]), bootstrap default_rng([93, data set number, 308]).', NF), ('', NF),
     ('Sheets', HF),
     ('Class_definition – median, classes and changes per series against the main analysis.', NF),
     ('Summary – the tests of H1 to H4 under the median rule next to the main analysis, and one row per data set.', NF),
     ('Comparison_per_series – main analysis against median rule per series: optimism, calibrated inflation factors, nested Q2 and AUC, H4 classification.', NF),
     ('All_results_long – every number of the data set sheets in long format.', NF),
     ('DS01 to DS16 – one sheet per data set, layout as in the main analysis.', NF),
     ('Robustness – R1, R2 and AUC instead of Q2, as in the main analysis, descriptive.', NF),
     ('Predictions – leave-one-out / OOB prediction of every compound by every model, with both classes and the label-free |z|.', NF),
     ('Exploratory – the descriptive groupings of the main analysis (band, diversity, activity type, A1u hit excess).', NF),
     ('Deviations – deviations from the addendum.', NF), ('', NF),
     ('Bands and diversity in the sheets refer to the cells of the design (main analysis). Under the median rule the minor class is 38 to 50 % in every series.', NF),
     ('SMILES, ChEMBL IDs, values, both class columns (y = median rule, y_ersatzregel = main analysis), descriptor matrices, raw permutation outputs and all scripts are in the code archive.', NF)]
for i, (t, f) in enumerate(L, 1):
    c = ws.cell(i, 1, t); c.font = f; c.alignment = Alignment(wrap_text=True, vertical='top')
ws.column_dimensions['A'].width = 160

# ---------------------------------------------------------------- Class_definition
rows = []
for d in ALL:
    I = INFO[d]; K = INFOK[d]; P = PROT[d]
    rows.append([d, *zelle(d), I['dokument'], UNIT.get(I['einheit'], I['einheit']), g4(P['median']), K['n'], K['aktiv'], K['anteil_kleinere'],
                 I['n'], I['aktiv'], I['nicht_aktiv'], I['anteil_kleinere'], P['geaendert'], len(P['ausgeschlossen']), ', '.join(P['ausgeschlossen']) or '-',
                 I['p_sauber'], I['gruppen'], I['redundanz_G_p']])
ws = sheet('Class_definition', ['Data set', 'Activity type', 'Diversity', 'Band (main analysis)', 'ChEMBL document', 'Unit', 'Median boundary',
                                 'n main', 'Active main', 'Minor share main', 'n median rule', 'Active', 'Inactive', 'Minor share',
                                 'Class changed', 'Excluded', 'Excluded compounds', 'p (cleaned)', 'G', 'G/p'], rows,
           [9, 13, 11, 13, 16, 8, 10, 7, 8, 9, 8, 7, 8, 9, 8, 8, 40, 9, 6, 8])
r = len(rows) + 3
tot_ch = sum(PROT[d]['geaendert'] for d in ALL); tot_n = sum(INFO[d]['n'] for d in ALL)
for t in [f'Class changed: compounds in both analyses whose class differs. Total {tot_ch} of {tot_n} compounds classified under the median rule.',
          'DS02, DS06 and DS12: many values lie exactly on the median (highest tested concentration or the same MIC), and values equal to the median are inactive, so the active share is below 50 %. In DS02 and DS06 the classes are identical to the main analysis.',
          'DS13: 6 compounds with ">" below the median are not decidable and are excluded. The data set was rebuilt (1,361 cleaned descriptors instead of 1,370, 58 groups instead of 66).']:
    r = text(ws, r, t)

# ---------------------------------------------------------------- Summary
ws = wb.create_sheet('Summary')
r = text(ws, 1, 'Tests of H1 to H4 under the median rule, next to the main analysis (16 series, descriptive, no decision, no Holm correction)', TF)


def res(h, x):
    if h == 'H1': return f"median difference {x['H1']['median_diff']:.4f}, positive in {x['H1']['positiv']}/16"
    if h == 'H1s': return ', '.join(f"{k}: median {v['median_diff']:.4f}, positive {v['positiv']}/16, p = {v['p']:.2e}" for k, v in x['H1']['sekundaer'].items())
    if h == 'H2': return f"L = {x['H2']['L']:.0f}, full order in {x['H2']['sekundaer']['ordnung_vollstaendig']}/16, median factors A1u {x['H2']['median_IF']['A1u']:.2f}, A1 {x['H2']['median_IF']['A1']:.2f}, A2 {x['H2']['median_IF']['A2']:.2f}"
    if h == 'H2s': return f"A2 > A1 p = {x['H2']['sekundaer']['A2_gt_A1']:.2e}, A1 > A1u p = {x['H2']['sekundaer']['A1_gt_A1u']:.2e}"
    if h == 'H3': return f"rho = {x['H3']['rho']:.3f}"
    if h == 'H3s': return ', '.join(f"{k}: rho = {v['rho']:.3f}, p = {v['p']:.3f}" for k, v in x['H3']['sekundaer'].items())
    if h == 'H4': return f"{x['H4']['anzahl']} of 16 ({', '.join(x['H4']['datensaetze']) or '-'}), binomial P(X >= {x['H4']['anzahl']}) = {x['H4']['binomial_p_ab_k']:.3f}"
    if h == 'H4s': return f"no selection: {len(x['H4']['sekundaer']['ohneAuswahl'])} ({', '.join(x['H4']['sekundaer']['ohneAuswahl']) or '-'}), random forest: {len(x['H4']['sekundaer']['RF'])} ({', '.join(x['H4']['sekundaer']['RF']) or '-'})"
    if h == 'H4e': return f"{x['H4']['median_beste_minus_kontrolle']:.4f}, 95 % interval {x['H4']['median_ci95'][0]:.4f} to {x['H4']['median_ci95'][1]:.4f}"


rows = [
    ['H1', 'Optimism by selection: apparent Q2 > nested Q2 under Route A1', 'Wilcoxon signed-rank, one-sided', res('H1', H), H['H1']['p'], res('H1', HK), HK['H1']['p']],
    ['H1 secondary', 'same for A1u and A2', 'Wilcoxon, one-sided', res('H1s', H), '', res('H1s', HK), ''],
    ['H2', 'Calibrated inflation factor ordered A2 > A1 > A1u', 'Page trend test, 16 blocks', res('H2', H), H['H2']['p'], res('H2', HK), HK['H2']['p']],
    ['H2 secondary', 'A2 > A1 and A1 > A1u', 'Wilcoxon, one-sided', res('H2s', H), '', res('H2s', HK), ''],
    ['H3', 'Smaller G/p, larger calibrated inflation factor (Route A2)', 'Spearman, one-sided negative', res('H3', H), H['H3']['p'], res('H3', HK), HK['H3']['p']],
    ['H3 secondary', 'same for A1 and A1u', 'Spearman, one-sided negative', res('H3s', H), '', res('H3s', HK), ''],
    ['H4', 'Series with a 95 % bootstrap interval (route - no filter) entirely above 0 (main analysis: refuted if >= 4 of 16)', 'count', res('H4', H), '', res('H4', HK), ''],
    ['H4 secondary', 'same count against no selection and against random forest', 'count', res('H4s', H), '', res('H4s', HK), ''],
    ['H4 effect', 'median of (best route - no filter), nested Q2', 'bootstrap over data sets, 20,000 draws', res('H4e', H), '', res('H4e', HK), ''],
]
r = table(ws, r, ['Hypothesis', 'Statement', 'Test', 'Median rule: result', 'Median rule: p (one-sided)', 'Main analysis: result', 'Main analysis: p (one-sided)'],
          rows, [14, 44, 30, 70, 14, 70, 14])
r = text(ws, r, 'Bootstrap over data sets for the H4 effect: default_rng([93, 0, 308]) under the median rule, default_rng([93, 0, 208]) in the main analysis (seed not specified, set analogously).')
r = text(ws, r, 'In the main analysis H1 and H2 were confirmed, H3 was not confirmed and H4 was not refuted (Holm, m = 4). These decisions are not repeated here.')
r += 1
r = text(ws, r, 'One row per data set (median rule)', HF)
rows = []
for d in ALL:
    t = T[d]; z = zelle(d)
    rows.append([d, *z, t['n'], t['aktiv'], t['p'], t['G'], t['Gp'],
                 *[t['app_' + q]['q2'] - t[q]['q2'] for q in ROUTES],
                 *[t[m]['q2'] for m in ['A1', 'A1u', 'A2', 'ohneAuswahl', 'ohneFilter', 'RF']],
                 *[t[m]['p_q2'] for m in ['A1', 'A1u', 'A2']],
                 *[t['T_' + q]['inf_kalibriert'] for q in ['A1u', 'A1', 'A2']],
                 ', '.join(q for q in ROUTES if t['bootstrap'][f'{q}-ohneFilter']['ganz_ueber_0']) or '-'])
table(ws, r, ['Data set', 'Activity type', 'Diversity', 'Band (main analysis)', 'n', 'Active', 'p (cleaned)', 'G', 'G/p',
              'Optimism A1 (apparent - nested Q2)', 'Optimism A1u', 'Optimism A2', 'Nested Q2 A1', 'Nested Q2 A1u', 'Nested Q2 A2',
              'Q2 no selection', 'Q2 no filter', 'Q2 random forest', 'p(Q2) A1', 'p(Q2) A1u', 'p(Q2) A2',
              'Calibrated inflation A1u', 'Calibrated inflation A1', 'Calibrated inflation A2', 'H4: routes with interval above 0'], rows)
ws.freeze_panes = 'A3'

# ---------------------------------------------------------------- Comparison_per_series
rows = []
for d in ALL:
    t, k = T[d], TK[d]
    h4 = lambda x: ', '.join(q for q in ROUTES if x['bootstrap'][f'{q}-ohneFilter']['ganz_ueber_0']) or '-'
    rows.append([d, PROT[d]['geaendert'], len(PROT[d]['ausgeschlossen']),
                 *[v for q in ROUTES for v in (k['app_' + q]['q2'] - k[q]['q2'], t['app_' + q]['q2'] - t[q]['q2'])],
                 *[v for q in ['A1u', 'A1', 'A2'] for v in (k['T_' + q]['inf_kalibriert'], t['T_' + q]['inf_kalibriert'])],
                 *[v for q in ROUTES + ['ohneFilter', 'RF'] for v in (k[q]['q2'], t[q]['q2'])],
                 *[v for q in ROUTES + ['ohneFilter', 'RF'] for v in (k[q]['auc'], t[q]['auc'])],
                 h4(k), h4(t)])
cols = ['Data set', 'Class changed', 'Excluded']
cols += [f'{lab} {q} {w}' for q in ROUTES for lab, w in [('Optimism', 'main'), ('Optimism', 'median')]]
cols += [f'Calibrated inflation {q} {w}' for q in ['A1u', 'A1', 'A2'] for w in ['main', 'median']]
QN = {'A1': 'A1', 'A1u': 'A1u', 'A2': 'A2', 'ohneFilter': 'no filter', 'RF': 'RF'}
cols += [f'Q2 {QN[q]} {w}' for q in ROUTES + ['ohneFilter', 'RF'] for w in ['main', 'median']]
cols += [f'AUC {QN[q]} {w}' for q in ROUTES + ['ohneFilter', 'RF'] for w in ['main', 'median']]
cols += ['H4 routes main', 'H4 routes median']
ws = sheet('Comparison_per_series', cols, rows, [9, 8, 8] + [10] * (len(cols) - 5) + [14, 14])
r = len(rows) + 3
r = text(ws, r, 'Optimism = apparent Q2 - nested Q2. Q2 and AUC of the routes, no filter and RF are nested leave-one-out (RF out-of-bag). H4 routes: routes whose 95 % bootstrap interval (route - no filter) lies entirely above 0.')
r = text(ws, r, 'DS02 and DS06: identical classes, results identical except for the bootstrap seed (208 against 308).')

# ---------------------------------------------------------------- All_results_long and DS sheets
long = []
for d in ALL:
    t = T[d]; z = zelle(d); base = [d, *z]
    for m in MODELS:
        for k in ['q2', 'auc', 'p_q2', 'p_auc', 'n_perm', 'sens', 'spez', 'ba', 'mcc']:
            long.append(base + ['model', MNAME[m], {'spez': 'spec', 'q2': 'Q2', 'auc': 'AUC', 'p_q2': 'p(Q2)', 'p_auc': 'p(AUC)', 'n_perm': 'permutations', 'ba': 'balanced accuracy', 'mcc': 'MCC', 'sens': 'sens'}[k], t[m][k]])
    for q in ROUTES:
        for k, kn in [('beobachtet', 'observed hits'), ('null_mittel', 'null mean'), ('null_sd', 'null SD'), ('tests', 'number of tests m'),
                      ('null_rate', 'null hit rate q'), ('p', 'p'), ('inf_nominal', 'inflation nominal (5 %)'), ('inf_kalibriert', 'inflation calibrated'),
                      ('perzentil_97_5', 'null 97.5th percentile')]:
            long.append(base + ['hit count', f'Route {q}', kn, t['T_' + q][k]])
    for key, v in t['bootstrap'].items():
        q, k = key.split('-')
        for kk, kn in [('diff', 'difference nested Q2'), ('lo', '2.5 %'), ('hi', '97.5 %'), ('ganz_ueber_0', 'interval above 0')]:
            long.append(base + ['paired bootstrap', f'Route {q} - {KNAME[k]}', kn, v[kk]])
sheet('All_results_long', ['Data set', 'Activity type', 'Diversity', 'Band (main analysis)', 'Block', 'Model / route / comparison', 'Quantity', 'Value'],
      long, [9, 16, 14, 14, 16, 48, 26, 14])

for d in ALL:
    I = INFO[d]; t = T[d]; z = zelle(d); c = BYASSAY.get(I['assay'], {}); P = PROT[d]
    ws = wb.create_sheet(d)
    r = text(ws, 1, f"{d} – {' / '.join(z)} – median rule", TF)
    desc = [['ChEMBL document', I['dokument']], ['DOI', I.get('doi') or c.get('doi')], ['Journal, year', f"{c.get('journal')}, {c.get('jahr')}"],
            ['Title', c.get('titel')], ['Assay', I['assay']], ['Assay description', c.get('beschreibung')], ['Target / organism / cell line', I.get('ziel')],
            ['Unit', UNIT.get(I['einheit'], I['einheit'])], ['Class boundary', f"median rule: active if value < {g4(P['median'])} {UNIT.get(I['einheit'], I['einheit'])}"],
            ['n', I['n']], ['Active / inactive', f"{I['aktiv']} / {I['nicht_aktiv']}"], ['Share of minor class', I.get('anteil_kleinere')],
            ['Class changed against main analysis', P['geaendert']], ['Excluded (censored, not decidable)', ', '.join(P['ausgeschlossen']) or 'none'],
            ['Descriptors raw / cleaned', f"{I.get('p_roh')} / {I.get('p_sauber')}"],
            ['Correlation groups G, G/p', f"{I.get('gruppen')}, {I.get('redundanz_G_p')}"], ['Top generic scaffold share', round(STRUK[d]['top'], 4)],
            ['Mean pairwise Tanimoto', round(STRUK[d]['mts'], 4)],
            ['Compounds removed by the supplementary rule (main analysis)', ', '.join(I.get('entfernt') or []) or 'none']]
    r = table(ws, r, ['Item', 'Value'], desc, [44, 70])
    r = text(ws, r, 'Models (leave-one-out, random forest out-of-bag)', HF)
    r = table(ws, r, ['Model', 'Q2', 'AUC', 'p(Q2)', 'p(AUC)', 'Permutations', 'Sensitivity', 'Specificity', 'Balanced accuracy', 'MCC'],
              [[MNAME[m], t[m]['q2'], t[m]['auc'], t[m]['p_q2'], t[m]['p_auc'], t[m]['n_perm'], t[m]['sens'], t[m]['spez'], t[m]['ba'], t[m]['mcc']] for m in MODELS],
              [44, 12, 12, 12, 12, 12, 12, 12, 12, 12])
    r = text(ws, r, 'Univariate hit count against its label-permutation null (1,000 permutations)', HF)
    r = table(ws, r, ['Route', 'Observed hits', 'Null mean', 'Null SD', 'Tests m', 'Null hit rate q', 'p', 'Inflation nominal (5 %)', 'Inflation calibrated', 'Null 97.5th percentile'],
              [[f'Route {q}', *[t['T_' + q][k] for k in ['beobachtet', 'null_mittel', 'null_sd', 'tests', 'null_rate', 'p', 'inf_nominal', 'inf_kalibriert', 'perzentil_97_5']]] for q in ROUTES])
    r = text(ws, r, 'Paired stratified bootstrap of nested Q2 (20,000 draws, seed 308)', HF)
    r = table(ws, r, ['Comparison', 'Difference', '2.5 %', '97.5 %', 'Interval entirely above 0'],
              [[f"Route {k.split('-')[0]} - {KNAME[k.split('-')[1]]}", v['diff'], v['lo'], v['hi'], 'yes' if v['ganz_ueber_0'] else 'no'] for k, v in t['bootstrap'].items()])
    for row in ws.iter_rows():
        for cc in row:
            if isinstance(cc.value, float): cc.number_format = '0.0000'

# ---------------------------------------------------------------- Robustness
ws = wb.create_sheet('Robustness')
r = text(ws, 1, 'Robustness checks as in the main analysis (descriptive, no decision)', TF)
r = text(ws, r, f"R1: series with at least one PLS prediction > 2 or < -1 in any model. Excluded under the median rule: {', '.join(R['ausgeschlossen_R1']) or '-'} (main analysis: {', '.join(RK['ausgeschlossen_R1'])}).")
r = text(ws, r, f"R2 (label-free): series in which a held-out compound reaches |z| > 100 on at least one descriptor. Excluded under the median rule: {', '.join(R['ausgeschlossen_R2']) or '-'} (main analysis: {', '.join(RK['ausgeschlossen_R2'])}).")
r = text(ws, r, 'R2 does not depend on the classes. It differs from the main analysis only where the compounds differ (DS13, 6 compounds excluded).')
r = text(ws, r, 'H4 limit of the main analysis for fewer series: smallest k with P(X >= k) < 0.05 under Binomial(n, 0.075). Shown for orientation only.')
rows = []
for key, lab in [('alle16', 'all 16 series'), ('R1', 'without R1 series'), ('R2', 'without R2 series')]:
    X = R[key]
    rows.append([lab, X['n'], f"{X['H1']['A1']['positiv']}/{X['n']}", X['H1']['A1']['p'], X['H1']['A1u']['p'], X['H1']['A2']['p'], X['H2']['p'],
                 f"{X['H2']['ordnung_vollstaendig']}/{X['n']}", X['H3']['rho'], X['H3']['p'], f"{X['H4']['anzahl']} ({', '.join(X['H4']['datensaetze']) or '-'})",
                 X['H4']['grenze_hauptanalyse'], X['H4']['median_beste_minus_kontrolle']])
r = table(ws, r + 1, ['Series', 'n', 'H1 A1 positive', 'H1 p A1', 'H1 p A1u', 'H1 p A2', 'H2 p', 'H2 full order', 'H3 rho', 'H3 p', 'H4 count (series)',
                      'H4 limit (main analysis)', 'H4 median best - control'], rows, [20, 5, 10, 10, 10, 10, 10, 10, 10, 10, 26, 10, 12])
r = text(ws, r, 'AUC instead of Q2, all 16 series', HF)
r = table(ws, r, ['Route', 'H1: apparent - nested AUC, positive', 'Median', 'Wilcoxon p (one-sided)'],
          [[f'Route {q}', f"{v['positiv']}/16", v['median'], v['p']] for q, v in R['AUC_H1'].items()])
a4 = R['AUC_H4']
r = text(ws, r, f"H4 with AUC (same bootstrap draws): {a4['anzahl']} of 16 ({', '.join(a4['datensaetze']) or '-'}). Median best route - control {a4['median_beste_minus_kontrolle']:.4f}. Main analysis: {RK['AUC_H4']['anzahl']} of 16 ({', '.join(RK['AUC_H4']['datensaetze']) or '-'}).")
r = table(ws, r + 1, ['Series and route', 'AUC difference (route - no filter)', '2.5 %', '97.5 %'], [[k, v['diff'], v['lo'], v['hi']] for k, v in a4['intervalle'].items()])

# ---------------------------------------------------------------- Predictions
ws = wb.create_sheet('Predictions')
r = text(ws, 1, 'Leave-one-out predictions of all PLS models and out-of-bag prediction of the random forest, every compound (median rule)', TF)
r = text(ws, r, 'Orange: at least one PLS prediction above 2 or below -1 (R1). Yellow: no such prediction, but label-free |z| above 100 (R2). Class main analysis: class of the same compound in the main analysis.')
head = ['Data set', 'ChEMBL ID', 'SMILES', 'Relation', 'Value', 'Unit', 'Class median rule', 'Class main analysis', 'Largest label-free |z|'] + [MNAME[m] for m in MODELS]
for j, h in enumerate(head, 1):
    x = ws.cell(r + 1, j, h); x.font = HF; x.fill = FILL; x.alignment = Alignment(wrap_text=True, vertical='top')
r += 2; r_head = r - 1
for d in ALL:
    tab = pd.read_csv(f'{DSD}/{d}_y.csv'); assert list(tab.y) == V[d]['y']
    for i in range(len(tab)):
        preds = [V[d]['pred'][m][i] for m in MODELS]
        unst = any(p > 2 or p < -1 for m, p in zip(MODELS, preds) if m != 'RF')
        row = [d, tab.chembl_id[i], tab.smiles[i], tab.relation[i], g4(tab.wert[i]), UNIT.get(INFO[d]['einheit'], INFO[d]['einheit']),
               'active' if tab.y[i] == 1 else 'inactive', 'active' if tab.y_ersatzregel[i] == 1 else 'inactive', round(Z[d][i], 1)] + [round(p, 4) for p in preds]
        for j, v in enumerate(row, 1):
            x = ws.cell(r, j, v); x.font = NF
            if unst: x.fill = ORA
            elif Z[d][i] > 100: x.fill = YEL
        r += 1
ws.freeze_panes = ws.cell(r_head + 1, 3).coordinate
ws.auto_filter.ref = f'A{r_head}:{get_column_letter(len(head))}{r - 1}'
for j, w in enumerate([9, 15, 40, 8, 10, 8, 10, 10, 10] + [12] * len(MODELS), 1):
    ws.column_dimensions[get_column_letter(j)].width = w

# ---------------------------------------------------------------- Exploratory
ws = wb.create_sheet('Exploratory')
r = text(ws, 1, 'Exploratory groupings as in the main analysis (descriptive, no test claim)', TF)


def grp(keyfun, label):
    global r
    r = text(ws, r, label, HF)
    g = {}
    for d in PRIM: g.setdefault(keyfun(d), []).append(d)
    rows = []
    for k, ds in sorted(g.items()):
        f = lambda fn: float(np.median([fn(T[x]) for x in ds]))
        rows.append([k, len(ds), ', '.join(ds), f(lambda t: t['app_A1']['q2'] - t['A1']['q2']), f(lambda t: t['app_A1u']['q2'] - t['A1u']['q2']),
                     f(lambda t: t['app_A2']['q2'] - t['A2']['q2']), f(lambda t: t['T_A1u']['inf_kalibriert']), f(lambda t: t['T_A1']['inf_kalibriert']),
                     f(lambda t: t['T_A2']['inf_kalibriert']), f(lambda t: t['A1']['q2']), f(lambda t: max(t[q]['q2'] for q in ROUTES) - t['ohneFilter']['q2']),
                     sum(any(T[x]['bootstrap'][f'{q}-ohneFilter']['ganz_ueber_0'] for q in ROUTES) for x in ds)])
    r = table(ws, r, ['Group', 'Series', 'Members', 'Median optimism A1', 'Median optimism A1u', 'Median optimism A2', 'Median calibrated inflation A1u',
                      'Median calibrated inflation A1', 'Median calibrated inflation A2', 'Median nested Q2 A1', 'Median best route - no filter', 'H4 counts'],
              rows, [24, 7, 60, 12, 12, 12, 12, 12, 12, 12, 12, 9])


grp(lambda d: zelle(d)[2], '1. Band of the minor class in the main analysis (under the median rule all series have 38-50 %)')
grp(lambda d: zelle(d)[1], '2. Structural diversity')
grp(lambda d: zelle(d)[0], '3. Activity type')
r = text(ws, r, '4. Does an excess of A1u hits over the label-free null coincide with a nominally significant nested Q2? (criterion as in the main analysis)', HF)
rows = []; ct = Counter()
for d in PRIM:
    t = T[d]; exc = t['T_A1u']['beobachtet'] > t['T_A1u']['perzentil_97_5']; sig = t['A1u']['p_q2'] <= 0.05
    ct[(exc, sig)] += 1
    rows.append([d, t['T_A1u']['beobachtet'], t['T_A1u']['perzentil_97_5'], t['T_A1u']['p'], 'yes' if exc else 'no', t['A1u']['q2'], t['A1u']['p_q2'], 'yes' if sig else 'no'])
r = table(ws, r, ['Series', 'A1u observed hits', 'Null 97.5th percentile', 'p (hit count)', 'Excess (observed > 97.5th percentile)', 'Nested Q2 A1u', 'p(Q2) A1u', 'p(Q2) <= 0.05'], rows)
r = table(ws, r, ['', 'Nested Q2 A1u p <= 0.05', 'p > 0.05'], [['Hit excess', ct[(True, True)], ct[(True, False)]], ['No hit excess', ct[(False, True)], ct[(False, False)]]])

# ---------------------------------------------------------------- Deviations
DEV = [
    ['30.10.2026', 'The bootstrap over data sets (interval of the H4 effect estimate) uses default_rng([93, 0, 308]). The addendum names bootstrap seed 308 but not this seed explicitly.', 'Gap in the addendum, set analogously to the main analysis (208).', 'Affects only the interval of the H4 effect estimate.'],
    ['30.10.2026', 'Exploratory groupings and the Predictions sheet were added to the workbook. The addendum names the same key figures and tests as for H1 to H4, R1, R2, AUC and the per-series comparison.', 'Same layout as the main analysis workbook.', 'None on results. Descriptive only.'],
    ['note', 'The run was interrupted once by a restart of the compute container (27.10.2026, block 16) and resumed on 30.10.2026. The script writes a package only when it is complete and skips finished packages on restart.', 'Technical.', 'None. All 2,320 packages complete, indices without gaps, 6 randomly chosen permutations recomputed and identical.'],
    ['note', 'The four size subsets DS21 to DS24 were not computed under the median rule.', 'As specified in the addendum: the pre-specification names the median rule only for the 16 data sets.', '-'],
]
ws = sheet('Deviations', ['Date', 'Deviation', 'Reason', 'Effect'], DEV, [12, 100, 50, 60])
for row in ws.iter_rows(min_row=2):
    for cc in row: cc.alignment = Alignment(wrap_text=True, vertical='top')

order = ['README', 'Class_definition', 'Summary', 'Comparison_per_series', 'All_results_long'] + ALL + ['Robustness', 'Predictions', 'Exploratory', 'Deviations']
wb._sheets = [wb[n] for n in order]
wb.save(OUT)
print('saved', OUT, len(wb.sheetnames), 'sheets, long', len(long))
