"""Results workbook of the ChEMBL benchmark following the reporting duties of the pre-specification
(SHA-256 a13bb2c0...ae419): sheets README, Selection, Summary, All_results_long, DS01-DS16, DS21-DS24,
plus Robustness, Exploratory, Deviations (and for the main analysis Comparison and Excluded_compounds).
Usage: python3 ergebnismappe.py haupt|vorab <output file>   (haupt = main analysis, vorab = analysis under the original selection rule)
Nothing is recomputed except descriptive key figures (scaffold share, mean Tanimoto similarity),
which are determined from the SMILES of the data sets, and the exploratory summaries.
Some cell texts are written in German and translated by revision2/tables/translate_workbooks.py into the
published Additional files 12, 17 and 18.
"""
import json, sys, importlib.util
from collections import Counter
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog('rdApp.*')

ANALYSE, OUT = sys.argv[1], sys.argv[2]
HAUPT = ANALYSE == 'haupt'
ERG = 'ergebnisse_k' if HAUPT else 'ergebnisse'
DSD = 'datasets_k' if HAUPT else 'datasets'
A = json.load(open(f'{ERG}/auswertung.json')); T = A['tabelle']; H = A['hypothesen']
R = json.load(open(f'{ERG}/robustheit.json'))
CAND = json.load(open('cand_ranked.json')); BYASSAY = {x['assay']: x for x in CAND}
VERD = json.load(open('verdikte.json'))
PRIM = [f'DS{i:02d}' for i in range(1, 17)]; TEIL = ['DS21', 'DS22', 'DS23', 'DS24']; ALL = PRIM + TEIL
INFO = {d: json.load(open(f'{DSD}/{d}_info.json'))['info'] for d in ALL}
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
F = 'Arial'; HF = Font(name=F, bold=True); NF = Font(name=F); TF = Font(name=F, bold=True, size=13)
FILL = PatternFill('solid', fgColor='D9E1F2'); GRN = PatternFill('solid', fgColor='C6EFCE'); RED = PatternFill('solid', fgColor='FFC7CE')
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
    return (ART[z[0]], DIV[z[1]], BAND[z[2]]) if z[0] != 'Teilmenge' else ('Subset', ART[z[1]], z[2])


# descriptive structural key figures from the SMILES (as in the selection: generic Bemis-Murcko scaffold, Morgan r2 2048)
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
STRUK = {}
for d in ALL:
    tab = pd.read_csv(f'{DSD}/{d}_y.csv'); ms = [Chem.MolFromSmiles(s) for s in tab.smiles]
    sc = [Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m))) for m in ms]   # as candidates.py / auswahl_neu.py
    fps = [GEN.GetFingerprint(m) for m in ms]
    sims = [s for i in range(len(fps)) for s in DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])]
    STRUK[d] = dict(top=Counter(sc).most_common(1)[0][1] / len(sc), mts=float(np.mean(sims)), n_scaf=len(set(sc)))


def schwelle(dok, art):
    v = VERD.get(dok, {}); s = v.get('schwelle', '')
    if s and s not in ('keine', 'offen', ''):
        return 'stage 1 (authors)', s
    return 'stage 2 (fallback rule)', {'MIC': 'active if MIC <= 8 ug/mL', 'Enzym': 'active if IC50 <= 1 uM (pChEMBL >= 6)',
                                       'Rezeptor': 'active if Ki <= 100 nM (pChEMBL >= 7)', 'Zytotox': 'active if IC50/GI50/CC50 <= 10 uM (pChEMBL >= 5)'}[art]


# ---------------------------------------------------------------- README
ws = wb.active; ws.title = 'README'
titel = ('MAIN ANALYSIS (supplementary rule excluding reference and control compounds)' if HAUPT
         else 'ANALYSIS UNDER THE ORIGINAL PRE-SPECIFICATION ONLY (reported in full, superseded as main analysis)')
L = [(f'ChEMBL benchmark – results workbook – {titel}', TF), ('', NF),
     ('Manuscript: Werle J, Rondevaldova J, Werle E. Honest performance estimation for small-sample QSAR (Journal of Cheminformatics, revision 2).', NF),
     ('Pre-specification frozen 22.10.2026 18:10, SHA-256 a13bb2c0b513289b160de744bd8d6bc47662d46553ab5a2fe3db6d336a0ae419.', NF)]
if HAUPT:
    L += [('Supplementary specification frozen 24.10.2026 17:20, before any model of this analysis was fitted, SHA-256 36668775daafd821c1871a15c6a697ea5e3c25afebbb3d6b6a4db818c902df87.', NF),
          ('It adds one rule that the pre-specification lacked: reference drugs, clinical candidates and tool compounds of other classes measured for comparison or as controls, and values taken from other publications, do not belong to the investigated series. The data set selection was repeated mechanically with this rule. Two of the 16 series were replaced (cells MIC/diverse/B and Receptor/diverse/A), and so was the MIC parent pool of the size subsets. See sheets Deviations, Excluded_compounds and Comparison.', NF)]
    L += [('This is the main analysis of the manuscript. The analysis under the original rule is reported in full in a separate workbook with identical layout.', HF)]
else:
    L += [('This workbook reports the analysis under the original pre-specification. A gap in the selection rule (reference and control compounds were not excluded) was found after these results were known. The main analysis applies a supplementary rule and is reported in a separate workbook with identical layout. All four hypothesis decisions are the same in both analyses.', HF)]
L += [('', NF), ('Data source: ChEMBL 37 (release 01.05.2026). One data set = one publication and one assay.', NF),
      ('Models: PLS-DA with one latent variable, leave-one-out, autoscaling on the training compounds of each fold, columns with training SD < 1e-6 dropped in the fold. Univariate filter Mann-Whitney p <= 0.05 and |r_rb| > 0.3, recomputed in every fold. Correlation grouping |Pearson| or |Spearman| > 0.7, most-connected-first, ties by largest variance. Folds without descriptor or with one class predict 0.5. Random forest 300 trees, min_samples_leaf = 1, OOB, Morgan fingerprints radius 2, 2,048 bit.', NF),
      ('Permutations: nested models (A1, A1u, A2, no selection, no filter) 2,000, apparent models 1,000, hit counts 1,000, random forest 500. Permutation p = (1 + #{null >= observed - 1e-9}) / (1 + N).', NF),
      ('Inflation factor: SD of the null hit count (ddof = 1) divided by the binomial SD, nominal at 5 % and calibrated at the own null hit rate q = mean(null)/m. m = all cleaned descriptors for A2, number of correlation groups G for A1 and A1u.', NF),
      ('Paired comparisons: paired bootstrap of the fixed LOO/OOB predictions, stratified by class, 20,000 draws, same draws for all models of a data set, 95 % percentile interval.', NF),
      ('Seeds: numpy default_rng([93, data set number, analysis code' + (' + 200' if HAUPT else '') + ', permutation index]). Analysis codes A1 1, A1u 2, A2 3, no selection 4, no filter 5, RF 6, hit count 7, bootstrap 8, subset draw 9, apparent A1/A1u/A2 11/12/13.', NF),
      ('Hypotheses H1 to H4 on the 16 natural series, one-sided, alpha = 0.05, Holm over H1 to H4. The four size subsets DS21 to DS24 are exploratory and do not enter the tests.', NF), ('', NF),
      ('Sheets', HF),
      ('Selection – cell assignment of the 20 data sets, excluded publications, full selection log (every examined assay with decision), all candidates.', NF),
      ('Summary – hypothesis tests H1 to H4 and one row per data set.', NF),
      ('All_results_long – every number of this workbook in long format (one value per row).', NF),
      ('DS01 to DS16, DS21 to DS24 – one sheet per data set, identical layout: description, models, hit counts, paired bootstrap.', NF),
      ('Robustness – H1 to H4 without series with numerically unstable predictions (R1), without series failing a label-free criterion (R2), and with AUC instead of Q2.', NF),
      ('Unstable_predictions – every compound of the series with numerically unstable predictions (R1), with predictions of all PLS models, label-free |z| and the reason.', NF),
      ('Exploratory – descriptive analyses named in the pre-specification (band, diversity, activity type, size, A1u hit excess, origin of the class boundary).', NF),
      ('Deviations – every deviation from the specification with date, reason and effect.', NF)]
if HAUPT:
    L += [('Excluded_compounds – compounds removed by the supplementary rule.', NF),
          ('Comparison – main analysis against the analysis under the original rule, hypothesis by hypothesis.', NF)]
L += [('', NF), ('SMILES, ChEMBL IDs, raw values, classes, descriptor matrices, raw permutation outputs and all scripts are in the code archive.', NF),
      ('Column names in German in the code archive: aktiv = active, nicht_aktiv = inactive, ohneAuswahl = no selection, ohneFilter = no filter, Treffer = hit count, scheinbar = apparent.', NF)]
for i, (t, f) in enumerate(L, 1):
    c = ws.cell(i, 1, t); c.font = f; c.alignment = Alignment(wrap_text=True, vertical='top')
ws.column_dimensions['A'].width = 160

# ---------------------------------------------------------------- Selection
rows = []
for d in ALL:
    I = INFO[d]; c = BYASSAY.get(I['assay'], {}); z = zelle(d)
    st, regel = schwelle(I['dokument'], I['zelle'][0] if I['zelle'][0] != 'Teilmenge' else I['zelle'][1])
    rows.append([d, *z, I['dokument'], I.get('doi') or c.get('doi'), c.get('journal'), c.get('jahr'), I['assay'], I.get('ziel'),
                 c.get('standard_type'), I.get('einheit'), I['n'], I['aktiv'], I['nicht_aktiv'], I.get('anteil_kleinere'),
                 round(STRUK[d]['top'], 4), round(STRUK[d]['mts'], 4), st, regel,
                 len(I.get('entfernt') or []) if HAUPT else '', I.get('eltern_n') or '', str(I.get('seed') or '')])
ws = wb.create_sheet('Selection')
r = text(ws, 1, 'A. The 20 data sets (16 natural series in the cells of the design, 4 size subsets)', HF)
r = table(ws, r, ['Data set', 'Activity type', 'Diversity', 'Band of minor class', 'ChEMBL document', 'DOI', 'Journal', 'Year', 'Assay',
                  'Target / organism / cell line', 'Activity type (ChEMBL)', 'Unit', 'n', 'Active', 'Inactive', 'Share of minor class',
                  'Top generic scaffold share', 'Mean pairwise Tanimoto', 'Class boundary', 'Rule'] +
          (['Compounds removed (supplementary rule)'] if HAUPT else ['']) + ['Parent n (subsets)', 'Seed (subsets)'], rows,
          [10, 14, 12, 14, 16, 28, 18, 6, 14, 30, 10, 10, 6, 7, 8, 9, 10, 10, 20, 34, 12, 10, 18])
r = text(ws, r, 'Diversity: congeneric = top generic Bemis-Murcko scaffold >= 50 % of compounds, diverse = top scaffold < 25 % and mean Tanimoto < 0.30 (Morgan r2, 2,048 bit). Values recomputed from the SMILES of the final data set.')
r = text(ws, r, 'Class boundary: none of the publications of the 16 natural series states an activity boundary, so the fallback rule of the specification applies to all 16 (stage 2). The MIC parent pool of the main analysis (Spink et al. 2015) states "Active compounds (MIC <= 8 ug/mL)" (stage 1, identical to the fallback rule).')
r += 1
r = text(ws, r, 'B. Publications excluded on reading (review, values compiled from other publications, compounds preselected by activity, method papers)', HF)
TXT = {'review': 'review', 'keine_serie': 'values compiled from other publications', 'vorselektiert': 'compounds preselected by activity (screening hits)',
       'methoden': 'method or assay paper without own series', 'unlesbar': 'full text not available'}
ex = [[d, i.get('doi'), TXT.get(i['art'], i['art']), i.get('rolle', '')] for d, i in sorted(VERD.items()) if i['art'] != 'original']
r = table(ws, r, ['ChEMBL document', 'DOI', 'Reason', 'Note'], ex)
ws.freeze_panes = 'A3'
# selection log
if HAUPT:
    prot = json.load(open('kontrollen/auswahl_neu.json'))['protokoll']
else:
    spec = importlib.util.spec_from_file_location('sel', 'select.py'); sel = importlib.util.module_from_spec(spec); spec.loader.exec_module(sel)
    _, prot, _, _ = sel.run(CAND, VERD)
DE = [('GEWAEHLT', 'SELECTED'), ('Klassengrenze verfehlt', 'class boundary missed'), ('Diversitaet dazwischen', 'diversity between congeneric and diverse'),
      ('Band schon besetzt', 'band already occupied'), ('Ziel schon belegt', 'target already used'), ('ausserhalb', 'outside'),
      ('Groesse', 'size'), ('kleinere Klasse', 'minor class'), ('Tanimoto', 'Tanimoto'), ('top', 'top scaffold')]
def en(s):
    for a, b in DE: s = s.replace(a, b)
    return s
sheet('Selection_log', ['ChEMBL document', 'Assay', 'Decision (in hash order)'], [[p[0], p[1] or '', en(p[2])] for p in prot], [18, 16, 90])
sheet('Candidates', ['Hash (first 16)', 'Activity type', 'Assay', 'ChEMBL document', 'DOI', 'Target', 'Organism', 'n', 'Active (provisional)',
                     'Inactive (provisional)', 'Share minor class (provisional)', 'Diversity', 'Exclusion reason (automatic stage)'],
      [[x['hash'][:16], ART[x['art']], x['assay'], x['dokument'], x.get('doi'), x['ziel'], x['organismus'], x['n_verbindungen'], x['aktiv_vorl'],
        x['nichtaktiv_vorl'], x['anteil_kleinere_klasse_vorl'], DIV.get(x['diversitaet'], x['diversitaet']), x['grund']] for x in CAND],
      [16, 12, 14, 16, 26, 26, 20, 6, 8, 8, 10, 12, 30])

# ---------------------------------------------------------------- Summary
ws = wb.create_sheet('Summary')
r = text(ws, 1, 'Hypothesis tests on the 16 natural series' + (' (main analysis)' if HAUPT else ' (original rule)'), TF)
h1, h2, h3, h4 = H['H1'], H['H2'], H['H3'], H['H4']
holm = H['holm_m4']
rows = [
    ['H1', 'Optimism by selection: apparent Q2 > nested Q2 under Route A1', 'Wilcoxon signed-rank, one-sided', f"median difference {h1['median_diff']:.4f}, positive in {h1['positiv']}/16", h1['p'], 'confirmed' if holm['H1'] else 'not confirmed'],
    ['H1 secondary', 'same for A1u and A2', 'Wilcoxon, one-sided', '; '.join(f"{k}: median {v['median_diff']:.4f}, positive {v['positiv']}/16, p = {v['p']:.2e}" for k, v in h1['sekundaer'].items()).replace(';', ','), '', 'descriptive'],
    ['H2', 'Calibrated inflation factor ordered A2 > A1 > A1u', 'Page trend test, 16 blocks', f"L = {h2['L']:.0f}, full order in {h2['sekundaer']['ordnung_vollstaendig']}/16, median factors A1u {h2['median_IF']['A1u']:.2f}, A1 {h2['median_IF']['A1']:.2f}, A2 {h2['median_IF']['A2']:.2f}", h2['p'], 'confirmed' if holm['H2'] else 'not confirmed'],
    ['H2 secondary', 'A2 > A1 and A1 > A1u', 'Wilcoxon, one-sided', f"A2 > A1 p = {h2['sekundaer']['A2_gt_A1']:.2e}, A1 > A1u p = {h2['sekundaer']['A1_gt_A1u']:.2e}", '', 'descriptive'],
    ['H3', 'Smaller G/p, larger calibrated inflation factor (Route A2)', 'Spearman, one-sided negative', f"rho = {h3['rho']:.3f}", h3['p'], 'confirmed' if holm['H3'] else 'not confirmed'],
    ['H3 secondary', 'same for A1 and A1u', 'Spearman, one-sided negative', ', '.join(f"{k}: rho = {v['rho']:.3f}, p = {v['p']:.3f}" for k, v in h3['sekundaer'].items()), '', 'descriptive'],
    ['H4', 'No supervised route beats the no-filter control beyond sampling uncertainty', 'count of series with a 95 % bootstrap interval (route - no filter) entirely above 0, refuted if >= 4 of 16',
     f"{h4['anzahl']} of 16 ({', '.join(h4['datensaetze'])}), binomial P(X >= {h4['anzahl']}) = {h4['binomial_p_ab_k']:.3f}", '', 'refuted' if h4['widerlegt'] else 'not refuted'],
    ['H4 secondary', 'same count against no selection and against random forest', 'count', f"no selection: {len(h4['sekundaer']['ohneAuswahl'])} ({', '.join(h4['sekundaer']['ohneAuswahl'])}), random forest: {len(h4['sekundaer']['RF'])} ({', '.join(h4['sekundaer']['RF'])})", '', 'descriptive'],
    ['H4 effect', 'median of (best route - no filter), nested Q2', 'bootstrap over data sets, 20,000 draws', f"{h4['median_beste_minus_kontrolle']:.4f}, 95 % interval {h4['median_ci95'][0]:.4f} to {h4['median_ci95'][1]:.4f}", '', 'descriptive'],
]
r = table(ws, r, ['Hypothesis', 'Statement', 'Test', 'Result', 'p (one-sided)', 'Decision (Holm, m = 4)'], rows, [14, 48, 34, 70, 14, 20])
fmt = lambda d: ', '.join(f"{k} {'confirmed' if v else 'not confirmed'}" for k, v in sorted(d.items()))
r = text(ws, r, f"Holm over H1 to H3 (H4 has its own counting criterion): m = 3: {fmt(H['holm_m3'])}. m = 4: {fmt(H['holm_m4'])}.")
r = text(ws, r, 'Bootstrap over data sets for the H4 effect: seed not specified in the pre-specification, used default_rng([93, 0, ' + ('208' if HAUPT else '8') + ']).')
r += 1
r = text(ws, r, 'One row per data set', HF)
rows = []
for d in ALL:
    t = T[d]; z = zelle(d)
    rows.append([d, *z, t['n'], t['p'], t['G'], t['Gp'],
                 *[t['app_' + q]['q2'] - t[q]['q2'] for q in ROUTES],
                 *[t[m]['q2'] for m in ['A1', 'A1u', 'A2', 'ohneAuswahl', 'ohneFilter', 'RF']],
                 *[t[m]['p_q2'] for m in ['A1', 'A1u', 'A2']],
                 *[t['T_' + q]['inf_kalibriert'] for q in ['A1u', 'A1', 'A2']],
                 ', '.join(q for q in ROUTES if t['bootstrap'][f'{q}-ohneFilter']['ganz_ueber_0']) or '-'])
table(ws, r, ['Data set', 'Activity type / subset', 'Diversity / type', 'Band / size', 'n', 'p (cleaned)', 'G', 'G/p',
              'Optimism A1 (apparent - nested Q2)', 'Optimism A1u', 'Optimism A2', 'Nested Q2 A1', 'Nested Q2 A1u', 'Nested Q2 A2',
              'Q2 no selection', 'Q2 no filter', 'Q2 random forest', 'p(Q2) A1', 'p(Q2) A1u', 'p(Q2) A2',
              'Calibrated inflation A1u', 'Calibrated inflation A1', 'Calibrated inflation A2', 'H4: routes with interval above 0'], rows)
ws.freeze_panes = 'A3'

# ---------------------------------------------------------------- All_results_long and DS sheets
long = []
for d in ALL:
    t = T[d]; z = zelle(d)
    base = [d, *z]
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
sheet('All_results_long', ['Data set', 'Activity type / subset', 'Diversity / type', 'Band / size', 'Block', 'Model / route / comparison', 'Quantity', 'Value'],
      long, [9, 16, 14, 14, 16, 48, 26, 14])

for d in ALL:
    I = INFO[d]; t = T[d]; z = zelle(d); c = BYASSAY.get(I['assay'], {})
    ws = wb.create_sheet(d)
    r = text(ws, 1, f"{d} – {' / '.join(z)}", TF)
    st, regel = schwelle(I['dokument'], I['zelle'][0] if I['zelle'][0] != 'Teilmenge' else I['zelle'][1])
    desc = [['ChEMBL document', I['dokument']], ['DOI', I.get('doi') or c.get('doi')], ['Journal, year', f"{c.get('journal')}, {c.get('jahr')}"],
            ['Title', c.get('titel')], ['Assay', I['assay']], ['Assay description', c.get('beschreibung')], ['Target / organism / cell line', I.get('ziel')],
            ['Unit', I.get('einheit')], ['Class boundary', f'{st}: {regel}'], ['n', I['n']], ['Active / inactive', f"{I['aktiv']} / {I['nicht_aktiv']}"],
            ['Share of minor class', I.get('anteil_kleinere')], ['Descriptors raw / cleaned', f"{I.get('p_roh')} / {I.get('p_sauber')}"],
            ['Correlation groups G, G/p', f"{I.get('gruppen')}, {I.get('redundanz_G_p')}"], ['Top generic scaffold share', round(STRUK[d]['top'], 4)],
            ['Mean pairwise Tanimoto', round(STRUK[d]['mts'], 4)]]
    if HAUPT:
        desc.append(['Compounds removed by the supplementary rule', ', '.join(I.get('entfernt') or []) or 'none'])
    if d in TEIL:
        desc.append(['Parent pool n, seed of the draw', f"{I.get('eltern_n')}, {I.get('seed')}"])
    r = table(ws, r, ['Item', 'Value'], desc, [44, 70])
    r = text(ws, r, 'Models (leave-one-out, random forest out-of-bag)', HF)
    r = table(ws, r, ['Model', 'Q2', 'AUC', 'p(Q2)', 'p(AUC)', 'Permutations', 'Sensitivity', 'Specificity', 'Balanced accuracy', 'MCC'],
              [[MNAME[m], t[m]['q2'], t[m]['auc'], t[m]['p_q2'], t[m]['p_auc'], t[m]['n_perm'], t[m]['sens'], t[m]['spez'], t[m]['ba'], t[m]['mcc']] for m in MODELS],
              [44, 12, 12, 12, 12, 12, 12, 12, 12, 12])
    r = text(ws, r, 'Univariate hit count against its label-permutation null (1,000 permutations)', HF)
    r = table(ws, r, ['Route', 'Observed hits', 'Null mean', 'Null SD', 'Tests m', 'Null hit rate q', 'p', 'Inflation nominal (5 %)', 'Inflation calibrated', 'Null 97.5th percentile'],
              [[f'Route {q}', *[t['T_' + q][k] for k in ['beobachtet', 'null_mittel', 'null_sd', 'tests', 'null_rate', 'p', 'inf_nominal', 'inf_kalibriert', 'perzentil_97_5']]] for q in ROUTES])
    r = text(ws, r, 'Paired stratified bootstrap of nested Q2 (20,000 draws)', HF)
    r = table(ws, r, ['Comparison', 'Difference', '2.5 %', '97.5 %', 'Interval entirely above 0'],
              [[f"Route {k.split('-')[0]} - {KNAME[k.split('-')[1]]}", v['diff'], v['lo'], v['hi'], 'yes' if v['ganz_ueber_0'] else 'no'] for k, v in t['bootstrap'].items()])
    for row in ws.iter_rows():
        for cc in row:
            if isinstance(cc.value, float): cc.number_format = '0.0000'

# ---------------------------------------------------------------- Robustness
ws = wb.create_sheet('Robustness')
r = text(ws, 1, 'Robustness of H1 to H4 against numerically unstable predictions (decision EW 24.10.2026, after the results under the original rule were known)', TF)
r = text(ws, r, 'Cause: autoscaling of descriptors that are almost constant among the training compounds but differ for the held-out compound.')
r = text(ws, r, f"R1: series with at least one PLS prediction > 2 or < -1 in any model. Excluded: {', '.join(R['ausgeschlossen_R1'])}.")
r = text(ws, r, f"R2 (label-free, stricter): series in which a held-out compound reaches |z| > 100 on at least one descriptor (descriptors constant in training not counted). Excluded: {', '.join(R['ausgeschlossen_R2'])}.")
r = text(ws, r, 'H4 limit for fewer series: smallest k with P(X >= k) < 0.05 under Binomial(n, 0.075), same derivation as in the specification.')
rows = []
for key, lab in [('alle16', 'all 16 series'), ('R1', 'without R1 series'), ('R2', 'without R2 series')]:
    X = R[key]
    rows.append([lab, X['n'], f"{X['H1']['A1']['positiv']}/{X['n']}", X['H1']['A1']['p'], X['H1']['A1u']['p'], X['H1']['A2']['p'], X['H2']['p'],
                 f"{X['H2']['ordnung_vollstaendig']}/{X['n']}", X['H3']['rho'], X['H3']['p'], f"{X['H4']['anzahl']} ({', '.join(X['H4']['datensaetze']) or '-'})",
                 X['H4']['grenze'], 'yes' if X['H4']['widerlegt'] else 'no', X['H4']['median_beste_minus_kontrolle'],
                 ', '.join(h for h, ok in X['holm_m4'].items() if ok) or '-'])
r = table(ws, r + 1, ['Series', 'n', 'H1 A1 positive', 'H1 p A1', 'H1 p A1u', 'H1 p A2', 'H2 p', 'H2 full order', 'H3 rho', 'H3 p', 'H4 count (series)',
                      'H4 limit', 'H4 refuted', 'H4 median best - control', 'Holm (m = 4) confirmed'], rows, [20, 5, 10, 10, 10, 10, 10, 10, 10, 10, 26, 8, 10, 12, 18])
r = text(ws, r, 'AUC instead of Q2, all 16 series', HF)
r = table(ws, r, ['Route', 'H1: apparent - nested AUC, positive', 'Median', 'Wilcoxon p (one-sided)'],
          [[f'Route {q}', f"{v['positiv']}/16", v['median'], v['p']] for q, v in R['AUC_H1'].items()])
a4 = R['AUC_H4']
r = text(ws, r, f"H4 with AUC (same bootstrap draws): {a4['anzahl']} of 16 ({', '.join(a4['datensaetze']) or '-'}), limit {a4['grenze']}, refuted: {'yes' if a4['widerlegt'] else 'no'}. Median best route - control {a4['median_beste_minus_kontrolle']:.4f}.")
r = table(ws, r + 1, ['Series and route', 'AUC difference (route - no filter)', '2.5 %', '97.5 %'], [[k, v['diff'], v['lo'], v['hi']] for k, v in a4['intervalle'].items()])

# ---------------------------------------------------------------- Exploratory
ws = wb.create_sheet('Exploratory')
r = text(ws, 1, 'Exploratory analyses named in the pre-specification (descriptive, no test claim)', TF)
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
grp(lambda d: zelle(d)[2], '1. Band of the minor class (A 25-35 % against B >35-50 %)')
grp(lambda d: zelle(d)[1], '2. Structural diversity')
grp(lambda d: zelle(d)[0], '3. Activity type')
r = text(ws, r, '4. Size effect: subsets of n = 25 and n = 60 drawn from the same parent pool (30 % minor class)', HF)
rows = []
for d in TEIL:
    t = T[d]; I = INFO[d]
    rows.append([d, ART[I['zelle'][1]], I['dokument'], I['n'], *[t['app_' + q]['q2'] - t[q]['q2'] for q in ROUTES], *[t[q]['q2'] for q in ROUTES],
                 *[t[q]['p_q2'] for q in ROUTES], *[t['T_' + q]['inf_kalibriert'] for q in ['A1u', 'A1', 'A2']], t['ohneFilter']['q2'], t['RF']['q2']])
r = table(ws, r, ['Subset', 'Parent', 'Parent document', 'n', 'Optimism A1', 'Optimism A1u', 'Optimism A2', 'Nested Q2 A1', 'Nested Q2 A1u', 'Nested Q2 A2',
                  'p(Q2) A1', 'p(Q2) A1u', 'p(Q2) A2', 'Calibrated inflation A1u', 'Calibrated inflation A1', 'Calibrated inflation A2', 'Q2 no filter', 'Q2 random forest'], rows)
r = text(ws, r, 'Nested Q2 values far below 0 are affected by numerically unstable predictions (see Robustness). For these, AUC is the more informative measure (sheets DS21 to DS24).')
r += 1
r = text(ws, r, '5. Does an excess of A1u hits over the label-free null coincide with a nominally significant nested Q2 (as for S. aureus and S. pneumoniae in the manuscript)?', HF)
rows = []; ct = Counter()
for d in PRIM:
    t = T[d]; exc = t['T_A1u']['beobachtet'] > t['T_A1u']['perzentil_97_5']; sig = t['A1u']['p_q2'] <= 0.05
    ct[(exc, sig)] += 1
    rows.append([d, t['T_A1u']['beobachtet'], t['T_A1u']['perzentil_97_5'], t['T_A1u']['p'], 'yes' if exc else 'no', t['A1u']['q2'], t['A1u']['p_q2'], 'yes' if sig else 'no'])
r = table(ws, r, ['Series', 'A1u observed hits', 'Null 97.5th percentile', 'p (hit count)', 'Excess (observed > 97.5th percentile)', 'Nested Q2 A1u', 'p(Q2) A1u', 'p(Q2) <= 0.05'], rows)
r = text(ws, r, 'Excess is defined here as observed hits above the 97.5th percentile of the label-permutation null. The specification names the question but not the criterion, so this definition was set when the analysis was written.')
r = table(ws, r, ['', 'Nested Q2 A1u p <= 0.05', 'p > 0.05'], [['Hit excess', ct[(True, True)], ct[(True, False)]], ['No hit excess', ct[(False, True)], ct[(False, False)]]])
r = text(ws, r, '6. Origin of the class boundary', HF)
r = text(ws, r, 'All 16 natural series use the fallback rule (stage 2), because none of the publications states a boundary. The factor therefore has no variation and cannot be analysed.')
r += 1
r = text(ws, r, '7. Median rule (active = better than the median of the data set) as sensitivity analysis', HF)
r = text(ws, r, 'Not yet computed. The full analysis with permutations requires a separate compute run (decision EW pending). See Deviations.')

# ---------------------------------------------------------------- Deviations
DEV = [
    ['22.10.2026', 'Run time for the choice of the permutation number was measured with permuted labels on all 20 data sets, not only on the one with the largest descriptor pool. The code was then vectorised and the run time measured again.', 'More information for the decision. No result with observed labels was computed.', 'None on results.'],
    ['23.10.2026', 'Nested models were run with 2,000 permutations (decision EW, based on run times measured with permuted labels).', 'Largest number of the rule within the available compute time.', 'Monte Carlo half-width of a p-value near 0.05 is 0.010.'],
    ['23.10.2026', 'Column "n" of the design table counted compounds before removing censored values that cannot be classified (DS01 51 -> 50, DS04 43 -> 40, DS13 36 -> 35, original rule).', 'Rule-conform, documentation only.', 'None.'],
    ['24.10.2026', 'The seed of the bootstrap over data sets (effect estimate of H4) was not specified. Used: default_rng([93, 0, 8]) under the original rule and default_rng([93, 0, 208]) in the main analysis.', 'Gap in the specification.', 'Affects only the interval of the H4 effect estimate.'],
    ['24.10.2026', 'Numerically unstable leave-one-out predictions were found in several series (Q2 far below 0). Robustness checks were added after the results were known: R1, R2 and AUC instead of Q2 (decision EW 04:59). The original analysis was not changed.', 'Instability by autoscaling of almost constant descriptors, not a code error. The pipeline follows the specification.', 'All four hypothesis decisions unchanged in every check (sheet Robustness).'],
    ['24.10.2026', 'Supplementary rule: reference drugs, clinical candidates and tool compounds of other classes measured for comparison or as controls, and values taken from other publications, are not part of the series. Frozen 24.10.2026 17:20 (SHA-256 36668775...), before any model of the new analysis. Decisions per compound from the article text (Excluded_compounds). Data set selection repeated mechanically: MIC/diverse/B Hansen 2019 -> Ogata 1987, Receptor/diverse/A Bonifazi 2023 -> Suckling 2007, MIC parent pool Carradori 2017 -> Spink 2015. Seeds analysis code + 200.', 'Gap in the selection rule, found while examining the unstable predictions, i.e. after the results under the original rule were known.', 'This analysis is the main analysis (decision EW 24.10.2026 14:23, before it was computed). All four hypothesis decisions are the same in both analyses (sheet Comparison).'],
    ['25.10.2026', '408 of the 2,900 permutation packages of the main analysis were computed in parallel by the automatic continuation job, in a separate container with the same code and seeds. They were merged after checking.', 'Technical.', 'None. All 2,432 packages computed by both are byte-identical, 12 randomly chosen permutations of the merged packages were recomputed and agree exactly.'],
    ['25.10.2026', 'For DS01, DS03, DS04 and DS15, which are identical in both analyses, the results were taken over as specified. Whether the bootstrap there uses seed 8 or 208 is not unambiguous in the supplementary specification. 208 was used.', 'Ambiguity of the specification.', 'None: the H4 classification of these four series is identical with both seeds.'],
    ['25.10.2026', 'The script for the label-free criterion R2 (maxz.py) was rewritten, because the original computation was not stored as a file.', 'Documentation.', 'None: it reproduces the values of 24.10.2026 up to rounding (relative 4e-13).'],
    ['open', 'Sensitivity analysis with the median rule not yet computed.', 'Requires a separate compute run with permutations, decision EW pending.', 'Reported descriptively only, does not enter H1 to H4.'],
    ['note', 'Origin of the class boundary as exploratory factor: all 16 series use the fallback rule, so the factor has no variation.', '-', 'Exploratory analysis not possible.'],
]
ws = sheet('Deviations', ['Date', 'Deviation', 'Reason', 'Effect'], DEV, [12, 100, 50, 60])
for row in ws.iter_rows(min_row=2):
    for cc in row: cc.alignment = Alignment(wrap_text=True, vertical='top')

# ---------------------------------------------------------------- main analysis only
if HAUPT:
    CR = json.load(open('kontrollen/compound_records.json'))
    names = {}
    for dok, recs in CR.items():
        for rec in recs:
            names[(dok, rec[0])] = (rec[1], rec[2])
    rows = []
    for d in ALL:
        I = INFO[d]
        for cid in I.get('entfernt') or []:
            k, n = names.get((I['dokument'], cid), ('', ''))
            rows.append([d, I['dokument'], I.get('doi'), cid, k, n])
    ws = sheet('Excluded_compounds', ['Data set', 'ChEMBL document', 'DOI', 'ChEMBL compound', 'Compound key in the article', 'Name (ChEMBL)'], rows, [10, 16, 28, 16, 24, 50])
    text(ws, len(rows) + 3, 'Each decision was taken from the article text and is documented with the source passage in the table of control and reference compounds per data set (Kontrollsubstanzen_je_Datensatz, revision archive).')
    A0 = json.load(open('ergebnisse/auswertung.json'))['hypothesen']; R0 = json.load(open('ergebnisse/robustheit.json'))
    def zeile(h, x):
        if h == 'H1': return f"p = {x['H1']['p']:.5f}, positive {x['H1']['positiv']}/16, median {x['H1']['median_diff']:.3f}"
        if h == 'H2': return f"p = {x['H2']['p']:.2e}, full order {x['H2']['sekundaer']['ordnung_vollstaendig']}/16"
        if h == 'H3': return f"rho = {x['H3']['rho']:.3f}, p = {x['H3']['p']:.3f}"
        return f"{x['H4']['anzahl']} of 16 ({', '.join(x['H4']['datensaetze'])}), median best - control {x['H4']['median_beste_minus_kontrolle']:.3f}"
    def ent(h, x):
        if h == 'H4': return 'refuted' if x['H4']['widerlegt'] else 'not refuted'
        return 'confirmed' if x['holm_m4'][h] else 'not confirmed'
    rows = [[h, zeile(h, A0), ent(h, A0), zeile(h, H), ent(h, H)] for h in ['H1', 'H2', 'H3', 'H4']]
    rows += [['Series excluded by R1', ', '.join(R0['ausgeschlossen_R1']), '', ', '.join(R['ausgeschlossen_R1']), ''],
             ['Series excluded by R2', ', '.join(R0['ausgeschlossen_R2']), '', ', '.join(R['ausgeschlossen_R2']), ''],
             ['H4 with AUC', f"{R0['AUC_H4']['anzahl']} of 16 ({', '.join(R0['AUC_H4']['datensaetze'])})", 'not refuted' if not R0['AUC_H4']['widerlegt'] else 'refuted',
              f"{R['AUC_H4']['anzahl']} of 16 ({', '.join(R['AUC_H4']['datensaetze'])})", 'not refuted' if not R['AUC_H4']['widerlegt'] else 'refuted']]
    ws = sheet('Comparison', ['Hypothesis', 'Original rule: result', 'Original rule: decision', 'Main analysis: result', 'Main analysis: decision'], rows, [22, 60, 18, 60, 18])


# ---------------------------------------------------------------- Unstable_predictions (decision EW 25.10.2026)
V_ = json.load(open(f'{ERG}/vorhersagen.json'))
ENT_ = json.load(open('kontrollen/entscheidungen.json'))
CR_ = json.load(open('kontrollen/compound_records.json'))
spec_ds = importlib.util.spec_from_file_location('dsmod', 'datasets.py'); dsm = importlib.util.module_from_spec(spec_ds); spec_ds.loader.exec_module(dsm)
PMOD = ['app_A1', 'app_A1u', 'app_A2', 'A1', 'A1u', 'A2', 'ohneAuswahl', 'ohneFilter']
GREY = PatternFill('solid', fgColor='D9D9D9'); ORA = PatternFill('solid', fgColor='F8CBAD'); YEL = PatternFill('solid', fgColor='FFF2CC')
ws = wb.create_sheet('Unstable_predictions')
r = text(ws, 1, 'Compounds of the series with numerically unstable predictions (criterion R1)', TF)
for t in ['One block per affected series, listing every compound that ChEMBL holds for the assay of the publication.',
          'Orange: at least one leave-one-out prediction of a PLS model above 2 or below -1 (R1). Unstable predictions in bold red.',
          'Yellow: no unstable prediction, but largest label-free |z| above 100 (criterion R2).',
          'Grey: reference or control compound excluded by the supplementary rule, not part of the analysis.' if HAUPT else
          'Column "Excluded in main analysis": reference or control compound that the supplementary rule removes in the main analysis. Under the original rule it is part of the series.',
          'Largest |z|: distance of the held-out compound from the mean of the other compounds of the series in units of their SD, over all cleaned descriptors, descriptors constant among the other compounds (SD < 1e-6) not counted. It does not depend on the class labels.',
          'Cause of the instability: autoscaling of a descriptor that is almost constant among the training compounds. The held-out compound then lies extremely far from them on this descriptor, and the linear PLS prediction overshoots.',
          'Generic scaffold: Bemis-Murcko scaffold with all atoms as carbon and all bonds single, as in the selection. "Compounds sharing scaffold" counts the compounds of the series with the same generic scaffold.']:
    r = text(ws, r, t)
r += 1
for d in R['ausgeschlossen_R1']:
    I = INFO[d]; yk = pd.read_csv(f'{DSD}/{d}_y.csv'); Xk = pd.read_csv(f'{DSD}/{d}_X.csv')
    idx = {c: i for i, c in enumerate(yk.chembl_id)}; Am = Xk.values.astype(float); cols_ = np.array(Xk.columns)
    mz = {}
    for i in range(len(Am)):
        tr = np.delete(Am, i, 0); mu = tr.mean(0); sd = tr.std(0); ok = sd >= 1e-6
        zz = np.abs((Am[i, ok] - mu[ok]) / sd[ok]); j = int(np.argmax(zz)); mz[yk.chembl_id[i]] = (float(zz[j]), cols_[ok][j], float(sd[ok][j]))
    aid = dsm.q("select assay_id from assays where chembl_id=?", (I['assay'],))[0][0]
    acts = dsm.q("""select molregno, standard_relation, standard_value, standard_units from activities where assay_id=? and src_id='1'
                    and data_validity_comment is null and potential_duplicate='0' and standard_value is not null and molregno is not null""", (str(aid),))
    per = {}
    for m, rel, val, unit in acts: per.setdefault(dsm.parent_of(int(m)), []).append((rel, float(val), unit))
    keys = {x[0]: (x[1], x[2]) for x in CR_.get(I['dokument'], [])}
    weg_main = set(ENT_.get(I['dokument'], {}).get('ausschluss', []))
    mols = {}
    for m, lst in per.items():
        mol, cid = dsm.mol_of(m); mols[cid] = (mol, lst)
    scaf = {cid: Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(mol))) for cid, (mol, _) in mols.items() if mol is not None and cid in idx}
    scc = Counter(scaf.values())
    ha_med = float(np.median([Chem.MolFromSmiles(x).GetNumHeavyAtoms() for x in yk.smiles]))
    c_ = BYASSAY.get(I['assay'], {})
    r = text(ws, r, f"{d} – {' / '.join(zelle(d))} – {I['dokument']} ({c_.get('journal')}, {c_.get('jahr')}), assay {I['assay']}, {I.get('ziel')}. In the analysis: {len(yk)} compounds, median heavy atoms {ha_med:g}.", HF)
    head = ['ChEMBL ID', 'Key in article', 'Name (ChEMBL)', 'SMILES', 'Value(s)', 'Unit', 'Class', 'Status'] + (['Excluded in main analysis'] if not HAUPT else []) + \
           ['Heavy atoms', 'Compounds sharing scaffold', 'Largest |z|', 'Descriptor with largest |z|', 'SD of this descriptor among the others'] + [MNAME[m] for m in PMOD] + ['Reason']
    for j, h in enumerate(head, 1):
        x = ws.cell(r, j, h); x.font = HF; x.fill = FILL; x.alignment = Alignment(wrap_text=True, vertical='top')
    r += 1
    for cid in sorted(mols, key=lambda c: (c not in idx, idx.get(c, 0))):
        mol, lst = mols[cid]; drin = cid in idx; key, name = keys.get(cid, ('', ''))
        excl_here = HAUPT and cid in (I.get('entfernt') or [])
        status = 'in the analysis' if drin else ('excluded: reference or control compound (supplementary rule)' if excl_here else
                  ('not used: structure not usable' if mol is None else 'not used: censored value cannot be classified'))
        preds = [V_[d]['pred'][m][idx[cid]] if drin else None for m in PMOD]
        unst = drin and any(p > 2 or p < -1 for p in preds)
        z, dn, dsd = mz.get(cid, (None, '', None))
        nsc = scc.get(scaf.get(cid)) if drin else None
        reason = ''
        if unst:
            bad = [MNAME[m] for m, p in zip(PMOD, preds) if p > 2 or p < -1]
            hz = mol.GetNumHeavyAtoms()
            reason = (f"{'Only compound with its generic scaffold in the series' if nsc == 1 else f'Shares its generic scaffold with {nsc - 1} other compounds of the series'}"
                      f", {hz} heavy atoms (series median {ha_med:g}). Descriptor {dn} has SD {dsd:.2g} among the other compounds, the held-out compound lies at |z| = {z:,.0f}. "
                      f"Unstable in {len(bad)} model(s): {' / '.join(bad)}.")
        elif drin and z is not None and z > 100:
            reason = f'No unstable prediction, but |z| = {z:,.0f} above 100 (R2).'
        row = [cid, key, name, Chem.MolToSmiles(mol) if mol is not None else '', ', '.join(f'{a}{b:g}' for a, b, _ in lst), lst[0][2],
               ('active' if yk.y[idx[cid]] == 1 else 'inactive') if drin else '', status] + \
              (['yes' if cid in weg_main else ''] if not HAUPT else []) + \
              [mol.GetNumHeavyAtoms() if mol is not None else '', nsc or '', round(z, 1) if z is not None else '', dn, dsd if dsd is not None else ''] + \
              [round(p, 3) if p is not None else '' for p in preds] + [reason]
        p0 = len(head) - len(PMOD)
        for j, v in enumerate(row, 1):
            x = ws.cell(r, j, v); x.font = NF; x.alignment = Alignment(wrap_text=(j in (4, len(row))), vertical='top')
            # Change of 16.11.2026 (item G13): openpyxl stores a text beginning with "=" (column Value(s), e.g. "=7000")
            # as a formula. It is text and is now stored as such. The workbooks in workbooks/ were written before this
            # change and hold these cells as formulas (see revision2/tables/aenderungen_cg2_workbooks.py).
            if isinstance(v, str) and v.startswith('='): x.data_type = 's'
            if excl_here: x.fill = GREY
            elif unst: x.fill = ORA
            elif drin and z is not None and z > 100: x.fill = YEL
        if unst:
            for k, p in enumerate(preds):
                if p > 2 or p < -1: ws.cell(r, p0 + k).font = Font(name=F, bold=True, color='C00000')
        r += 1
    r += 2
for j, w in enumerate([15, 12, 30, 40, 16, 9, 9, 26] + ([10] if not HAUPT else []) + [8, 9, 10, 18, 12] + [11] * 8 + [70], 1):
    ws.column_dimensions[get_column_letter(j)].width = w

order = ['README', 'Selection', 'Summary', 'All_results_long'] + ALL + ['Robustness', 'Unstable_predictions', 'Exploratory', 'Deviations'] + (['Excluded_compounds', 'Comparison'] if HAUPT else []) + ['Selection_log', 'Candidates']
wb._sheets = [wb[n] for n in order]
wb.save(OUT)
print('saved', OUT, len(wb.sheetnames), 'sheets, long', len(long))
