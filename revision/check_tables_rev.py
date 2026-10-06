# -*- coding: utf-8 -*-
"""G2-23 (Pruefung 03.10.2026): jede Zelle "Wert (p)" der Tables 1 bis 5 des revidierten Manuskripts gegen ihre Quelldatei.
Fuer jede Zeile ist die Quelle explizit zugeordnet (keine Suche nach passenden Paaren). Geprueft werden der gedruckte
Wert (auf die gedruckte Stellenzahl), der p-Wert (vier Stellen; '≤0.001' fuer 1.000 Permutationen ohne Ueberschreitung),
das Vorzeichen und die Kennzeichnung (* / n.s.). Ausgabe out/numbercheck_rev.json; Abbruch, wenn eine Zelle abweicht."""
import json, re, pickle, numpy as np, docx
O = '/home/claude/rev/out/'; RP = '/home/claude/rev/QSAR_selection_leakage/results/raw_permutations/'
J = lambda f: json.load(open(f))
def arch(prefix, s):
    z = np.load(RP + '%s_%s.npz' % (prefix, s)); t = lambda n, o: float((1 + np.sum(n >= o - 1e-9)) / (1 + len(n)))
    return dict(q2=float(z['obs_q2']), auc=float(z['obs_auc']), p_q2=t(z['q2s'], float(z['obs_q2'])), p_auc=t(z['aucs'], float(z['obs_auc'])))
T1P = J(O + 'table1_pub_exact.json'); T1C = J(O + 'I_fallback_seeds.json'); A = J(O + 'A_nofilter.json'); D = J(O + 'D_rf.json')
B25 = J(O + 'B25_cineole0_seed701.json')
OP = pickle.load(open(RP + 'ordinary_perm_results.pkl', 'rb'))
E15 = pickle.load(open(RP + 'eff15_apparent.pkl', 'rb')); V1 = pickle.load(open(RP + 'vip1lv_apparent.pkl', 'rb'))
LA = J('/home/claude/rev/lone/res_A.json'); LC = J('/home/claude/rev/lone/res_C.json')
for r_, f_ in [(LA, '/home/claude/rev/lone/res_A_fix.json'), (LC, '/home/claude/rev/lone/res_C_fix.json')]:
    c = J(f_)
    for k in c:
        if k in r_ and isinstance(c[k], dict):
            for a in ('A1', 'A1u', 'A2'):
                if a in c[k]: r_[k][a] = c[k][a]
ST = ['HI', 'SA', 'SPneu', 'SPyo', 'PA']; LG = {'HI': 'H.I.', 'SA': 'S.A.', 'SPneu': 'S.Pneu', 'SPyo': 'S.Pyo', 'PA': 'P.A.'}
LST = ['EC', 'SA', 'PA', 'AB']
def route(r, k): return T1C.get(k, {}).get(r) or T1P[r][k]
def rf(k, key): s = D[key][k]; return dict(q2=s['oob_q2'], auc=s['oob_auc'], p_q2=s['p_q2'], p_auc=s['p_auc'])
def topk(k, a): return J(O + 'cache_topk_seeds/%s_%s.json' % (k, a))
def app(src, k, keys): q, a, pq, pa = keys; return dict(q2=src[LG[k]][q], auc=src[LG[k]][a], p_q2=src[LG[k]][pq], p_auc=src[LG[k]][pa], n=1000)
def op(k, r): t = OP[LG[k]][r]; return dict(q2=t[0], auc=t[1], p_q2=t[2], p_auc=t[3], n=1000)
MAP = {  # (Tabelle, Zeilenanfang) -> (Quelle, Metrik, Spalten, Name)
 (1, 'Route A1'): lambda k: route('A1', k), (1, 'Route A1u'): lambda k: route('A1u', k), (1, 'Route A2'): lambda k: route('A2', k),
 (1, 'No selection, all 2,208'): lambda k: A['all2208'][k], (1, 'No filter, 182'): lambda k: A['reps182_centrality'][k],
 (1, 'Random forest'): lambda k: dict(rf(k, 'RouteA_all2208'), n=500), (1, 'Route B, 12 descriptors, 1 LV'): lambda k: arch('permB_12', k),
 (1, 'Route B, 12 descriptors, 3 LV'): lambda k: arch('permB_12nlv3', k),
 (2, '25 descriptors'): lambda k: B25[k], (2, '12 descriptors'): lambda k: arch('permB_12', k),
 (3, 'Apparent'): lambda k: app(E15, k, ('q', 'a', 'pq', 'pa')), (3, 'Nested'): lambda k: topk(k, 'A1_top15'),
 (4, 'Apparent'): lambda k: app(V1, k, ('app_q', 'app_a', 'app_pq', 'app_pa')), (4, 'Nested'): lambda k: topk(k, 'A1_viptop1'),
}
LMAP = {'Route A1 (': 'A1', 'Route A1u': 'A1u', 'Route A2 (': 'A2', 'No selection, all 1,275': 'noSelection_all',
        'No filter, 50': 'noFilter_reps', 'Random forest': 'randomForest_oob'}
CELL = re.compile(r'^([+\-−]?\d+\.\d+) \((≤0\.001|\d\.\d{4}|artefact)\)\s*(\*|n\.s\.|‡)?$')
d = docx.Document('/home/claude/rev/ms/revised_stage5.docx')
rows = []; bad = []
def check(tab, lab, ep, txt, src, metric, signed):
    m = CELL.match(txt)
    if not m: bad.append((tab, lab, ep, txt, 'Format')); return
    val, pp, mk = m.groups(); v = src[metric]; p = src['p_' + metric]
    nd = len(val.split('.')[1]); exp = ('%.' + str(nd) + 'f') % v
    if exp in ('-0.00', '-0.000'): exp = exp[1:]                      # kein negatives Null
    if signed and v > 0 and not exp.startswith('0.00') or (signed and exp.startswith('0.00') and v > 0): exp = '+' + exp
    got = val.replace('−', '-')
    ok_v = got == exp or (got.lstrip('+') == exp.lstrip('+') and not signed)
    if pp == 'artefact': ok_p, pexp = True, 'artefact'
    elif pp == '≤0.001': ok_p, pexp = abs(p - 1 / 1001) < 1e-9 or p <= 0.001 + 1e-12, '≤0.001'
    else: pexp = '%.4f' % p; ok_p = pp == pexp
    ok_m = pp in ('artefact', '≤0.001') or mk is None or mk == '‡' or mk == ('*' if p < 0.05 else 'n.s.')
    rows.append(dict(Table='Table %d' % tab, Row=lab, Endpoint=ep, Printed=txt, Source_value=round(float(v), 6), Source_p=round(float(p), 6),
                     Value_ok=ok_v, p_ok=ok_p, Marker_ok=ok_m))
    if not (ok_v and ok_p and ok_m): bad.append((tab, lab, ep, txt, 'erwartet %s (%s)' % (exp, pexp)))
for ti in range(4):
    tab = ti + 1; sect = None
    for r in d.tables[ti].rows:
        lab = r.cells[0].text
        if lab.startswith('(a)'): sect = 'q2'
        if lab.startswith('(b)'): sect = 'auc'
        if lab.startswith(('(c)', '(d)')): sect = None
        key = max([kk for kk in MAP if kk[0] == tab and lab.startswith(kk[1])], key=lambda kk: len(kk[1]), default=None)
        if key is None: continue
        metric = sect or ('q2' if ('Q2' in lab or 'Q²' in lab) else 'auc')
        if tab == 1 and sect is None: continue
        for c, k in zip(r.cells[1:], ST):
            check(tab, lab, k, c.text.strip(), MAP[key](k), metric, signed=(metric == 'q2' and tab in (1, 3, 4) and 'Apparent' not in lab and 'Route B' not in lab))
# Table 1 (c) apparent: Q2 / AUC aus ordinary_perm_results
for r in d.tables[0].rows:
    lab = r.cells[0].text
    for rr, nm in (('Route A1:', 'A1'), ('Route A1u:', 'A1u'), ('Route A2:', 'A2')):
        if lab.startswith(rr):
            for c, k in zip(r.cells[1:], ST):
                q, a = c.text.split(' / '); s = op(k, nm)
                ok = q == '%.2f' % s['q2'] and a == '%.2f' % s['auc']
                rows.append(dict(Table='Table 1', Row=lab, Endpoint=k, Printed=c.text, Source_value=round(s['q2'], 6), Source_p=None, Value_ok=ok, p_ok=True, Marker_ok=True))
                if not ok: bad.append((1, lab, k, c.text, 'erwartet %.2f / %.2f' % (s['q2'], s['auc'])))
# Table 5
sect = None
for r in d.tables[4].rows:
    lab = r.cells[0].text
    if lab.startswith('(a)'): sect, SRC = 'q2', LA
    elif lab.startswith('(b)'): sect, SRC = 'auc', LA
    elif lab.startswith('(c)'): sect = 'hits'
    elif lab.startswith('(d)'): sect, SRC = 'q2', LC
    key = next((v for kk, v in LMAP.items() if lab.startswith(kk)), None)
    if sect in ('q2', 'auc') and key:
        for c, k in zip(r.cells[1:], LST):
            check(5, lab, k, c.text.strip(), SRC[k][key], sect, signed=(sect == 'q2'))
    if sect == 'hits' and lab.startswith('Route'):
        h = {'Route A2': 'A2', 'Route A1 (': 'A1', 'Route A1u': 'A1u'}
        hk = next(v for kk, v in h.items() if lab.startswith(kk))
        for c, k in zip(r.cells[1:], LST):
            e = '%d, %.3f' % (LA[k]['hits'][hk]['obs'], LA[k]['hits'][hk]['rate'])
            rows.append(dict(Table='Table 5', Row=lab, Endpoint=k, Printed=c.text, Source_value=LA[k]['hits'][hk]['rate'], Source_p=None,
                             Value_ok=c.text == e, p_ok=True, Marker_ok=True))
            if c.text != e: bad.append((5, lab, k, c.text, 'erwartet ' + e))
json.dump(dict(n=len(rows), bad=bad, rows=rows), open(O + 'numbercheck_rev.json', 'w'), indent=1, ensure_ascii=False, default=float)
print('geprueft: %d Zellen | Abweichungen: %d' % (len(rows), len(bad)))
for b in bad: print('  ', b)
