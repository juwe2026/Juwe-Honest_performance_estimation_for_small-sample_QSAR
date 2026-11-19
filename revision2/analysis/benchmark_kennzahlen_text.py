"""Key figures for the results section 'Generalisation across 20 ChEMBL data sets'
(Revision 2). Reads only the sheet Summary of the results workbook of the main analysis.
Pre-specified quantities are already there. Computed here AFTERWARDS (post hoc, descriptive):
medians of the difference of each single route from the control without
filter, comparison of the control without filter with PLS-DA without selection, ranges."""
import sys, json, numpy as np, openpyxl
wb = openpyxl.load_workbook(sys.argv[1], read_only=True, data_only=True)
rows = list(wb['Summary'].iter_rows(values_only=True))
h = [i for i, r in enumerate(rows) if r and r[0] == 'Data set'][0]
H = list(rows[h]); D = [r for r in rows[h+1:] if r and r[0] and str(r[0]).startswith('DS')]
nat = [r for r in D if int(r[0][2:]) <= 16]; sub = [r for r in D if int(r[0][2:]) >= 21]
c = lambda n, R=nat: np.array([r[H.index(n)] for r in R], float)
out = {'n_series': len(nat)}
for n in ['n', 'p (cleaned)', 'G', 'G/p', 'Calibrated inflation A1u', 'Calibrated inflation A1', 'Calibrated inflation A2',
          'Optimism A1 (apparent - nested Q2)', 'Optimism A1u', 'Optimism A2', 'Q2 no selection', 'Q2 no filter', 'Q2 random forest']:
    v = c(n); out[n] = dict(min=round(v.min(), 4), median=round(float(np.median(v)), 4), max=round(v.max(), 4))
for r in ['A1', 'A1u', 'A2']:
    out[f'optimism {r} positive'] = int((c('Optimism A1 (apparent - nested Q2)' if r == 'A1' else f'Optimism {r}') > 0).sum())
    out[f'p(Q2) {r} <= 0.05'] = int((c(f'p(Q2) {r}') <= 0.05).sum())
nf = c('Q2 no filter'); ns = c('Q2 no selection')
for r in ['A1', 'A1u', 'A2']:
    d = c(f'Nested Q2 {r}') - nf
    out[f'POSTHOC median {r} - no filter'] = round(float(np.median(d)), 4)
    out[f'POSTHOC {r} > no filter (point estimate)'] = int((d > 0).sum())
best = np.max([c('Nested Q2 A1'), c('Nested Q2 A1u'), c('Nested Q2 A2')], axis=0) - nf
out['median best route - no filter (pre-specified effect, check)'] = round(float(np.median(best)), 4)
out['POSTHOC no filter > no selection'] = int((nf > ns).sum())
out['POSTHOC median no filter - no selection'] = round(float(np.median(nf - ns)), 4)
out['subsets'] = {r[0]: dict(n=r[H.index('n')], optimism_A1=round(r[H.index('Optimism A1 (apparent - nested Q2)')], 3),
                  calibrated_A2=round(r[H.index('Calibrated inflation A2')], 2)) for r in sub}
json.dump(out, open(sys.argv[2], 'w'), indent=1, ensure_ascii=False); print(json.dumps(out, indent=1, ensure_ascii=False))
