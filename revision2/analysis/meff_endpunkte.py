"""Effective number of independent endpoints (Li & Ji 2005, Nyholt 2004) from the
binary activity vectors. Independent recomputation of Table S19.
Inputs: Additional file 2 (Activity_vectors_31) and Additional file 4
(Peptides_MIC_classes) of the revision of 06.10.2026."""
import sys, numpy as np, openpyxl, json
AF = sys.argv[1]
def lade_mono():
    w = openpyxl.load_workbook(f"{AF}/Additional_file_2_machine_readable_tables 06.10.2026.xlsx", read_only=True, data_only=True)
    r = list(w['Activity_vectors_31'].iter_rows(values_only=True))
    return list(r[0][1:]), np.array([x[1:] for x in r[1:]], float)
def lade_pep(regel):
    w = openpyxl.load_workbook(f"{AF}/Additional_file_4_second_MIC_data_set 06.10.2026.xlsx", read_only=True, data_only=True)
    r = list(w['Peptides_MIC_classes'].iter_rows(values_only=True)); h = r[0]
    idx = [h.index(f'Class rule {regel} {s}') for s in ['EC','SA','PA','AB']]
    X = np.array([[1.0 if x[i] == 'active' else 0.0 for i in idx] for x in r[1:]])
    assert all(x[i] in ('active','non-active') for x in r[1:] for i in idx)
    return ['E. coli','S. aureus','P. aeruginosa','A. baumannii'], X
def meff(X):
    C = np.corrcoef(X, rowvar=False)          # Pearson = phi for binary vectors
    lam = np.sort(np.linalg.eigvalsh(C))[::-1]
    M = len(lam)
    lj = float(np.sum((lam >= 1).astype(float) + (lam - np.floor(lam))))
    ny = float(1 + (M - 1) * (1 - np.var(lam, ddof=1) / M))
    off = C[np.triu_indices(M, 1)]
    return C, lam, lj, ny, off.min(), off.max()
out = {}
nm, X = lade_mono()
fall = [('Monoterpenoide, 5 Endpunkte', nm, X),
        ('Monoterpenoide, 4 auswertbare (ohne P. aeruginosa)', nm[:4], X[:, :4])]
for rg in ['A', 'C']:
    n2, X2 = lade_pep(rg); fall.append((f'Peptide, Regel {rg}, 4 Endpunkte', n2, X2))
for name, n, X in fall:
    C, lam, lj, ny, lo, hi = meff(X)
    print(f'\n== {name}  (n = {X.shape[0]}, actives per endpoint {X.sum(0).astype(int).tolist()})')
    print('   phi pairwise', round(lo, 2), 'to', round(hi, 2))
    print('   eigenvalues', np.round(lam, 4).tolist())
    print(f'   Meff Li & Ji {lj:.2f}   Nyholt {ny:.2f}')
    out[name] = dict(n=int(X.shape[0]), aktive=X.sum(0).astype(int).tolist(), endpunkte=n,
                     phi=np.round(C, 4).tolist(), eigenwerte=np.round(lam, 4).tolist(),
                     meff_li_ji=round(lj, 2), meff_nyholt=round(ny, 2), phi_min=round(lo, 4), phi_max=round(hi, 4))
json.dump(out, open(sys.argv[2], 'w'), ensure_ascii=False, indent=1)
