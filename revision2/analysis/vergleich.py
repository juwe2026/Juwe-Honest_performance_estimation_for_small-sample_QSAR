"""Comparison: factors and optimism of the data sets of the study, code of Revision 1 against code of the benchmark."""
import sys, json, openpyxl, numpy as np
G = json.load(open(sys.argv[1])); S = json.load(open(sys.argv[2]))
hc = list(openpyxl.load_workbook(sys.argv[3], read_only=True, data_only=True)['HitCounts_all_data_sets'].iter_rows(values_only=True)); hh = list(hc[0])
key = {"Route A2": "A2", "Route A1": "A1", "Route A1u": "A1u"}
EPm = {'H. influenzae': 'HI', 'S. aureus': 'SA', 'S. pneumoniae': 'SPneu', 'S. pyogenes': 'SPyo'}
EPp = {'E. coli ATCC 25922': 'EC', 'S. aureus ATCC 29213': 'SA', 'P. aeruginosa ATCC 27853': 'PA', 'A. baumannii ATCC 19606': 'AB'}
rev1 = {}
for r in hc[1:]:
    if not r or r[0] is None: continue
    d = dict(zip(hh, r))
    if d['Data_set'].startswith('Monoterpenoids') and d['Endpoint'] in EPm: k = ('Monoterpenoide', EPm[d['Endpoint']])
    elif d['Data_set'].startswith('Octapeptides, main'): k = ('Octapeptide', EPp[d['Endpoint']])
    elif d['Data_set'].startswith('Gaucher'): k = ('Gaucher', 'Gaucher')
    else: continue
    rev1.setdefault(k, {})[key[d['Route']]] = d['Inflation_vs_own_rate']
gq = {('Gaucher','Gaucher'): {'A1': (0.709, 0.387), 'A1u': (0.759, 0.445), 'A2': (0.737, 0.377)}}
zeilen = []; out = {}
for ds in ['Monoterpenoide', 'Octapeptide', 'Gaucher']:
    print(f"\n{ds}: G Benchmark-Code {G[ds]['G']}, p {G[ds]['p']}")
    for ep, e in G[ds]['endpunkte'].items():
        for rt in ['A1', 'A1u', 'A2']:
            f_bm = e['faktoren'][rt]['kalibriert']; f_r1 = rev1[(ds, ep)][rt]
            if ds == 'Gaucher': o_r1 = gq[(ds, ep)][rt][0] - gq[(ds, ep)][rt][1]; n_r1 = gq[(ds, ep)][rt][1]
            else:
                s = S['Monoterpenoide' if ds == 'Monoterpenoide' else 'Octapeptide'][ep][rt]
                o_r1 = s['scheinbar'] - s['genestet']; n_r1 = s['genestet']
            o_bm = e['q2'][rt]['optimismus']; n_bm = e['q2'][rt]['genestet']
            zeilen.append((ds, ep, rt, f_r1, f_bm, o_r1, o_bm, n_r1, n_bm))
            print(f"  {ep:7s} {rt:4s} factor R1 {f_r1:6.2f} BM {f_bm:6.2f} ({100*(f_bm/f_r1-1):+5.1f} %) | optimism R1 {o_r1:+.3f} BM {o_bm:+.3f} | nested R1 {n_r1:+.3f} BM {n_bm:+.3f}")
z = np.array([[a[3], a[4], a[5], a[6], a[7], a[8]] for a in zeilen], float)
for rt in ['A1', 'A1u', 'A2']:
    idx = [i for i, a in enumerate(zeilen) if a[2] == rt]
    rel = z[idx, 1] / z[idx, 0] - 1
    print(f"{rt}: factor, relative deviation, median {100*np.median(rel):+.1f} %, range {100*rel.min():+.1f} to {100*rel.max():+.1f} %; "
          f"optimism, difference, median {np.median(z[idx,3]-z[idx,2]):+.3f}, max |d| {np.max(np.abs(z[idx,3]-z[idx,2])):.3f}")
json.dump([dict(datensatz=a[0], endpunkt=a[1], route=a[2], faktor_rev1=a[3], faktor_benchmarkcode=a[4], optimismus_rev1=round(a[5],4),
                optimismus_benchmarkcode=a[6], genestet_rev1=a[7], genestet_benchmarkcode=a[8]) for a in zeilen], open(sys.argv[4], 'w'), indent=1)
