"""Route B, both descriptor sets (13 and 12 representatives), one script: permutation tests,
VIP tables and confusion matrices, written to results/recalc/results.json by default; override
with the QSAR_RECALC_DIR environment variable to write elsewhere.
"""
import sys, os, json, time

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
import core

OUT_DIR = os.environ.get("QSAR_RECALC_DIR", os.path.join(_REPO_ROOT, "results", "recalc"))
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PATH = os.path.join(OUT_DIR, "results.json")

t0 = time.time(); out = {}
SETS = {"D13": core.RED13, "D12": core.RED12}

# ---- Permutationstests, identische Seeds fuer beide Saetze -> Differenz ist rein der Deskriptor
out["perm"] = {}
for tag, cols in SETS.items():
    out["perm"][tag] = {}
    for nlv in (1, 3):
        out["perm"][tag][str(nlv)] = {}
        for s in core.STRAINS:
            r = core.perm_test(cols, s, nlv)
            out["perm"][tag][str(nlv)][s] = r
            print("%s nLV=%d %-7s Q2 %+.3f (p=%.4f)  AUC %.3f (p=%.4f)   [%.0f s]"
                  % (tag, nlv, s, r["q2"], r["p_q2"], r["auc"], r["p_auc"], time.time()-t0), flush=True)
        json.dump(out, open(OUT_PATH, "w"), indent=1)

# ---- VIP-Tabellen mit dem publizierten nLV je Stamm
out["vip"] = {tag: core.vip_table(cols) for tag, cols in SETS.items()}
# ---- VIP zusaetzlich bei nLV = 1 und 3, wie in Fig. S2 und S3 verwendet
out["vip_nlv"] = {tag: {str(n): core.vip_table(cols, {s: n for s in core.STRAINS})
                        for n in (1, 3)} for tag, cols in SETS.items()}
# ---- Konfusionsmatrizen
out["conf"] = {tag: {str(n): {s: core.confusion(cols, s, n) for s in core.STRAINS}
                     for n in (1, 3)} for tag, cols in SETS.items()}
out["descriptors"] = {tag: cols for tag, cols in SETS.items()}
json.dump(out, open(OUT_PATH, "w"), indent=1)
print("FERTIG nach %.0f s -> %s" % (time.time()-t0, OUT_PATH))
