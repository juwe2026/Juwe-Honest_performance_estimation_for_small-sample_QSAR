"""Loading helpers: every figure reads from the archived result objects."""
import json, os, pickle
import numpy as np

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Corrected 2026-08-27: these defaults previously pointed at an absolute, machine-specific
# upload path (/mnt/user-data/uploads/...) that is not part of this archive and only happened
# to resolve on the machine the figures were last regenerated on. They now default to the
# archive's own bundled locations; QSAR_RAW / QSAR_XLS still override them if the raw
# permutation outputs or source spreadsheets live elsewhere.
RAW = os.environ.get("QSAR_RAW", os.path.join(_REPO_ROOT, "results", "raw_permutations"))
XLS = os.environ.get("QSAR_XLS", os.path.join(_REPO_ROOT, "data"))

SHORT = {"H.I.": "HI", "S.A.": "SA", "S.Pneu": "SPneu", "S.Pyo": "SPyo", "P.A.": "PA"}

def npy(name):
    return np.load(os.path.join(RAW, name), allow_pickle=True).item()

def pkl(name):
    with open(os.path.join(RAW, name), "rb") as fh:
        return pickle.load(fh)

def jsn(name):
    with open(os.path.join(RAW, name)) as fh:
        return json.load(fh)

def perm(prefix, strain):
    """Null distributions of one model. prefix maps to the route:
       perm_alt = Route A1 full, perm = Route A2 full, perm_uns = Route A1u,
       perm_a1e = Route A1 top-15 by effect size, perm15 = Route A2 top-15,
       perm_a1vip1 = Route A1 top-15 by VIP, permB_13 / permB_25 = Route B."""
    z = np.load(os.path.join(RAW, "%s_%s.npz" % (prefix, SHORT[strain])))
    return dict(q2s=z["q2s"], aucs=z["aucs"], obs_q2=float(z["obs_q2"]), obs_auc=float(z["obs_auc"]))

def xlsx(name):
    import openpyxl
    return openpyxl.load_workbook(os.path.join(XLS, name), data_only=True)

def sheet_rows(wb, title):
    ws = wb[title]
    rows = list(ws.iter_rows(values_only=True))
    hdr = list(rows[0])
    return hdr, [r for r in rows[1:] if r[0] is not None]


RECALC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_recalc")

def routeB12():
    """Route B after removal of the perfectly collinear C_Count: 12 descriptors,
    2,000 permutations per model. Replaces routeB_final.json / routeB_nlv3.pkl."""
    with open(os.path.join(RECALC, "routeB12.json")) as fh:
        return json.load(fh)


def vip_routeB():
    """Route B VIP values (two latent variables), apparent and mean over the LOO
    folds, for the 13-descriptor control and the 12-descriptor modelling set."""
    with open(os.path.join(RECALC, "vip_routeB.json")) as fh:
        return json.load(fh)
