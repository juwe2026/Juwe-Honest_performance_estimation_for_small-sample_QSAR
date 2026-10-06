# Permutations-p-Wert mit Gleichstandsbehandlung (Revision 03.10.2026).
# p = (1 + #{Nullwerte >= beobachteter Wert}) / (1 + Zahl der Permutationen).
# AUC-Werte sind Brueche k/(n1*n0); gleiche Brueche werden im Rechner je nach
# Rechenweg mit Abweichungen um 1e-16 gespeichert. Der bisherige strikte Vergleich
# zaehlte solche Gleichstaende teilweise nicht mit und machte p etwas zu klein.
# Jetzt gilt ein Nullwert als >= beobachtet, wenn er hoechstens 1e-9 darunter liegt.
# Jeder Aufruf legt die Nullverteilung zur Pruefung unter out/nulls/ ab.
import os as _os, sys as _sys, numpy as _np
_TIE_EPS = 1e-9
_NULLDIR = '/home/claude/rev/out/nulls'
_CNT = [0]
def pval(null, obs, tag=None):
    null = _np.asarray(null, float); obs = float(obs)
    strict = float((1 + _np.sum(null >= obs)) / (1 + len(null)))
    tie = float((1 + _np.sum(null >= obs - _TIE_EPS)) / (1 + len(null)))
    try:
        _os.makedirs(_NULLDIR, exist_ok=True); _CNT[0] += 1
        nm = '%s_%s_%d_%04d' % (_os.path.basename(_sys.argv[0]).replace('.py', ''),
                               '-'.join(_sys.argv[1:]) or 'x', _os.getpid(), _CNT[0])
        _np.savez_compressed(_os.path.join(_NULLDIR, nm + '.npz'), null=null, obs=obs, strict=strict, tie=tie,
                             tag=str(tag), argv=' '.join(_sys.argv))
    except Exception:
        pass
    return tie
