"""Sicherheitsnetz: meldet Abbildungen, die praktisch leer sind.

Gemessen wird der Anteil nicht weisser Pixel, also alles Gezeichnete einschliesslich
grauer Histogramme, und zusaetzlich der Anteil farbiger Pixel. Faellt einer der beiden
Werte unter die Schwelle, greift vermutlich ein Datenfilter ins Leere und es werden nur
noch Achsen gezeichnet.
"""
import sys, os, glob
import numpy as np
from PIL import Image

INK_MIN = 0.012      # Anteil nicht weisser Pixel
COLOUR_MIN = 0.0015  # Anteil farbiger Pixel


def fractions(path):
    im = Image.open(path).convert("RGB")
    im.thumbnail((900, 900))
    a = np.asarray(im, dtype=int)
    mx = a.max(axis=2)
    mn = a.min(axis=2)
    ink = mx < 245                      # alles Gezeichnete
    colour = (mx - mn) > 25             # farbig, nicht grau
    return float(ink.mean()), float(colour.mean())


def main(folder):
    bad = []
    for f in sorted(glob.glob(os.path.join(folder, "*.png"))):
        ink, colour = fractions(f)
        ok = ink >= INK_MIN and colour >= COLOUR_MIN
        print("%-46s Zeichnung %5.2f %%  Farbe %5.2f %%%s"
              % (os.path.basename(f), 100 * ink, 100 * colour,
                 "" if ok else "   <-- VERDAECHTIG LEER"))
        if not ok:
            bad.append(f)
    if bad:
        print("\n%d Abbildung(en) pruefen!" % len(bad))
        return 1
    print("\nAlle Abbildungen enthalten Daten.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
