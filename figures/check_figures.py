"""Safety net: reports figures that are practically empty.

Measured is the share of non-white pixels, that is everything drawn including
grey histograms, and in addition the share of coloured pixels. If one of the two
values falls below the threshold, a data filter presumably matches nothing and only
the axes are drawn.
"""
import sys, os, glob
import numpy as np
from PIL import Image

INK_MIN = 0.012      # share of non-white pixels
COLOUR_MIN = 0.0015  # share of coloured pixels


def fractions(path):
    im = Image.open(path).convert("RGB")
    im.thumbnail((900, 900))
    a = np.asarray(im, dtype=int)
    mx = a.max(axis=2)
    mn = a.min(axis=2)
    ink = mx < 245                      # everything drawn
    colour = (mx - mn) > 25             # coloured, not grey
    return float(ink.mean()), float(colour.mean())


def main(folder):
    bad = []
    for f in sorted(glob.glob(os.path.join(folder, "*.png"))):
        ink, colour = fractions(f)
        ok = ink >= INK_MIN and colour >= COLOUR_MIN
        print("%-46s drawing %5.2f %%  colour %5.2f %%%s"
              % (os.path.basename(f), 100 * ink, 100 * colour,
                 "" if ok else "   <-- SUSPICIOUSLY EMPTY"))
        if not ok:
            bad.append(f)
    if bad:
        print("\n%d figure(s) to check!" % len(bad))
        return 1
    print("\nAll figures contain data.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
