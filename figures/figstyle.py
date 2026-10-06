"""Shared style for all figures of the manuscript.

Colour convention, used identically in every figure:

    Route A1   orange   #E69F00    supervised representative choice
    Route A1u  green    #009E73    label-free representative choice
    Route A2   blue     #0072B2    filter first, grouping second
    Route B    purple   #CC79A7    a priori literature set

    apparent estimates are drawn as a pale tint of the same route colour,
    nested (honest) estimates as the solid colour.

    Null distributions are grey, the observed value is a solid black line,
    the 95th percentile of the null a dashed black line.
    P. aeruginosa is an artefact and is always hatched.

The palette is the Okabe-Ito colourblind-safe set. Adjacent-pair separation was
checked for protanopia, deuteranopia and tritanopia.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

DPI = 300

ROUTE = {
    "A1":  "#E69F00",
    "A1u": "#009E73",
    "A2":  "#0072B2",
    "B":   "#CC79A7",
}
ROUTE_LABEL = {"A1": "Route A1", "A1u": "Route A1u", "A2": "Route A2", "B": "Route B"}

NULL_FILL = "#D9D9D9"
NULL_EDGE = "#8C8C8C"
OBSERVED = "#1A1A1A"
PCT95 = "#1A1A1A"
INK = "#1A1A1A"
MUTED = "#5A5A5A"
GRID = "#D8D8D8"
ARTEFACT_HATCH = "///"

STRAINS = ["H.I.", "S.A.", "S.Pneu", "S.Pyo", "P.A."]
STRAIN_LABEL = {
    "H.I.": "H. influenzae", "S.A.": "S. aureus", "S.Pneu": "S. pneumoniae",
    "S.Pyo": "S. pyogenes", "P.A.": "P. aeruginosa †",
}
ARTEFACT = "P.A."


def tint(hex_colour, amount=0.62):
    """Pale version of a route colour, used for apparent estimates."""
    r, g, b = to_rgb(hex_colour)
    return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)


def rcparams():
    plt.rcParams.update({
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.labelsize": 11.5,
        "axes.titlesize": 12.5,
        "axes.titleweight": "bold",
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.7,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "legend.fontsize": 10.5,
        "legend.frameon": False,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def tidy(ax, grid_axis="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(True, axis=grid_axis, zorder=0)
        ax.set_axisbelow(True)


# Revision (October 2026, EW): Fig. 6 and Fig. 7 share one canvas and one horizontal crop, so that
# at the same printed width their y-axes lie on the same vertical lines. The span is the content
# width of Fig. 7 (inches, before the 0.1 in padding of the tight crop); each figure checks that
# none of its content lies outside it.
XSPAN_FIG6_FIG7 = (-0.1961129032258066, 6.268208333333334)   # measured from Fig. 7


def save(fig, stem, outdir=".", xspan=None):
    import os
    from matplotlib.transforms import Bbox
    # dpi is passed explicitly (not just via rcParams) so every PNG carries a pHYs
    # resolution chunk; a PNG without it is read by most viewers as 72 dpi, which
    # made Fig. 1 (1967 px wide) render at 694 mm instead of 166 mm (round 4 report).
    bbox = "tight"
    if xspan is not None:
        tb = fig.get_tightbbox(fig.canvas.get_renderer())
        assert tb.x0 >= xspan[0] - 1e-3 and tb.x1 <= xspan[1] + 1e-3, (stem, tb.x0, tb.x1, xspan)
        pad = plt.rcParams["savefig.pad_inches"]
        bbox = Bbox.from_extents(xspan[0] - pad, tb.y0 - pad, xspan[1] + pad, tb.y1 + pad)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, "%s.%s" % (stem, ext)), dpi=DPI, bbox_inches=bbox)
    plt.close(fig)
    print("  wrote %s.png and %s.pdf" % (stem, stem))


def stars(p):
    return "*" if (p is not None and p < 0.05) else ""

# Compound classes in the score displays. Deliberately not red/green, which is the
# worst pair for red-green colour vision deficiency.
CLASS = {"active": "#D55E00", "non-active": "#56B4E9"}
