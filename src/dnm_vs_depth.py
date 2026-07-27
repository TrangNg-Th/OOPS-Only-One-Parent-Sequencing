#!/usr/bin/env python3
"""
De novo mutation calling vs read depth (PacBio).

One figure, three stacked panels sharing the read-depth x-axis:
  Panel 1 : # mutation calls vs depth   -- dad (blue), mom (red)
  Panel 2 : # false positives vs depth  -- dad (blue), mom (red)
  Panel 3 : estimated mutation rate vs depth, per parent, each with its own
            dashed ground-truth line and a 95% Poisson band (count CI / callable).

A soft-yellow vertical band marks 10x (where the experiment started).

EDIT ONLY THE "DATA" BLOCK BELOW.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

# ===========================================================================
# DATA  -- fill these in.  All arrays must line up with DEPTHS, same length.
#          Use np.nan for any (depth, parent) you haven't run yet.
# ===========================================================================

DEPTHS = np.array([5, 7, 10, 15, 37])              # x-axis (coverage)

# --- Panel 1: number of de novo calls, by parent of origin ---
CALLS_DAD = np.array([8, 8, 12, 13, 16])          # paternal DNM calls
CALLS_MOM = np.array([1, 2, 3, 3, 3])    # maternal DNM calls



# --- Panel 2: number of false positives, by parent ---
FALSE_POS_DAD = np.array([0, 0, 0, 0, 0])
FALSE_POS_MOM = np.array([0, 0, 0, 0, 0])

# --- Panel 3: estimated mutation rate, by parent (per bp/gen) ---
RATE_EST_DAD = np.array([1.658423e-08, 1.28e-8, 1.7e-8, 1.28e-8, 1.50e-8])
RATE_EST_MOM = np.array([5.894931e-09, 0.31e-8, 0.38e-8, 0.28e-8, 0.28e-08])


# Callable / accessible genome size per depth, PER PARENT (drives Poisson band).
# (From your pipeline's "Callable genome size" line for each parent run.)
# >>> PLACEHOLDER NUMBERS -- replace with your real per-parent values <<<
# CALLABLE_DAD = np.array([1.50e8, 2.00e8, 2.75e8, 3.20e8, 4.00e8])
# CALLABLE_MOM = np.array([1.40e8, 1.90e8, 2.60e8, np.nan, np.nan])

# --- Ground truth, per parent ---
GROUND_TRUTH_RATE_DAD = 1.54e-8
GROUND_TRUTH_RATE_MOM = 0.37e-8

CONF = 0.95
DATA_TYPE = "PacBio HiFi"   # for title only

# Where the experiment started (yellow vertical band).
START_DEPTH = 10

# ===========================================================================
# Colors  (mom/dad consistent across ALL panels, incl. their ground truth)
# ===========================================================================
C_DAD       = "#185FA5"   # blue  -- dad line, dad ground truth, dad band
C_MOM       = "#A32D2D"   # red   -- mom line, mom ground truth, mom band
C_DAD_BAND  = "#85B7EB"   # light blue band fill (dad)
C_MOM_BAND  = "#F09595"   # light red  band fill (mom)
C_START     = "#FAC775"   # soft yellow
C_GRID      = "#D3D1C7"

# ===========================================================================
# Poisson interval helper (exact, chi-square / Garwin form)
# ===========================================================================
def poisson_ci(k, conf=0.95):
    """Exact two-sided Poisson CI for an observed count k (scalar)."""
    a = 1.0 - conf
    lo = 0.0 if k == 0 else 0.5 * chi2.ppf(a / 2.0, 2 * k)
    hi = 0.5 * chi2.ppf(1.0 - a / 2.0, 2 * k + 2)
    return lo, hi

def rate_band(calls, callable_genome, conf=0.95):
    """Per-parent rate band: Poisson CI on the count / callable genome.
    Returns (lo, hi) arrays, np.nan where calls or callable are missing."""
    lo = np.full(len(calls), np.nan)
    hi = np.full(len(calls), np.nan)
    for i, (k, c) in enumerate(zip(calls, callable_genome)):
        if np.isnan(k) or np.isnan(c) or c <= 0:
            continue
        l, u = poisson_ci(int(round(k)), conf)
        lo[i], hi[i] = l / c, u / c
    return lo, hi

# ===========================================================================
# Plot
# ===========================================================================
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.edgecolor": "#888780",
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": C_GRID,
    "grid.linewidth": 0.6,
    "grid.alpha": 0.6,
})

fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1, figsize=(8.5, 10), sharex=True,
    gridspec_kw={"height_ratios": [1.1, 0.8, 1.1], "hspace": 0.12},
)

xpos = np.arange(len(DEPTHS))   # categorical spacing so 5/7/10/15/37 don't crush
m_dad = ~np.isnan(CALLS_DAD)
m_mom = ~np.isnan(CALLS_MOM)
mr_dad = ~np.isnan(RATE_EST_DAD)
mr_mom = ~np.isnan(RATE_EST_MOM)

def start_band(ax):
    xc = np.interp(START_DEPTH, DEPTHS, xpos)
    ax.axvspan(xc - 0.18, xc + 0.18, color=C_START, alpha=0.55, zorder=0)

# ---- Panel 1: calls ----
# start_band(ax1)
ax1.plot(xpos[m_dad], CALLS_DAD[m_dad], "-o", color=C_DAD, lw=2, ms=7,
         label="Paternal", zorder=3)
ax1.plot(xpos[m_mom], CALLS_MOM[m_mom], "-s", color=C_MOM, lw=2, ms=7,
         label="Maternal", zorder=3)
ax1.set_yticks(np.arange(0, max(CALLS_DAD) + 5, 2))
ax1.set_ylabel("De novo calls")
ax1.legend(frameon=False, loc="upper left", fontsize=10)
ax1.set_title(f"De novo mutation calling vs {DATA_TYPE} read depth",
              fontsize=13, fontweight="medium", loc="left", pad=10)

# ---- Panel 2: false positives ----
# start_band(ax2)
mfp_dad = ~np.isnan(FALSE_POS_DAD)
mfp_mom = ~np.isnan(FALSE_POS_MOM)
ax2.plot(xpos[mfp_dad], FALSE_POS_DAD[mfp_dad], "-o", color=C_DAD, lw=2, ms=7,
         label="False positives (dad)", zorder=3)
ax2.plot(xpos[mfp_mom], FALSE_POS_MOM[mfp_mom], "-s", color=C_MOM, lw=2, ms=7,
         label="False positives (mom)", zorder=3)
ax2.set_ylabel("False positives")
ax2.legend(frameon=False, loc="upper right", fontsize=10)
ax2.set_ylim(bottom=0)

# ---- Panel 3: rates, per parent, each with ground truth ----
# start_band(ax3)

# ground-truth lines, matching parent colors
ax3.axhline(GROUND_TRUTH_RATE_DAD, color=C_DAD, lw=2, ls="--", zorder=2,
            label="Ground truth (dad)")
ax3.axhline(GROUND_TRUTH_RATE_MOM, color=C_MOM, lw=2, ls="--", zorder=2,
            label="Ground truth (mom)")

# estimated rate lines
ax3.plot(xpos[mr_dad], RATE_EST_DAD[mr_dad], "-o", color=C_DAD, lw=2, ms=7,
         zorder=3, label="Estimated rate (dad)")
ax3.plot(xpos[mr_mom], RATE_EST_MOM[mr_mom], "-s", color=C_MOM, lw=2, ms=7,
         zorder=3, label="Estimated rate (mom)")

ax3.set_ylabel("Mutation rate (per bp/gen)")
ax3.set_ylim(bottom=0)
ax3.legend(frameon=False, loc="best", fontsize=9)

# ---- shared x-axis cosmetics ----
ax3.set_xticks(xpos)
ax3.set_xticklabels([f"{int(d)}x" for d in DEPTHS])
ax3.set_xlabel("Read depth (coverage)")

xc = np.interp(START_DEPTH, DEPTHS, xpos)
# ax1.annotate("experiment\nstart (10x)", xy=(xc, ax1.get_ylim()[1]),
            #  xytext=(xc + 0.25, ax1.get_ylim()[1] * 0.92),
            #  fontsize=9, color="#854F0B", ha="left", va="top")

# for ax in (ax1, ax2, ax3):
#     ax.margins(x=0.04)
#     ax.spines["top"].set_visible(False)
#     ax.spines["right"].set_visible(False)

# ===========================================================================
# Save
# ===========================================================================
out_dir = os.path.join(os.getcwd(), "plots")
os.makedirs(out_dir, exist_ok=True)
fig.savefig(os.path.join(out_dir, f"dnm_vs_depth_{DATA_TYPE}.png"), dpi=200, bbox_inches="tight")
fig.savefig(os.path.join(out_dir, f"dnm_vs_depth_{DATA_TYPE}.pdf"), bbox_inches="tight")

print("saved to", out_dir)