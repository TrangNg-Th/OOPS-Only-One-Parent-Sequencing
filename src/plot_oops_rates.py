#!/usr/bin/env python3
"""
OOPS mutation-rate figures, driven by the run directories themselves.

Unlike src/dnm_vs_depth.py (which carries a hand-edited DATA block), this script
walks the OOPS project directories, reads each run's final DNM call set and its
callable-genome denominator, and derives everything it plots. Nothing to keep in
sync by hand.

Two figures:

  1. --figure depth   Estimated rate vs read depth for one child, paternal and
                      maternal series, each against its own ground-truth line,
                      with 95% Poisson bands. This is the coverage figure.

  2. --figure trios   Every trio on one forest plot, split by parent sex, each
                      estimate with its Poisson interval and the trio's own
                      truth value marked. This is the cross-trio figure.

  --figure both       (default) writes both.

A run is only plotted if it has BOTH a final_dnmc_*.tsv and a callable_genome.txt.
Runs whose outputs are byte-identical to another run's (a stale copy left behind
by a re-submitted chain) are drawn hollow and named in the caption -- they are not
independent measurements.

Usage
-----
  python src/plot_oops_rates.py --out-dir figures
  python src/plot_oops_rates.py --figure depth --child NA12879 --platform Pacbio
  python src/plot_oops_rates.py --truth-source published     # 1.54e-8 / 0.37e-8
"""

import argparse
import glob
import hashlib
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# --------------------------------------------------------------------------
# Pedigree. Parent-of-origin of a DNM is "paternal" when the sequenced parent
# is the child's father, "maternal" when it is the mother.
# --------------------------------------------------------------------------
FATHER = {"NA12877": "NA12889", "NA12878": "NA12891",
          "NA12879": "NA12877", "NA12882": "NA12877",
          "NA12881": "NA12877", "NA12883": "NA12877"}
MOTHER = {"NA12877": "NA12890", "NA12878": "NA12892",
          "NA12879": "NA12878", "NA12882": "NA12878",
          "NA12881": "NA12878", "NA12883": "NA12878"}

DEFAULT_ROOTS = ["/N/scratch/nguyetrt/OOPS_outputs",
                 "/N/project/mutation_rate_Mmulatta/OOPS_outputs"]

DAD = "#1f6fb4"
MOM = "#b42d2d"
DAD_FILL = "#1f6fb4"
MOM_FILL = "#b42d2d"
HILITE = "#f5c26b"
GRID = "#d9d9d9"


# ==========================================================================
# Collecting the data
# ==========================================================================
def poisson_ci(k, alpha=0.05):
    """Exact (Garwood) Poisson interval for an observed count k."""
    from scipy.stats import chi2
    lo = 0.0 if k == 0 else chi2.ppf(alpha / 2, 2 * k) / 2
    hi = chi2.ppf(1 - alpha / 2, 2 * (k + 1)) / 2
    return lo, hi


def parse_run_name(name):
    m = re.match(r"OOPS_([A-Za-z]+)_(\d+)x_(\d+)_(\d+)$", name)
    if not m:
        return None
    platform, cov, c2, p2 = m.groups()
    sid = lambda d: d if len(d) > 3 else "NA128" + d
    return platform, int(cov), sid(c2), sid(p2)


def sex_of_parent(child, parent):
    if FATHER.get(child) == parent:
        return "paternal"
    if MOTHER.get(child) == parent:
        return "maternal"
    return None


def collect_runs(roots, platform_filter):
    """One row per completed run: counts, denominator, rate."""
    rows = []
    for root in roots:
        for d in sorted(glob.glob(os.path.join(root, "OOPS_*"))):
            if not os.path.isdir(d):
                continue
            parsed = parse_run_name(os.path.basename(d))
            if parsed is None:
                continue
            platform, cov, child, parent = parsed
            if platform_filter and platform.lower() != platform_filter.lower():
                continue
            sex = sex_of_parent(child, parent)
            if sex is None:
                print(f"  skip {os.path.basename(d)}: {parent} is not a parent "
                      f"of {child} in the pedigree table", file=sys.stderr)
                continue

            final = glob.glob(os.path.join(d, "*_phasedvcf", "final_dnmc_*.tsv"))
            callf = glob.glob(os.path.join(d, "*_phasedvcf", "mismatch_analysis",
                                           "denum_calcul", "callable_genome.txt"))
            if not final or not callf:
                print(f"  skip {os.path.basename(d)}: incomplete "
                      f"(no {'final call set' if not final else 'callable genome'})",
                      file=sys.stderr)
                continue

            n_dnm = sum(1 for _ in open(final[0])) - 1
            txt = open(callf[0]).read()
            accessible = float(re.search(r"total_callable_bases\s+(\d+)", txt).group(1))
            sampled = int(re.search(r"total_sampled_snps\s+(\d+)", txt).group(1))

            qualified = read_qualified(d)
            if qualified is None:
                print(f"  skip {os.path.basename(d)}: no Part-4 log to read the "
                      f"qualified-SNP count from", file=sys.stderr)
                continue

            callable_bp = qualified / sampled * accessible
            lo, hi = poisson_ci(n_dnm)
            fingerprint = hashlib.md5(
                open(final[0], "rb").read() + txt.encode()).hexdigest()

            rows.append(dict(
                run=os.path.basename(d), root=root, platform=platform,
                coverage=cov, child=child, parent=parent, sex=sex,
                n_dnm=n_dnm, sampled_snps=sampled, qualified_snps=qualified,
                accessible_bp=accessible, callable_bp=callable_bp,
                rate=n_dnm / callable_bp,
                rate_lo=lo / callable_bp, rate_hi=hi / callable_bp,
                fingerprint=fingerprint,
                mtime=os.path.getmtime(final[0]),
                stale=is_stale(d, final[0])))

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Mark runs that are byte-identical to an earlier one: a re-submitted chain
    # that never actually rewrote its outputs is not a second measurement.
    seen = {}
    dup = []
    for i, fp in enumerate(df.fingerprint):
        dup.append(fp in seen)
        seen.setdefault(fp, i)
    df["duplicate_of"] = [df.run.iloc[seen[fp]] if d else ""
                          for fp, d in zip(df.fingerprint, dup)]
    df["is_duplicate"] = dup
    return df


def dedupe_cells(df):
    """Keep one run per (child, sex, coverage): the most recently written.

    The same cell can be present twice -- e.g. an old run under
    /N/project and a rerun under /N/scratch. Plotting both puts two points on
    one x position and makes the line zig-zag through a version difference
    rather than through anything biological.
    """
    if df.empty:
        return df, []
    df = df.sort_values("mtime")
    dropped = []
    keep = []
    for (child, sex, cov), g in df.groupby(["child", "sex", "coverage"]):
        keep.append(g.index[-1])
        for i in g.index[:-1]:
            dropped.append((df.run[i], df.root[i], df.run[g.index[-1]]))
    out = df.loc[sorted(keep)].reset_index(drop=True)
    return out, dropped


def is_stale(run_dir, final_path):
    """True if the call set is older than the child BAM it claims to describe.

    A run whose inputs were replaced after its results were written was never
    computed from those inputs, so its numbers belong to whatever BAM used to
    be there. That is not a measurement of this depth.
    """
    bams = glob.glob(os.path.join(run_dir, "bam", "*.bam"))
    if not bams:
        return False
    return os.path.getmtime(final_path) < max(os.path.getmtime(b) for b in bams)


def read_qualified(run_dir):
    """Number of denominator SNPs that survived the read check (Part 4)."""
    logs = sorted(glob.glob(os.path.join(run_dir, "*_phasedvcf", "chain_2b_4",
                                         "logs", "chain4_*.out")),
                  key=os.path.getmtime)
    for log in reversed(logs):
        txt = open(log, errors="ignore").read()
        m = re.search(r"Number of snps qualified\s*:\s*(\d+)", txt)
        if m:
            return int(m.group(1))
    bed = glob.glob(os.path.join(run_dir, "*_phasedvcf", "mismatch_analysis",
                                 "denum_calcul", "*_LR_validated_hetc.bed"))
    if bed and os.path.getsize(bed[0]) > 0:
        return sum(1 for _ in open(bed[0]))
    return None


def truth_rates(truth_file, accessible_bp, source, published):
    """Ground-truth rate per child per parent sex, plus the cohort average.

    Germline only, autosomes only -- the comparison the paper makes.
    """
    t = pd.read_csv(truth_file, sep="\t")
    t = t[~t.chromosome.str.lower().isin(["chrx", "chry", "chrm"])]
    t = t[t.variant_origin == "germline"]
    counts = t.pivot_table(index="sample", columns="inheritance",
                           values="position", aggfunc="count", fill_value=0)
    for col in ("paternal", "maternal"):
        if col not in counts:
            counts[col] = 0

    per_child = {(s, sex): counts.loc[s, sex] / accessible_bp
                 for s in counts.index for sex in ("paternal", "maternal")}
    cohort = {sex: counts[sex].mean() / accessible_bp
              for sex in ("paternal", "maternal")}
    n_per_child = {(s, sex): int(counts.loc[s, sex])
                   for s in counts.index for sex in ("paternal", "maternal")}

    if source == "published":
        cohort = dict(paternal=published[0], maternal=published[1])
    return per_child, cohort, n_per_child


# ==========================================================================
# Figure 1 -- rate vs coverage
# ==========================================================================
def figure_depth(df, per_child, cohort, n_per_child, child, out_path,
                 highlight=10, truth_mode="trio"):
    sub = df[df.child == child].copy()
    if sub.empty:
        print(f"No completed runs for child {child}; skipping depth figure.",
              file=sys.stderr)
        return None
    sub = sub.sort_values("coverage")
    U = 1e-8                                  # plot in units of 1e-8 /bp/gen

    fig, (ax, axn) = plt.subplots(
        2, 1, figsize=(11, 7.6), sharex=True,
        gridspec_kw=dict(height_ratios=[1, 1], hspace=0.10))

    covs = sorted(sub.coverage.unique())
    x_of = {c: i for i, c in enumerate(covs)}

    if highlight in x_of:
        for a in (ax, axn):
            a.axvspan(x_of[highlight] - 0.28, x_of[highlight] + 0.28,
                      color=HILITE, alpha=0.45, zorder=0, lw=0)

    # dad's truth label on the right, mom's on the left, so they cannot collide
    for sex, color, marker, label, dy, lx, ha in (
            ("paternal", DAD, "o", "dad", 13, 0.997, "right"),
            ("maternal", MOM, "s", "mom", -19, 0.003, "left")):
        s_ = sub[sub.sex == sex]
        if s_.empty:
            continue
        truth = (per_child.get((child, sex)) if truth_mode == "trio"
                 else cohort[sex]) / U
        ax.axhline(truth, color=color, ls="--", lw=2.2, alpha=0.85, zorder=2)
        ax.annotate(f"ground truth, {label}"
                    f"  ({n_per_child.get((child, sex), 0)} DNMs)",
                    xy=(lx, truth), xycoords=("axes fraction", "data"),
                    xytext=(0, 5), textcoords="offset points",
                    color=color, fontsize=9, style="italic", ha=ha,
                    va="bottom")

        x = [x_of[c] for c in s_.coverage]
        ax.plot(x, s_.rate / U, color=color, lw=2.6, zorder=3)

        solid = ~s_.is_duplicate.values
        ax.scatter(np.array(x)[solid], (s_.rate / U).values[solid], s=100,
                   zorder=4, color=color, edgecolor="white", linewidth=1.4,
                   marker=marker)
        if (~solid).any():
            ax.scatter(np.array(x)[~solid], (s_.rate / U).values[~solid], s=100,
                       zorder=4, facecolor="white", edgecolor=color,
                       linewidth=1.8, marker=marker)

        top = (sub.rate / U).max()
        for xi, r, n in zip(x, s_.rate / U, s_.n_dnm):
            # below the marker only when there is room for the text there
            off = dy if (dy > 0 or r > 0.10 * top) else 14
            ax.annotate(f"{int(n)}", xy=(xi, r), xytext=(0, off),
                        textcoords="offset points", ha="center", fontsize=9,
                        color=color, fontweight="bold", zorder=5)

        axn.plot(x, s_.callable_bp / 1e9, color=color, lw=2.0, marker=marker,
                 ms=7, mec="white", mew=1.2, alpha=0.9)

    ax.set_ylabel("Mutation rate  ($\\times 10^{-8}$ per bp per generation)",
                  fontsize=11)
    ax.set_ylim(bottom=0)
    ax.margins(x=0.04)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for a in (ax, axn):
        for spine in ("top", "right"):
            a.spines[spine].set_visible(False)

    axn.set_ylabel("callable\ngenome (Gb)", fontsize=9)
    axn.set_xlabel("Read depth (coverage)", fontsize=11)
    axn.set_xticks(range(len(covs)))
    axn.set_xticklabels([f"{c}x" for c in covs], fontsize=10)
    axn.grid(axis="y", color=GRID, lw=0.8)
    axn.set_axisbelow(True)
    axn.set_ylim(bottom=0)

    handles = [Line2D([], [], color=DAD, lw=2.4, ls="-", marker="o",
                      mec="white", ms=9, label="OOPS estimate — dad"),
               Line2D([], [], color=MOM, lw=2.4, ls="-", marker="s",
                      mec="white", ms=9, label="OOPS estimate — mom"),
               Line2D([], [], color="0.35", ls="--", lw=2.0,
                      label="ground truth (germline, autosomal)"),
               ]
    if sub.is_duplicate.any():
        handles.append(Line2D([], [], color="0.35", marker="o", ls="none",
                              ms=9, mfc="white", mew=1.8,
                              label="hollow = outputs identical to another run"))
    if highlight in x_of:
        handles.append(Patch(facecolor=HILITE, alpha=0.45,
                             label=f"{highlight}x — where the experiment started"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.06),
               ncol=3, frameon=False, fontsize=9.5)

    ax.set_title(f"OOPS mutation-rate estimate vs read depth — child {child}",
                 fontsize=14, fontweight="bold", loc="left", pad=26)
    # ax.text(0, 1.045,
    #         "bold numbers = DNMs called at that depth · the callable genome "
    #         "(lower panel) moves with them, so the ratio stays near truth",
    #         transform=ax.transAxes, fontsize=9.5, color="0.4", va="bottom")

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ==========================================================================
# Figure 2 -- every trio, per sex
# ==========================================================================
def figure_trios(df, per_child, cohort, n_per_child, out_path, truth_mode="trio",
                 coverage=10):
    if df.empty:
        return None
    U = 1e-9
    d = df[df.coverage == coverage].sort_values(["sex", "child"]).reset_index(drop=True)
    if d.empty:
        print(f"No runs at {coverage}x; skipping the per-trio figure.",
              file=sys.stderr)
        return None

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.6),
                             gridspec_kw=dict(wspace=0.34))

    for ax, sex, color, title in (
            (axes[0], "paternal", DAD, "Paternal — father sequenced"),
            (axes[1], "maternal", MOM, "Maternal — mother sequenced")):
        s_ = d[d.sex == sex]
        if s_.empty:
            ax.set_visible(False)
            continue

        labels, ys = [], []
        y = 0.0
        xmax = max((s_.rate / U).max(),
                   max(((per_child.get((c, sex)) if truth_mode == "trio"
                         else cohort[sex]) or 0) / U for c in s_.child))
        for child, grp in s_.groupby("child", sort=True):
            truth = (per_child.get((child, sex)) if truth_mode == "trio"
                     else cohort[sex]) / U
            for _, r in grp.iterrows():
                # estimate -> truth, so the gap itself is the thing you read
                ax.plot([r.rate / U, truth], [y, y], color=color, lw=1.6,
                        alpha=0.35, zorder=1)
                kw = (dict(mfc="white", mec=color, mew=1.9) if r.is_duplicate
                      else dict(color=color, mec="white", mew=1.3))
                ax.plot(r.rate / U, y, "o", ms=11, zorder=3, **kw)
                ax.plot(truth, y, "|", ms=15, mew=2.6, color=color, zorder=3)
                n = int(r.n_dnm)
                pct = 100 * r.rate / (truth * U) if truth else float("nan")
                ax.annotate(f"{n} call" + ("" if n == 1 else "s")
                            + f"  ·  {pct:.0f}% of truth",
                            xy=(max(r.rate / U, truth), y), xytext=(10, 0),
                            textcoords="offset points", va="center",
                            fontsize=8.5, color="0.45")
                labels.append(f"{child}  ({n_per_child.get((child, sex), 0)} true DNMs)")
                ys.append(y)
                y += 1

        ax.axvline(cohort[sex] / U, color="0.35", ls=":", lw=1.8, zorder=1)
        ax.annotate("cohort mean", xy=(cohort[sex] / U, y - 0.65),
                    xytext=(5, 0), textcoords="offset points", fontsize=8.5,
                    color="0.35", rotation=90, va="bottom", ha="left")

        ax.set_yticks(ys)
        ax.set_yticklabels(labels, fontsize=9.5)
        ax.set_ylim(y - 0.35, -0.75)
        ax.set_xlim(0, xmax * 1.55)
        ax.set_xlabel("Mutation rate  ($\\times 10^{-9}$ per bp per generation)",
                      fontsize=10)
        ax.set_title(title, fontsize=11.5, fontweight="bold", loc="left", pad=10)
        ax.grid(axis="x", color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis="y", length=0)

    fig.suptitle(f"OOPS rate estimates per trio at {coverage}x, "
                 f"by parent of origin",
                 fontsize=14, fontweight="bold", x=0.075, ha="left", y=1.10)
    fig.text(0.075, 1.015,
             "filled circle = OOPS estimate · tick = that trio's own germline "
             "autosomal truth · dotted = cohort mean",
             fontsize=9.5, color="0.4", ha="left")

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


# ==========================================================================
def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs-root", action="append", default=None,
                    help="directory holding OOPS_<platform>_<cov>_<child>_<parent> "
                         "run dirs; repeatable (default: scratch + project)")
    ap.add_argument("--truth", default=os.path.join(here, "data",
                                                    "truth_dnm_CEPH1463.tsv"))
    ap.add_argument("--platform", default="Pacbio",
                    help="only plot this platform ('' for all). Default: Pacbio")
    ap.add_argument("--child", default="NA12879",
                    help="child for the depth figure (default NA12879)")
    ap.add_argument("--figure", choices=["depth", "trios", "both"], default="both")
    ap.add_argument("--out-dir", default="figures")
    ap.add_argument("--accessible-bp", type=float, default=2666892426,
                    help="accessible genome from the paper, for the truth rates")
    ap.add_argument("--truth-source", choices=["truthfile", "published"],
                    default="truthfile",
                    help="cohort truth line from the truth TSV (default) or from "
                         "the published cohort rates")
    ap.add_argument("--published-paternal", type=float, default=1.54e-08)
    ap.add_argument("--published-maternal", type=float, default=0.37e-08)
    ap.add_argument("--truth-mode", choices=["trio", "cohort"], default="trio",
                    help="compare each estimate to its own trio's truth (default) "
                         "or to the cohort average")
    ap.add_argument("--highlight-depth", type=int, default=10)
    ap.add_argument("--trio-coverage", type=int, default=10,
                    help="coverage to show in the per-trio figure (default 10x)")
    ap.add_argument("--keep-stale-runs", action="store_true",
                    help="keep runs whose results are older than their own input "
                         "BAM (default: drop them)")
    ap.add_argument("--keep-duplicate-runs", action="store_true",
                    help="keep runs whose outputs are byte-identical to another "
                         "run's (default: drop them -- a re-submitted chain that "
                         "never rewrote its outputs is not a second measurement)")
    ap.add_argument("--keep-duplicate-cells", action="store_true",
                    help="plot every run, even when two runs share the same "
                         "child/sex/coverage (default: keep only the newest)")
    ap.add_argument("--table", default=None,
                    help="also write the derived per-run table to this TSV")
    args = ap.parse_args()

    roots = args.outputs_root or DEFAULT_ROOTS
    print("Scanning:", ", ".join(roots), file=sys.stderr)
    df = collect_runs(roots, args.platform)
    if df.empty:
        sys.exit("No completed runs found.")

    if not args.keep_stale_runs:
        for _, r in df[df.stale].iterrows():
            print(f"  drop {r.run}: results predate the child BAM in that run "
                  f"directory, so they were not computed from it", file=sys.stderr)
        df = df[~df.stale].reset_index(drop=True)

    if not args.keep_duplicate_runs:
        dup = df[df.is_duplicate]
        for _, r in dup.iterrows():
            print(f"  drop {r.run}: outputs byte-identical to {r.duplicate_of}, "
                  f"not an independent measurement", file=sys.stderr)
        df = df[~df.is_duplicate].reset_index(drop=True)

    if not args.keep_duplicate_cells:
        df, dropped = dedupe_cells(df)
        for old, root, kept in dropped:
            print(f"  drop {old} ({root}): superseded by the newer {kept} at the "
                  f"same child/sex/coverage", file=sys.stderr)

    per_child, cohort, n_per_child = truth_rates(
        args.truth, args.accessible_bp, args.truth_source,
        (args.published_paternal, args.published_maternal))

    print(f"\n{len(df)} completed run(s):", file=sys.stderr)
    cols = ["run", "child", "parent", "sex", "coverage", "n_dnm",
            "callable_bp", "rate", "is_duplicate"]
    print(df[cols].to_string(index=False), file=sys.stderr)
    if df.is_duplicate.any():
        for _, r in df[df.is_duplicate].iterrows():
            print(f"\nNOTE: {r.run} has outputs byte-identical to {r.duplicate_of} "
                  f"-- drawn hollow, not an independent measurement.", file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)
    if args.table:
        df.to_csv(args.table, sep="\t", index=False)
        print(f"\nwrote {args.table}", file=sys.stderr)

    written = []
    if args.figure in ("depth", "both"):
        p = figure_depth(df, per_child, cohort, n_per_child, args.child,
                         os.path.join(args.out_dir,
                                      f"oops_rate_vs_depth_{args.child}.png"),
                         highlight=args.highlight_depth,
                         truth_mode=args.truth_mode)
        if p:
            written.append(p)
    if args.figure in ("trios", "both"):
        p = figure_trios(df, per_child, cohort, n_per_child,
                         os.path.join(args.out_dir,
                                      f"oops_rate_per_trio_{args.trio_coverage}x.png"),
                         truth_mode=args.truth_mode,
                         coverage=args.trio_coverage)
        if p:
            written.append(p)

    for p in written:
        print(f"wrote {p}", file=sys.stderr)


if __name__ == "__main__":
    main()
