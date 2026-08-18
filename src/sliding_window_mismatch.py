#!/usr/bin/env python3
"""
Sliding-window haplotype-vs-parent mismatch profile, per PS block.

For every phase set (PS_child) we walk a sliding window along the block and
count, inside each window, how many child variants have a haplotype allele
that is absent from the parent genotype (same mismatch definition as
count_mismatches.py). A correctly phased block has one haplotype (the
parent-inherited one) with ~0 mismatches all the way through; a phase-switch
error shows up as the two curves crossing / swapping roles mid-block.

Usage:
  python sliding_window_mismatch.py <input_ps_tsv> <out_dir> [options]

Options:
  --window INT      window size in bp (default 10000)
  --step INT        step size in bp (default = window/2, i.e. 50% overlap)
  --min-dp INT      min DP in both samples (default 5)
  --max-dp INT      max DP in both samples (default 100)
  --gq INT          min GQ in both samples (default 20)
  --min-variants    min variants in a PS block to keep it (default 100)
  --top INT         number of blocks to plot, largest first (default 20)
  --chrom STR       restrict to one chromosome
  --ps STR          restrict to one PS_child (implies --top on that block only)
  --no-plots        write the TSV only
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(add_help=True)
    p.add_argument("input_tsv")
    p.add_argument("out_dir")
    p.add_argument("--window", type=int, default=10000)
    p.add_argument("--step", type=int, default=None)
    p.add_argument("--min-dp", type=int, default=5)
    p.add_argument("--max-dp", type=int, default=100)
    p.add_argument("--gq", type=int, default=10)
    p.add_argument("--min-variants", type=int, default=5)
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--smooth", type=int, default=21,
                   help="rolling-median width (in windows) for the phase signal")
    p.add_argument("--chrom", default=None)
    p.add_argument("--ps", default=None)
    p.add_argument("--no-plots", action="store_true")
    a = p.parse_args()
    if a.step is None:
        a.step = max(1, a.window // 2)
    return a


def load(args):
    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str)
    if args.chrom:
        df = df[df["chrom"] == args.chrom]
    if args.ps:
        df = df[df["PS_child"] == args.ps]

    # same QC as count_mismatches.py
    df = df[
        df["DP_mom"].astype(int).between(args.min_dp, args.max_dp)
        & df["DP_child"].astype(int).between(args.min_dp, args.max_dp)
    ]
    df = df[
        (df["GQ_mom"].astype(int) >= args.gq)
        & (df["GQ_child"].astype(int) >= args.gq)
    ]
    biallelic01 = re.compile(r"^[01][\/|][01]$")
    df = df[
        df["GT_mom"].str.match(biallelic01) & df["GT_child"].str.match(biallelic01)
    ]
    # only phased child sites carry haplotype information
    df = df[df["GT_child"].str.contains(r"\|")]
    df = df[df["PS_child"] != "."]
    df["pos"] = df["pos"].astype(int)

    # per-site mismatch flags, computed once for the whole table
    ch = df["GT_child"].values.astype(str)
    mo = df["GT_mom"].values.astype(str)
    h0 = np.array([s[0] for s in ch])
    h1 = np.array([s[2] for s in ch])
    m0 = np.array([s[0] for s in mo])
    m1 = np.array([s[2] for s in mo])
    df["mm_h0"] = ~((m0 == h0) | (m1 == h0))
    df["mm_h1"] = ~((m0 == h1) | (m1 == h1))
    return df


def per_site_mismatch(g):
    """Return (mm_h0, mm_h1) boolean arrays: allele absent from parent GT."""
    return g["mm_h0"].values, g["mm_h1"].values


def window_profile(g, window, step):
    """Sliding-window counts over one PS block. g must be sorted by pos."""
    pos = g["pos"].values
    mm_h0, mm_h1 = per_site_mismatch(g)
    start0, end0 = pos[0], pos[-1]
    starts = np.arange(start0, max(start0 + 1, end0 - window + 1 + step), step)
    # cumulative sums for O(1) window queries
    c0 = np.concatenate([[0], np.cumsum(mm_h0)])
    c1 = np.concatenate([[0], np.cumsum(mm_h1)])
    li = np.searchsorted(pos, starts, side="left")
    ri = np.searchsorted(pos, starts + window, side="left")
    n = ri - li
    n0 = c0[ri] - c0[li]
    n1 = c1[ri] - c1[li]
    keep = n > 0
    return pd.DataFrame(
        {
            "win_start": starts[keep],
            "win_end": (starts + window)[keep],
            "win_mid": (starts + window // 2)[keep],
            "n_variants": n[keep],
            "n_mismatch_h0": n0[keep],
            "n_mismatch_h1": n1[keep],
            "frac_mismatch_h0": n0[keep] / n[keep],
            "frac_mismatch_h1": n1[keep] / n[keep],
        }
    )


def add_phase_signal(prof, smooth):
    """
    phase_signal = frac_mismatch_h1 - frac_mismatch_h0, rolling-median smoothed.
    Positive => H0 is the parent-matching haplotype in this region.
    A sign flip along the block is a phase-switch error.
    """
    d = prof["frac_mismatch_h1"] - prof["frac_mismatch_h0"]
    prof["delta"] = d
    prof["phase_signal"] = (
        d.rolling(smooth, center=True, min_periods=max(1, smooth // 3)).median()
    )
    s = np.sign(prof["phase_signal"]).replace(0, np.nan).ffill().bfill()
    prof["phase_state"] = s

    # what actually gets plotted: the fraction of variants in the window whose
    # allele is NOT found in the parent, for each haplotype, smoothed the same
    # way. The parent-inherited haplotype sits near 0, the other one high.
    roll = dict(window=smooth, center=True, min_periods=max(1, smooth // 3))
    m0 = prof["frac_mismatch_h0"].rolling(**roll).median()
    m1 = prof["frac_mismatch_h1"].rolling(**roll).median()
    prof["mismatch_frac_h0_smoothed"] = m0
    prof["mismatch_frac_h1_smoothed"] = m1
    # the better-matching haplotype in this window, and its mismatch fraction.
    # phase_signal = frac_mismatch_h1 - frac_mismatch_h0, so s > 0 means H1
    # (= h0) has FEWER mismatches and is therefore the parent-matching one.
    prof["best_hap"] = np.where(s > 0, "H1", "H2")
    prof["lowest_mismatch_frac"] = np.where(s > 0, m0, m1)
    return prof


def switch_intervals(prof):
    """Consecutive runs of phase_state; boundaries between runs are switches."""
    s = prof["phase_state"].values
    if len(s) == 0:
        return []
    brk = np.flatnonzero(s[1:] != s[:-1])
    out = []
    for i in brk:
        out.append(
            {
                "left_end": int(prof["win_end"].values[i]),
                "right_start": int(prof["win_start"].values[i + 1]),
                "from_state": "H1-matches" if s[i] > 0 else "H2-matches",
                "to_state": "H1-matches" if s[i + 1] > 0 else "H2-matches",
            }
        )
    return out


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    df = load(args)
    if df.empty:
        sys.exit("No sites left after filtering.")

    wk = f"{args.window // 1000}kb"
    rows = []
    blocks = []
    switches = []
    for (chrom, ps), g in df.groupby(["chrom", "PS_child"], sort=False):
        if g.shape[0] < args.min_variants:
            continue
        g = g.sort_values("pos")
        prof = add_phase_signal(window_profile(g, args.window, args.step), args.smooth)
        prof.insert(0, "chrom", chrom)
        prof.insert(1, "PS_child", ps)
        rows.append(prof)

        sw = switch_intervals(prof)
        for s in sw:
            s.update({"chrom": chrom, "PS_child": ps})
        switches.extend(sw)

        mm_h0, mm_h1 = per_site_mismatch(g)
        blocks.append(
            {
                "chrom": chrom,
                "PS_child": ps,
                "block_start": g["pos"].iloc[0],
                "block_end": g["pos"].iloc[-1],
                "block_length": g["pos"].iloc[-1] - g["pos"].iloc[0] + 1,
                "n_variants": g.shape[0],
                "n_mismatches_h0": int(mm_h0.sum()),
                "n_mismatches_h1": int(mm_h1.sum()),
                "n_switches": len(sw),
            }
        )

    if not rows:
        sys.exit("No PS block passed --min-variants.")

    prof_all = pd.concat(rows, ignore_index=True)
    blk = pd.DataFrame(blocks).sort_values("n_variants", ascending=False)

    prof_file = os.path.join(args.out_dir, f"sliding_mismatch.{wk}.tsv")
    blk_file = os.path.join(args.out_dir, f"sliding_mismatch.{wk}.blocks.tsv")
    sw_file = os.path.join(args.out_dir, f"sliding_mismatch.{wk}.switches.tsv")
    prof_all.to_csv(prof_file, sep="\t", index=False)
    blk.to_csv(blk_file, sep="\t", index=False)
    pd.DataFrame(
        switches,
        columns=["chrom", "PS_child", "left_end", "right_start",
                 "from_state", "to_state"],
    ).to_csv(sw_file, sep="\t", index=False)
    print(f"Wrote {prof_file} ({len(prof_all)} windows)")
    print(f"Wrote {blk_file} ({len(blk)} blocks)")
    print(f"Wrote {sw_file} ({len(switches)} candidate switch points)")

    if args.no_plots:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top = blk.head(args.top)
    for _, b in top.iterrows():
        p = prof_all[
            (prof_all["chrom"] == b["chrom"])
            & (prof_all["PS_child"] == b["PS_child"])
        ]
        fig, ax = plt.subplots(2, 1, figsize=(14, 7), sharex=True,
                               gridspec_kw={"height_ratios": [3, 1]})
        x = (p["win_mid"] / 1e6).values

        # panel 1: per-window mismatch fraction of each haplotype against the
        # parent, H1 upward and H2 downward on the same 0..1 scale. The
        # parent-inherited haplotype sits near 0 and the other one high, so a
        # phase-switch error is where the high side swaps.
        m0 = p["mismatch_frac_h0_smoothed"].values
        m1 = p["mismatch_frac_h1_smoothed"].values
        ax[0].axhline(0, color="k", lw=0.8)
        ax[0].fill_between(x, 0, m0, color="#1f77b4", label="H1 mismatch fraction")
        ax[0].fill_between(x, 0, -m1, color="#d62728", label="H2 mismatch fraction")
        ax[0].plot(x, p["frac_mismatch_h0"].values, color="0.25", lw=0.3, alpha=0.45)
        ax[0].plot(x, -p["frac_mismatch_h1"].values, color="0.25", lw=0.3,
                   alpha=0.45, label="raw (unsmoothed)")
        ax[0].set_ylabel(
            "Mismatch fraction\n"
            f"← H2      H1 →     ({wk} win, med{args.smooth})"
        )
        ax[0].set_ylim(-1.05, 1.05)
        ax[0].set_yticks([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
        ax[0].set_yticklabels(
            [f"{abs(t):.2f}".rstrip("0").rstrip(".") or "0"
             for t in ax[0].get_yticks()]
        )
        ax[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False,
                     fontsize=8, ncol=3)
        ax[0].set_title(
            f"{b['chrom']} PS={b['PS_child']},  "
            f"PS length: {b['block_length']/1e6:.2f} Mb, {b['n_variants']} variants, "
            f"Mismatches H1={b['n_mismatches_h0']} H2={b['n_mismatches_h1']}, "
            f"switches={b['n_switches']}",
            pad=22,
        )
        for s in [d for d in switches
                  if d["chrom"] == b["chrom"] and d["PS_child"] == b["PS_child"]]:
            for a in ax:
                a.axvline(s["left_end"] / 1e6, color="orange", lw=1.0, ls="--", label="phase switch")

        
        # panel 2: number of variants per window
        ax[1].plot(x, p["n_variants"], color="grey", lw=0.8)
        ax[1].set_ylabel("Variants\nper window")
        ax[1].set_xlabel(f"{b['chrom']} position (Mb)")
        fig.tight_layout()
        out = os.path.join(
            args.out_dir, f"sliding_mismatch.{wk}.{b['chrom']}_PS{b['PS_child']}.png"
        )
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
