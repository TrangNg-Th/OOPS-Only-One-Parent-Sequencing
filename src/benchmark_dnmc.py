#!/usr/bin/env python3
"""Score a final OOPS call set against an external ground-truth DNM table.

Standalone: the pipeline itself never reads this. It exists so a run whose
child has published truth data can be labelled TP/FP without giving OOPS any
knowledge of the second parent.

Usage:
  python benchmark_dnmc.py --calls <final_dnmc_*.tsv> --truth <truth.tsv> \
      --child NA12878 --parent NA12891 [--out-dir DIR] [--parent-role paternal]
"""
import argparse
import os
import sys

import pandas as pd

# CEPH 1463. Only used to decide which half of the truth table a run could
# possibly have detected; override with --parent-role for other pedigrees.
FATHER = {"NA12877": "NA12889", "NA12878": "NA12891",
          "NA12879": "NA12877", "NA12881": "NA12877", "NA12882": "NA12877",
          "NA12885": "NA12877", "NA12886": "NA12877"}
MOTHER = {"NA12877": "NA12890", "NA12878": "NA12892",
          "NA12879": "NA12878", "NA12881": "NA12878", "NA12882": "NA12878",
          "NA12885": "NA12878", "NA12886": "NA12878"}


def parent_role(child, parent, override):
    if override:
        return override
    if FATHER.get(child) == parent:
        return "paternal"
    if MOTHER.get(child) == parent:
        return "maternal"
    sys.exit(f"ERROR: cannot place {parent} as a parent of {child}; "
             "pass --parent-role paternal|maternal")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="final_dnmc_<child>-from-<parent>.tsv")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--child", required=True)
    ap.add_argument("--parent", required=True)
    ap.add_argument("--parent-role", choices=["paternal", "maternal"], default=None)
    ap.add_argument("--out-dir", default=None,
                    help="default: alongside --calls")
    ap.add_argument("--paper-paternal-rate", type=float, default=1.54e-08,
                    help="published paternal DNM rate to compare against "
                         "[1.54e-08]")
    ap.add_argument("--paper-maternal-rate", type=float, default=0.37e-08,
                    help="published maternal DNM rate to compare against "
                         "[0.37e-08]")
    ap.add_argument("--denum-dir", default=None,
                    help="dir holding callable_genome.txt and "
                         "<child>_LR_validated_hetc.bed; default: derived from "
                         "--calls (<calls dir>/mismatch_analysis/denum_calcul)")
    ap.add_argument("--exclude-chroms", default="chrX,chrY,chrM",
                    help="must match the pipeline run [chrX,chrY,chrM]")
    args = ap.parse_args()

    role = parent_role(args.child, args.parent, args.parent_role)
    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.calls))
    os.makedirs(out_dir, exist_ok=True)
    excluded = {c.strip() for c in args.exclude_chroms.split(",") if c.strip()}

    calls = pd.read_csv(args.calls, sep="\t", dtype={"chrom": str})
    if "chrom" not in calls.columns or "pos" not in calls.columns:
        # zero-candidate run: final file holds only the placeholder comment
        calls = pd.DataFrame({"chrom": [], "pos": [], "PS_child": []})
    calls = calls[calls["chrom"].astype(str).str.startswith("chr")]
    # Apply the same chromosome exclusion to the calls as to the truth table.
    # Filtering only one side turns any call on an excluded chromosome into a
    # guaranteed false positive, which is a scoring artefact, not a result.
    n_before_excl = len(calls)
    excluded_calls = calls[calls["chrom"].isin(excluded)]
    calls = calls[~calls["chrom"].isin(excluded)]
    if len(excluded_calls):
        print(f"NOTE: {n_before_excl - len(calls)} call(s) on excluded chroms "
              f"dropped from scoring: "
              + ", ".join(f"{c}:{p}" for c, p in
                          zip(excluded_calls["chrom"], excluded_calls["pos"])))
    calls["pos"] = calls["pos"].astype(int) if len(calls) else calls["pos"]

    truth = pd.read_csv(args.truth, sep="\t", dtype={"chromosome": str})
    truth = truth[truth["sample"] == args.child].copy()
    if truth.empty:
        sys.exit(f"ERROR: no truth rows for {args.child} in {args.truth}")
    truth["position"] = truth["position"].astype(int)
    truth = truth[~truth["chromosome"].isin(excluded)]

    # A run against one parent can only see mutations that arose on the
    # haplotype that parent transmitted. "cannot_determine" rows are counted
    # separately: they are detectable in principle but not attributable.
    detectable = truth[truth["inheritance"] == role]
    ambiguous = truth[truth["inheritance"] == "cannot_determine"]

    tkey = truth.set_index(["chromosome", "position"])
    labelled = []
    for _, r in calls.iterrows():
        k = (r["chrom"], r["pos"])
        if k in tkey.index:
            t = tkey.loc[k]
            if isinstance(t, pd.DataFrame):
                t = t.iloc[0]
            labelled.append({**r.to_dict(), "label": "TP",
                             "ref": t["ref"], "alt": t["alt"],
                             "variant_origin": t["variant_origin"],
                             "inheritance": t["inheritance"]})
        else:
            labelled.append({**r.to_dict(), "label": "FP",
                             "ref": ".", "alt": ".",
                             "variant_origin": ".", "inheritance": "."})
    lab = pd.DataFrame(labelled, columns=[
        "chrom", "pos", "PS_child", "label", "ref", "alt",
        "variant_origin", "inheritance"])

    found = set(zip(lab.loc[lab["label"] == "TP", "chrom"],
                    lab.loc[lab["label"] == "TP", "pos"]))
    missed = detectable[~detectable.apply(
        lambda r: (r["chromosome"], r["position"]) in found, axis=1)]

    lab_file = os.path.join(out_dir, f"benchmark_{args.child}-from-{args.parent}_labelled.tsv")
    miss_file = os.path.join(out_dir, f"benchmark_{args.child}-from-{args.parent}_missed.tsv")
    sum_file = os.path.join(out_dir, f"benchmark_{args.child}-from-{args.parent}_summary.txt")
    lab.to_csv(lab_file, sep="\t", index=False)
    missed.to_csv(miss_file, sep="\t", index=False)

    # Callable genome, recomputed exactly as final_summary() in main.sh does:
    #   callable = (qualified_snps / sampled_snps) * accessible_bases
    #   rate     = n_calls / callable
    denum_dir = args.denum_dir or os.path.join(
        os.path.dirname(os.path.abspath(args.calls)),
        "mismatch_analysis", "denum_calcul")
    callable_bp = None
    rate_note = ""
    try:
        cg = os.path.join(denum_dir, "callable_genome.txt")
        vals = {}
        with open(cg) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 2:
                    vals[parts[0]] = float(parts[1])
        accessible = vals["total_callable_bases"]
        sampled = vals["total_sampled_snps"]
        qual_bed = os.path.join(denum_dir, f"{args.child}_LR_validated_hetc.bed")
        with open(qual_bed) as f:
            qualified = sum(1 for _ in f)
        if sampled > 0:
            callable_bp = (qualified / sampled) * accessible
    except Exception as e:
        rate_note = f"  (callable genome unavailable: {e})"

    paper_rate = (args.paper_paternal_rate if role == "paternal"
                  else args.paper_maternal_rate)

    def rate(n):
        if not callable_bp:
            return float("nan")
        return n / callable_bp

    n_tp = int((lab["label"] == "TP").sum())
    n_fp = int((lab["label"] == "FP").sum())
    n_calls = len(lab)
    prec = n_tp / n_calls if n_calls else float("nan")
    rec = len(found & set(zip(detectable["chromosome"], detectable["position"]))) / len(detectable) \
        if len(detectable) else float("nan")

    lines = [
        "================ OOPS benchmark vs ground truth ================",
        f"Calls file               : {args.calls}",
        f"Truth file               : {args.truth}",
        f"Child / parent           : {args.child} / {args.parent} ({role})",
        f"Excluded chroms          : {','.join(sorted(excluded))}",
        "----------------------------------------------------------------",
        f"Calls made               : {n_calls}",
        f"  TP (in truth)          : {n_tp}",
        f"  FP (not in truth)      : {n_fp}",
        f"Precision                : {prec:.3f}",
        "----------------------------------------------------------------",
        f"Truth rows for child     : {len(truth)}",
        f"  detectable ({role:<9}): {len(detectable)}",
        f"  cannot_determine       : {len(ambiguous)}",
        f"  other parent           : {len(truth) - len(detectable) - len(ambiguous)}",
        f"Recall (vs detectable)   : {rec:.3f}",
        "----------------------------------------------------------------",
        f"Callable genome (bp)     : {callable_bp:.6e}" if callable_bp
        else f"Callable genome (bp)     : NA{rate_note}",
        f"Mutation rate            : {rate(n_calls):.6e}   ({n_calls} calls)",
        f"Published {role:<9} rate : {paper_rate:.6e}",
        f"  difference (ours-paper): {rate(n_calls) - paper_rate:+.6e}",
        f"  ratio (ours/paper)     : "
        + (f"{rate(n_calls) / paper_rate:.2f}x" if paper_rate else "n/a"),
        f"  percent of published   : "
        + (f"{100.0 * rate(n_calls) / paper_rate:.1f}%" if paper_rate else "n/a"),
        f"Truth rate, detectable   : {rate(len(detectable)):.6e}   "
        f"({len(detectable)} detectable DNMs over the same callable genome)",
        "----------------------------------------------------------------",
        "TP breakdown by truth annotation:",
    ]
    if n_tp:
        bd = lab[lab["label"] == "TP"].groupby(["variant_origin", "inheritance"]).size()
        for (vo, inh), n in bd.items():
            lines.append(f"  {vo:<12} {inh:<18} {n}")
    else:
        lines.append("  (none)")
    lines += [
        "----------------------------------------------------------------",
        f"Labelled calls -> {lab_file}",
        f"Missed truth   -> {miss_file}",
        "================================================================",
    ]
    out = "\n".join(lines)
    print(out)
    with open(sum_file, "w") as f:
        f.write(out + "\n")


if __name__ == "__main__":
    main()
