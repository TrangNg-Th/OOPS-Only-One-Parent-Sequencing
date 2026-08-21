#!/usr/bin/env python3
"""
Why each ground-truth DNM was missed, per trio and per parent sex.

Runs src/trace_missed_dnm.py over every OOPS run directory (in parallel), then
folds the per-site results into one table:

    trio (child x sequenced parent) x parent sex x loss stage -> count, reason

Only germline autosomal truth DNMs whose parent of origin IS the sequenced
parent are counted: those are the ones OOPS could in principle see. A DNM that
arose in the un-sequenced parent leaves no signal for the method at all, so
counting it as a "miss" would just restate the study design.

Usage
-----
  python src/missed_dnm_report.py --out-dir logs/missed_dnm_analysis
  python src/missed_dnm_report.py --platform Pacbio --jobs 8
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
TRACER = os.path.join(HERE, "trace_missed_dnm.py")

DEFAULT_ROOTS = ["/N/scratch/nguyetrt/OOPS_outputs",
                 "/N/project/mutation_rate_Mmulatta/OOPS_outputs"]

FATHER = {"NA12877": "NA12889", "NA12878": "NA12891",
          "NA12879": "NA12877", "NA12882": "NA12877"}

# stage -> (pipeline step, plain-language reason)
STAGES = {
    "DETECTED": ("-", "recovered by OOPS"),
    "S0_excluded_chrom": ("n/a", "on chrX/chrY/chrM, excluded from both numerator and denominator"),
    "S0_wrong_parent_of_origin": ("n/a", "arose in the un-sequenced parent: no signal by design"),
    "S0b_on_unsequenced_parent_haplotype": ("n/a", "sits on the un-sequenced parent's haplotype (truth label disagrees with the data)"),
    "S1_absent_from_phased_merged_vcf": ("upstream of Part 2b", "site missing from the DRAGEN VCF, or no FORMAT/PS after WhatShap"),
    "S2_child_unphased": ("Part 2b", "child genotype present but not phased"),
    "S2_child_not_het": ("Part 2b", "child not called heterozygous"),
    "S2_parent_not_homref": ("Part 2b", "parent not called homozygous reference"),
    "S3_site_DPGQAD_filter": ("Part 2b count_mismatches.py", "site-level DP / GQ / parent-AD cut"),
    "S4_not_a_haplotype_mismatch": ("Part 2b count_mismatches.py", "site is not a mismatch against the parent at all"),
    "S4_block_too_few_variants": ("Part 2b count_mismatches.py", "PS block carries fewer variants than the n_variants quantile cutoff"),
    "S4_phase_noise_extra_mm_on_parent_hap": ("Part 2b count_mismatches.py", "not a lone mismatch: phase-switch noise adds mismatches on the parent-matching haplotype"),
    "S4_not_lone_mismatch_on_its_hap": ("Part 2b count_mismatches.py", "not a lone mismatch on its haplotype"),
    "S4_other_hap_not_gt1_mismatch": ("Part 2b count_mismatches.py", "the other haplotype has <=1 mismatch, so the block is uninformative"),
    "S4_mismatch_difference_too_high": ("Part 2b count_mismatches.py", "mismatch_difference above threshold: the two haplotypes are not asymmetric enough"),
    "S4_block_length_1": ("Part 2b count_mismatches.py", "PS block spans a single position"),
    "S5_cluster_rejected": ("Part 2b count_mismatches.py", "another mismatch on the same haplotype within the cluster window"),
    "S5_homopolymer_rejected": ("Part 2b count_mismatches.py", "candidate inside a long homopolymer run"),
    "S5_dropped_at_candidate_step_other": ("Part 2b count_mismatches.py", "dropped at the candidate step for another reason"),
    "S6_parent_BAM_alt_reads": ("Part 2b parent_readcheck.py", "ALT-supporting reads present in the parent Illumina BAM"),
    "S7_longread_validation_failed": ("Part 3 dnmc_readcheck.py", "ALT not confined to one haplotype, or too few supporting long reads"),
    "S8_lost_in_20kb_rephase_recount": ("Part 3b re-phase + recount", "survived Part 3 but failed the +-20 kb local re-phase recount"),
}

# Losses that are properties of the study design rather than of the pipeline.
BY_DESIGN = {"S0_excluded_chrom", "S0_wrong_parent_of_origin",
             "S0b_on_unsequenced_parent_haplotype"}


def run_dirs(roots, platform):
    out = []
    for root in roots:
        for d in sorted(glob.glob(os.path.join(root, "OOPS_*"))):
            name = os.path.basename(d)
            m = re.match(r"OOPS_([A-Za-z]+)_\d+x_\d+_\d+$", name)
            if not m or not os.path.isdir(d):
                continue
            if platform and m.group(1).lower() != platform.lower():
                continue
            out.append(d)
    return out


def trace_all(dirs, truth, tmp_dir, jobs):
    os.makedirs(tmp_dir, exist_ok=True)

    def one(d):
        tag = os.path.basename(d) + ("__project" if "/project/" in d else "__scratch")
        out = os.path.join(tmp_dir, tag + ".tsv")
        p = subprocess.run([sys.executable, TRACER, d, truth, out],
                           capture_output=True, text=True)
        if p.returncode != 0:
            last = p.stderr.strip().splitlines()[-1] if p.stderr.strip() else "?"
            print(f"  skip {os.path.basename(d)}: {last}", file=sys.stderr)
            return None
        return out

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        return [r for r in ex.map(one, dirs) if r]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outputs-root", action="append", default=None)
    ap.add_argument("--truth", default=os.path.join(REPO, "data", "truth_dnm_CEPH1463.tsv"))
    ap.add_argument("--platform", default="Pacbio")
    ap.add_argument("--out-dir", default=os.path.join(REPO, "logs", "missed_dnm_analysis"))
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--coverage", type=int, default=None,
                    help="restrict the summary to one coverage (e.g. 10)")
    args = ap.parse_args()

    roots = args.outputs_root or DEFAULT_ROOTS
    dirs = run_dirs(roots, args.platform)
    print(f"tracing {len(dirs)} run director{'y' if len(dirs)==1 else 'ies'}",
          file=sys.stderr)
    os.makedirs(args.out_dir, exist_ok=True)
    files = trace_all(dirs, args.truth, os.path.join(args.out_dir, "per_run"), args.jobs)
    if not files:
        sys.exit("no run could be traced")

    a = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)

    # Split the "not a lone mismatch" bucket: when the DNM's own haplotype is the
    # one carrying hundreds of mismatches, that haplotype is the un-sequenced
    # parent's, and the miss is by design rather than a filter effect.
    m = a.stage == "S4_not_lone_mismatch_on_its_hap"
    wrong = m & (a.n_mm_same_hap > a.n_mm_other_hap)
    a.loc[wrong, "stage"] = "S0b_on_unsequenced_parent_haplotype"
    a.loc[m & ~wrong, "stage"] = "S4_phase_noise_extra_mm_on_parent_hap"

    a["sex"] = ["paternal" if p == FATHER.get(c) else "maternal"
                for c, p in zip(a.child, a.parent)]
    a["trio"] = a.child + " from " + a.parent
    a.to_csv(os.path.join(args.out_dir, "per_site_all_runs.tsv"), sep="\t", index=False)

    # The detectable set: germline, autosomal, parent of origin == sequenced parent
    d = a[(a.origin == "germline")
          & (~a.chrom.str.lower().isin(["chrx", "chry", "chrm"]))
          & (a.inheritance == a.sex)].copy()
    if args.coverage:
        d = d[d.coverage == f"{args.coverage}x"]

    # one run per trio: the deepest completed one
    d["cov_n"] = d.coverage.str.rstrip("x").astype(int)
    pick = d.sort_values("cov_n").groupby(["trio", "sex"]).run.last()
    d = d[[r == pick.get((t, s)) for t, s, r in zip(d.trio, d.sex, d.run)]]

    rows = []
    for (trio, sex, run, stage), g in d.groupby(["trio", "sex", "run", "stage"]):
        step, reason = STAGES.get(stage, ("?", "?"))
        rows.append(dict(trio=trio, sex=sex, run=run, stage=stage,
                         pipeline_step=step, reason=reason, n=len(g),
                         by_design=stage in BY_DESIGN,
                         sites="; ".join(f"{c}:{p}" for c, p in zip(g.chrom, g.pos))))
    s = pd.DataFrame(rows)
    tot = d.groupby(["trio", "sex"]).size().rename("truth_n").reset_index()
    s = s.merge(tot, on=["trio", "sex"])
    s["pct_of_truth"] = (100 * s.n / s.truth_n).round(1)
    s = s.sort_values(["sex", "trio", "n"], ascending=[True, True, False])

    out = os.path.join(args.out_dir, "missed_by_trio_sex_stage.tsv")
    s.to_csv(out, sep="\t", index=False)

    pd.set_option("display.width", 400, "display.max_rows", 400,
                  "display.max_colwidth", 78)
    print(s[["trio", "sex", "truth_n", "stage", "pipeline_step", "n",
             "pct_of_truth", "reason"]].to_string(index=False))
    print(f"\nwrote {out}", file=sys.stderr)
    print(f"wrote {os.path.join(args.out_dir, 'per_site_all_runs.tsv')}", file=sys.stderr)


if __name__ == "__main__":
    main()
