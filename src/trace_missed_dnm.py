#!/usr/bin/env python3
"""Trace, for one OOPS run, where each truth DNM of the child was lost."""
import sys, os, re, glob, json
import numpy as np, pandas as pd

RUN = sys.argv[1]                      # project dir
TRUTH = sys.argv[2]
OUT = sys.argv[3]

name = os.path.basename(RUN.rstrip('/'))
# OOPS_<platform>_<cov>_<child2>_<parent2>
m = re.match(r'OOPS_([A-Za-z]+)_([0-9]+x)_(\d+)_(\d+)$', name)
plat, cov, c2, p2 = m.groups()
sid = lambda d: d if len(d) > 3 else 'NA128' + d
CHILD, PARENT = sid(c2), sid(p2)
SEX = {'NA12889':'M','NA12891':'M','NA12877':'M','NA12890':'F','NA12892':'F','NA12878':'F'}
POO = {'M':'paternal','F':'maternal'}.get(SEX.get(PARENT), '?')

PV = os.path.join(RUN, f'{CHILD}_phasedvcf')
MA = os.path.join(PV, 'mismatch_analysis')

# ---- parameters from the most recent chain2b log ----
P = dict(min_dp=15, max_dp=50, gq=30, nvq=0.5, mmdiff=0.1, maxad=0,
         cluster=10000, homo=8, ref=None)
logs = sorted(glob.glob(os.path.join(PV, 'chain_2b_4', 'logs', 'chain2b_*.out')),
              key=os.path.getmtime)
if logs:
    txt = open(logs[-1], errors='ignore').read()
    def g(pat, cast, key):
        mm = re.search(pat, txt)
        if mm:
            P[key] = cast(mm.group(1))
    g(r'DP range:\s*(\d+)\s*-', int, 'min_dp')
    g(r'DP range:\s*\d+\s*-\s*(\d+)', int, 'max_dp')
    g(r'Genotype Quality threshold:\s*(\d+)', int, 'gq')
    g(r'Number of Variants Quantile threshold:\s*([\d.]+)', float, 'nvq')
    g(r'Mismatch Difference Minimum:\s*([\d.]+)', float, 'mmdiff')
    g(r'Max parent contradicting AD reads:\s*(\d+)', int, 'maxad')
    g(r'Cluster-rejection window \(bp\):\s*(\d+)', int, 'cluster')
    g(r'Min homopolymer run to reject:\s*(\d+)', int, 'homo')
    g(r'Reference for homopolymer mask:\s*(\S+)', str, 'ref')
    if 'contradicting AD reads' not in txt:
        P['maxad'] = 10**9; P['cluster'] = 0; P['homo'] = 0; P['legacy'] = True

# ---- truth ----
t = pd.read_csv(TRUTH, sep='\t')
t = t[t['sample'] == CHILD].copy()
t['pos'] = t['position'].astype(int)
EXCL = {'chrx', 'chry', 'chrm'}

# ---- per-site table ----
ps_f = os.path.join(PV, 'HTblocks', f'{PARENT}_{CHILD}_ps.tsv')
df = pd.read_csv(ps_f, sep='\t', dtype=str)
df['pos'] = df['pos'].astype(int)
bi = re.compile(r'^[01][/|][01]$')
gt_ok = df['GT_mom'].str.match(bi).fillna(False) & df['GT_child'].str.match(bi).fillna(False)
df_gt = df[gt_ok].copy()
df_gt = df_gt[~df_gt['chrom'].str.lower().isin(EXCL)]
unf = df_gt.copy()                                     # pre-DP/GQ snapshot

num = lambda s: pd.to_numeric(s, errors='coerce').fillna(0)
keep = (num(df_gt['DP_mom']).between(P['min_dp'], P['max_dp']) &
        num(df_gt['DP_child']).between(P['min_dp'], P['max_dp']) &
        (num(df_gt['GQ_mom']) >= P['gq']) & (num(df_gt['GQ_child']) >= P['gq']))
f = df_gt[keep].copy()
gm = f['GT_mom'].str.replace('|', '/', regex=False)
bad = (((gm == '0/0') & (num(f['ADalt_mom']) > P['maxad'])) |
       ((gm == '1/1') & (num(f['ADref_mom']) > P['maxad'])))
f = f[~bad]

def mmcounts(g):
    v = g['GT_child'].str.split(r'[/|]', expand=True).astype(int)
    h0, h1 = v[0].values, v[1].values
    mom = g['GT_mom'].str.split(r'[/|]', expand=True).astype(int).values
    return ~(mom == h0[:, None]).any(axis=1), ~(mom == h1[:, None]).any(axis=1)

blk = {}
for (ch, ps), g in f.groupby(['chrom', 'PS_child']):
    m0, m1 = mmcounts(g)
    blk[(ch, ps)] = dict(n_var=len(g), n0=int(m0.sum()), n1=int(m1.sum()),
                         blen=int(g['pos'].max() - g['pos'].min() + 1),
                         p0=set(g['pos'].values[m0]), p1=set(g['pos'].values[m1]))
sumdf = pd.DataFrame([dict(n_var=v['n_var']) for v in blk.values()])
nv_cut = sumdf['n_var'].quantile(P['nvq']) if len(sumdf) else 0

raw = {}
for (ch, ps), g in unf.groupby(['chrom', 'PS_child']):
    m0, m1 = mmcounts(g)
    raw[(ch, ps)] = (g['pos'].values[m0], g['pos'].values[m1])

fa = None
if P['ref'] and os.path.exists(P['ref']):
    try:
        import pysam; fa = pysam.FastaFile(P['ref'])
    except Exception:
        fa = None

def homopoly(ch, pos, flank=60):
    if fa is None or P['homo'] <= 0:
        return False
    st = max(0, pos - 1 - flank)
    try:
        seq = fa.fetch(ch, st, pos + flank).upper()
    except Exception:
        return False
    idx = pos - 1 - st
    for pr in (idx - 1, idx, idx + 1):
        if pr < 0 or pr >= len(seq) or seq[pr] not in 'ACGT':
            continue
        b = seq[pr]; lo = hi = pr
        while lo > 0 and seq[lo - 1] == b: lo -= 1
        while hi < len(seq) - 1 and seq[hi + 1] == b: hi += 1
        if hi - lo + 1 >= P['homo']:
            return True
    return False

def load(path, cols):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        try:
            d = pd.read_csv(path, sep='\t', dtype=str)
            return set(zip(d[cols[0]], d[cols[1]].astype(int)))
        except Exception:
            return set()
    return set()

dnmc = load(os.path.join(MA, f'{PARENT}_{CHILD}_dnmc.tsv'), ['#[1]CHROM', '[2]POS'])
if not dnmc:
    dnmc = load(os.path.join(MA, f'{PARENT}_{CHILD}_dnmc.tsv'), ['chrom', 'pos'])
pbam = os.path.join(MA, f'{PARENT}_{CHILD}_parent_BAM_check.tsv')
pbam = pbam if os.path.exists(pbam) else os.path.join(MA, f'{PARENT}_parent_BAM_check.tsv')
pb = pd.read_csv(pbam, sep='\t') if os.path.exists(pbam) else pd.DataFrame()
pb_fail = set()
if len(pb):
    q = pb[pb['decision'] != 'PASS']
    pb_fail = set(zip(q['chrom'], q['pos'].astype(int)))
lrb = os.path.join(MA, 'sliced', f'{CHILD}_LR_validated_dnmc.bed')
lr_ok = set()
if os.path.exists(lrb) and os.path.getsize(lrb) > 0:
    d = pd.read_csv(lrb, sep='\t', header=None, dtype=str)
    lr_ok = set(zip(d[0], d[2].astype(int)))
fin = os.path.join(PV, f'final_dnmc_{CHILD}-from-{PARENT}.tsv')
final = load(fin, ['chrom', 'pos'])

sites = f.set_index(['chrom', 'pos'])
unf_i = unf.set_index(['chrom', 'pos'])
gt_i = df_gt.set_index(['chrom', 'pos'])
all_i = df.set_index(['chrom', 'pos'])

rows = []
for _, r in t.iterrows():
    ch, pos = r['chromosome'], r['pos']
    rec = dict(run=name, platform=plat, coverage=cov, child=CHILD, parent=PARENT,
               chrom=ch, pos=pos, ref=r['ref'], alt=r['alt'],
               origin=r['variant_origin'], inheritance=r['inheritance'])
    poo = r['inheritance']
    detectable = (poo == POO) or (r['variant_origin'] == 'postzygotic') or (poo == 'cannot_determine')
    rec['parent_of_origin_detectable'] = detectable
    key = (ch, pos)

    if ch.lower() in EXCL:
        rec['stage'] = 'S0_excluded_chrom'
    elif not detectable:
        rec['stage'] = 'S0_wrong_parent_of_origin'
    elif key not in all_i.index:
        rec['stage'] = 'S1_absent_from_phased_merged_vcf'
    elif key not in gt_i.index:
        rec['stage'] = 'S2_genotype_not_biallelic_phased'
    else:
        g = gt_i.loc[key]
        g = g.iloc[0] if isinstance(g, pd.DataFrame) else g
        rec['GT_parent'], rec['GT_child'] = g['GT_mom'], g['GT_child']
        rec['DP_parent'], rec['DP_child'] = g['DP_mom'], g['DP_child']
        rec['GQ_parent'], rec['GQ_child'] = g['GQ_mom'], g['GQ_child']
        rec['AD_parent'] = f"{g['ADref_mom']},{g['ADalt_mom']}"
        rec['AD_child'] = f"{g['ADref_child']},{g['ADalt_child']}"
        gtc = str(g['GT_child'])
        if '|' not in gtc:
            rec['stage'] = 'S2_child_unphased'
        elif gtc not in ('0|1', '1|0'):
            rec['stage'] = 'S2_child_not_het'
        elif str(g['GT_mom']).replace('|', '/') != '0/0':
            rec['stage'] = 'S2_parent_not_homref'
        elif key not in sites.index:
            why = []
            if not (P['min_dp'] <= float(g['DP_mom'] or 0) <= P['max_dp']): why.append('DP_parent')
            if not (P['min_dp'] <= float(g['DP_child'] or 0) <= P['max_dp']): why.append('DP_child')
            if float(g['GQ_mom'] or 0) < P['gq']: why.append('GQ_parent')
            if float(g['GQ_child'] or 0) < P['gq']: why.append('GQ_child')
            if float(g['ADalt_mom'] or 0) > P['maxad']: why.append('parentAD_veto')
            rec['stage'] = 'S3_site_DPGQAD_filter'
            rec['detail'] = '+'.join(why) or 'unknown'
        else:
            s = sites.loc[key]
            s = s.iloc[0] if isinstance(s, pd.DataFrame) else s
            ps = s['PS_child']
            b = blk.get((ch, ps))
            rec['PS'] = ps; rec['block_n_var'] = b['n_var']
            rec['block_n_mm_h0'] = b['n0']; rec['block_n_mm_h1'] = b['n1']
            rec['nv_cutoff'] = nv_cut
            is0, is1 = pos in b['p0'], pos in b['p1']
            if not (is0 or is1):
                rec['stage'] = 'S4_not_a_haplotype_mismatch'
            else:
                hap, nsame, nother = (0, b['n0'], b['n1']) if is0 else (1, b['n1'], b['n0'])
                rec['hap'] = hap; rec['n_mm_same_hap'] = nsame; rec['n_mm_other_hap'] = nother
                mn, mx = min(b['n0'], b['n1']), max(b['n0'], b['n1'])
                mmd = (mn / mx) if mx > 0 else 0.0
                rec['mismatch_difference'] = round(mmd, 4)
                if b['n_var'] < nv_cut:
                    rec['stage'] = 'S4_block_too_few_variants'
                elif b['blen'] <= 1:
                    rec['stage'] = 'S4_block_length_1'
                elif nsame > 1:
                    rec['stage'] = 'S4_not_lone_mismatch_on_its_hap'
                elif nother <= 1:
                    rec['stage'] = 'S4_other_hap_not_gt1_mismatch'
                elif mmd >= P['mmdiff']:
                    rec['stage'] = 'S4_mismatch_difference_too_high'
                else:
                    others = raw.get((ch, ps), (np.array([]), np.array([])))[hap]
                    near = others[(np.abs(others - pos) <= P['cluster']) & (others != pos)]
                    if P['cluster'] > 0 and near.size:
                        rec['stage'] = 'S5_cluster_rejected'
                    elif homopoly(ch, pos):
                        rec['stage'] = 'S5_homopolymer_rejected'
                    elif key not in dnmc:
                        rec['stage'] = 'S5_dropped_at_candidate_step_other'
                    elif key in pb_fail:
                        rec['stage'] = 'S6_parent_BAM_alt_reads'
                    elif lr_ok and key not in lr_ok:
                        rec['stage'] = 'S7_longread_validation_failed'
                    elif key not in final:
                        rec['stage'] = 'S8_lost_in_20kb_rephase_recount'
                    else:
                        rec['stage'] = 'DETECTED'
    rows.append(rec)

out = pd.DataFrame(rows)
out.to_csv(OUT, sep='\t', index=False)
print(f'{name}\tdone\tn_truth={len(out)}\tdetected={(out["stage"]=="DETECTED").sum()}')
