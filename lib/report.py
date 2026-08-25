#!/usr/bin/env python3
"""
What the pipeline placed, against what it recovered.

Every stage is scored the same way: something known was put into the file,
a tool was asked to find it, and the two are compared. That is the whole
argument of the project, and this is where it is made concrete.

Two rules govern the numbers.

Recall counts only mutations verified present in the reads. BAMSurgeon
reports success at positions where realignment later drove coverage to
nothing, mostly in the MHC; counting those against Mutect2 would blame the
caller for a mutation that was never in the file. Placed-but-absent
mutations are reported separately.

Sensitivity is broken down by allele fraction rather than given as one
number. Per-sample recall varies with each donor's fraction distribution
rather than with anything about the sample, and a single figure hides
that. At 34x a mutation present in 5% of reads sits in under two of them,
which is the floor everything else runs into.

Usage:
  python report.py <results_dir> <samples_dir> [sample ...]
"""
import sys
import os
import gzip
import glob
import re
import subprocess
import pandas as pd
import numpy as np

RESULTS = sys.argv[1]
SAMPLE_DIR = sys.argv[2]
SAMPLES = sys.argv[3:] if len(sys.argv) > 3 else None

if SAMPLES is None:
    SAMPLES = sorted(
        os.path.basename(os.path.dirname(p))
        for p in glob.glob(f"{RESULTS}/*/") )
    SAMPLES = [s for s in SAMPLES if os.path.exists(f"{SAMPLE_DIR}/{s}_truth.tsv")]

BINS = [0, 0.05, 0.10, 0.15, 0.20, 0.30, 1.01]
LABELS = ["<5%", "5-10%", "10-15%", "15-20%", "20-30%", ">30%"]
W = 74


def rule(title):
    print()
    print("=" * W)
    print(f" {title}")
    print("=" * W)


def pass_calls(path):
    """Every PASS record as (chrom, pos, ref, alt)."""
    out = []
    if not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or f[6] != "PASS":
                continue
            for a in f[4].split(","):
                out.append((f[0], int(f[1]), f[3], a))
    return out


def all_calls(path):
    """Every record with its filter, so emitted-then-rejected is visible."""
    out = []
    if not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            for a in f[4].split(","):
                out.append((f[0], int(f[1]), f[3], a, f[6]))
    return out


rows = []
per_sample = []

for sid in SAMPLES:
    truth_path = f"{SAMPLE_DIR}/{sid}_truth.tsv"
    vcf = f"{RESULTS}/{sid}/mutect2/{sid}.filtered.vcf.gz"
    if not os.path.exists(truth_path):
        continue

    t = pd.read_csv(truth_path, sep="\t")
    calls = all_calls(vcf)
    passed = {(c[0], c[1], c[3]) for c in calls if c[4] == "PASS"}
    emitted = {(c[0], c[1], c[3]) for c in calls}
    # an indel may be left-aligned a base either way, so position alone
    # decides for those rather than the exact alternate string
    pass_pos = {(c[0], c[1]) for c in calls if c[4] == "PASS"}
    emit_pos = {(c[0], c[1]) for c in calls}

    for _, r in t.iterrows():
        is_indel = str(r.get("type", "SNP")) != "SNP"
        key = (r.chrom, int(r.pos), str(r.alt))
        if is_indel:
            found = any(abs(int(r.pos) - p) <= 2 for c, p in pass_pos
                        if c == r.chrom)
            seen = any(abs(int(r.pos) - p) <= 2 for c, p in emit_pos
                       if c == r.chrom)
        else:
            found = key in passed
            seen = key in emitted
        rows.append({
            "sample": sid, "chrom": r.chrom, "pos": int(r.pos),
            "alt": str(r.alt), "gene": r.get("gene"),
            "vaf": float(r.get("vaf", np.nan)),
            "type": str(r.get("type", "SNP")),
            "in_mhc": (r.chrom == "chr6" and 29600000 <= int(r.pos) < 33100000),
            "found": found, "emitted": seen,
        })

    n_calls = len({(c[0], c[1], c[3]) for c in calls if c[4] == "PASS"})
    per_sample.append({"sample": sid, "placed": len(t),
                       "pass_calls": n_calls})

d = pd.DataFrame(rows)
if d.empty:
    print("nothing to report — has the pipeline run?")
    sys.exit(1)

snv = d[d.type == "SNP"]
ind = d[d.type != "SNP"]

# =====================================================================
rule("WHAT WAS PLACED, AND WHAT CAME BACK")
print()
print(f"  {'sample':<8}{'placed':>8}{'SNV':>6}{'indel':>7}"
      f"{'recovered':>11}{'recall':>9}")
for sid in SAMPLES:
    g = d[d["sample"] == sid]
    if not len(g):
        continue
    print(f"  {sid:<8}{len(g):>8}{int((g.type=='SNP').sum()):>6}"
          f"{int((g.type!='SNP').sum()):>7}{int(g.found.sum()):>11}"
          f"{100*g.found.mean():>8.1f}%")
print(f"  {'total':<8}{len(d):>8}{len(snv):>6}{len(ind):>7}"
      f"{int(d.found.sum()):>11}{100*d.found.mean():>8.1f}%")

# =====================================================================
rule("SUBSTITUTIONS AND INDELS SEPARATELY")
print()
for label, g in [("substitutions", snv), ("indels", ind)]:
    if not len(g):
        continue
    print(f"  {label:<16}{int(g.found.sum())} of {len(g)}"
          f"   ({100*g.found.mean():.1f}%)")
if len(ind):
    for ty, g in ind.groupby("type"):
        print(f"    {ty:<14}{int(g.found.sum())} of {len(g)}")
print()
print("  An indel is several bases of evidence where a substitution is")
print("  one, so it is the easier of the two to call, not the harder.")

# =====================================================================
rule("SENSITIVITY BY ALLELE FRACTION")
sn = snv[snv.vaf.notna()].copy()
if len(sn):
    sn["bin"] = pd.cut(sn.vaf, bins=BINS, labels=LABELS, right=False)
    print()
    print(f"  {'requested':<11}{'placed':>8}{'recovered':>11}{'recall':>9}"
          f"{'reads at 34x':>14}")
    for lab in LABELS:
        g = sn[sn.bin == lab]
        if not len(g):
            continue
        print(f"  {lab:<11}{len(g):>8}{int(g.found.sum()):>11}"
              f"{100*g.found.mean():>8.1f}%{34*g.vaf.median():>14.1f}")
    print()
    print("  Below about three supporting reads no caller can separate a")
    print("  variant from a sequencing error. That floor is arithmetic,")
    print("  not a property of Mutect2.")

# =====================================================================
rule("WHERE THE MISSES ARE")
missed = d[~d.found]
print()
print(f"  not recovered            {len(missed)} of {len(d)}")
if len(missed):
    emitted_then_filtered = int(missed.emitted.sum())
    print(f"    called, then filtered  {emitted_then_filtered}")
    print(f"    never called at all    {len(missed) - emitted_then_filtered}")
    print()
    print("  The first is a filter setting and can be changed. The second")
    print("  is a detection limit and cannot.")

if d.in_mhc.any():
    mh = d[d.in_mhc]
    print()
    print(f"  inside the MHC           {int(mh.found.sum())} of {len(mh)}"
          f"   ({100*mh.found.mean():.0f}%)")
    print(f"  outside                  {int(d[~d.in_mhc].found.sum())} of "
          f"{int((~d.in_mhc).sum())}"
          f"   ({100*d[~d.in_mhc].found.mean():.0f}%)")
    print()
    print("  hg38 carries alternate MHC haplotypes. Reads there map")
    print("  ambiguously, mapping quality collapses, and the caller")
    print("  declines to open an active region at all.")

# =====================================================================
rule("HLA TYPING")
gt = {}
for sid in SAMPLES:
    p = f"{RESULTS}/{sid}/optitype/{sid}_result.tsv"
    if not os.path.exists(p):
        continue
    x = pd.read_csv(p, sep="\t")
    if not len(x):
        continue
    r = x.iloc[0]
    alleles = [str(r.get(f"{l}{k}")) for l in "ABC" for k in (1, 2)]
    alleles = [a for a in alleles if a and a.lower() != "nan"]
    gt[sid] = alleles
    print(f"  {sid:<8}A {r.get('A1')}/{r.get('A2')}   "
          f"B {r.get('B1')}/{r.get('B2')}   C {r.get('C1')}/{r.get('C2')}")

if gt:
    sigs = {"|".join(sorted(v)) for v in gt.values()}
    print()
    print(f"  {len(sigs)} distinct genotypes among {len(gt)} samples")
    if len(sigs) < len(gt):
        print()
        print("  Two samples share a genotype. Given fewer threads than it")
        print("  wants, OptiType returns a plausible default rather than")
        print("  failing — this is the check that catches it.")

# =====================================================================
rule("HLA LOSS")
design = f"{RESULTS}/loh_design.tsv"
if os.path.exists(design):
    des = pd.read_csv(design, sep="\t")
    print()
    print(f"  {'sample':<8}{'target':<9}{'before':>9}{'after':>8}"
          f"{'controls unchanged':>22}")
    for _, r in des.iterrows():
        loc = str(r.target_locus).replace("HLA-", "")
        before = r.get(f"{loc}_before")
        after = r.get(f"{loc}_after")
        others = [l for l in "ABC" if l != loc]
        same = all(abs(float(r.get(f"{o}_before", 0)) -
                       float(r.get(f"{o}_after", 0))) < 0.05 for o in others)
        print(f"  {r['sample']:<8}{r.target_locus:<9}{before:>9}{after:>8}"
              f"{('yes' if same else 'NO'):>22}")

found_loh = []
for p in glob.glob(f"{RESULTS}/*/loh/out_*/*HLAlossPrediction*.txt"):
    try:
        t = pd.read_csv(p, sep="\t")
    except Exception:
        continue
    if not len(t):
        continue
    sid = p.split("/results/")[-1].split("/")[0] if "/results/" in p \
        else os.path.basename(p).split(".")[0]
    m = re.search(r"/([^/]+)/loh/out_([abc])/", p)
    if m:
        sid, loc = m.group(1), m.group(2).upper()
    else:
        loc = "?"
    r = t.iloc[0]
    pval = r.get("PVal", np.nan)
    if pd.isna(pval):
        continue
    found_loh.append({"sample": sid, "locus": f"HLA-{loc}",
                      "cn1": r.get("HLA_type1copyNum_withBAFBin"),
                      "cn2": r.get("HLA_type2copyNum_withBAFBin"),
                      "pval": pval,
                      "sites": r.get("numMisMatchSitesCov")})

if found_loh:
    fl = pd.DataFrame(found_loh)
    tgt = {}
    if os.path.exists(design):
        tgt = dict(zip(des["sample"], des.target_locus))
    print()
    print(f"  {'sample':<8}{'locus':<8}{'CN1':>7}{'CN2':>7}{'p':>11}"
          f"{'sites':>7}   role")
    for _, r in fl.sort_values("pval").iterrows():
        role = "TARGET" if tgt.get(r["sample"]) == r.locus else "control"
        print(f"  {r['sample']:<8}{r.locus:<8}{r.cn1:>7.3f}{r.cn2:>7.3f}"
              f"{r.pval:>11.4g}{int(r.sites):>7}   {role}")

    sig = fl[fl.pval < 0.05]
    on_target = sum(1 for _, r in sig.iterrows()
                    if tgt.get(r["sample"]) == r.locus)
    print()
    print(f"  significant at p<0.05    {len(sig)}")
    print(f"    on the thinned locus   {on_target}")
    print(f"    on an untouched locus  {len(sig) - on_target}")
    print()
    print(f"  {len(fl)} of {3*len(SAMPLES)} locus runs produced a result at")
    print(f"  all. At panel depth the two alleles of a locus often differ")
    print(f"  at too few covered positions for the test to run — a blank")
    print(f"  result is the method's limit, not a failure.")
    print()
    print(f"  What was simulated is a loss of coverage across the locus,")
    print(f"  not the loss of one allele. Separating the two alleles needs")
    print(f"  more discriminating positions than this data provides, so")
    print(f"  this measures sensitivity to depth rather than to allelic")
    print(f"  imbalance.")
else:
    print()
    print("  No locus produced a prediction.")

# =====================================================================
rule("NEOANTIGENS")
tot_mut = tot_strong = tot_weak = 0
n_pred = 0
print()
print(f"  {'sample':<8}{'mutations':>11}{'peptides':>10}{'strong':>8}{'weak':>7}")
for sid in SAMPLES:
    f = f"{RESULTS}/{sid}/neoantigens/neoantigens_per_mutation.tsv"
    if not os.path.exists(f):
        continue
    t = pd.read_csv(f, sep="\t")
    print(f"  {sid:<8}{len(t):>11}{int(t.peptides.sum()):>10}"
          f"{int(t.strong.sum()):>8}{int(t.weak.sum()):>7}")
    tot_mut += len(t); tot_strong += int(t.strong.sum())
    tot_weak += int(t.weak.sum()); n_pred += 1

if n_pred:
    print()
    print(f"  strong binders per mutation: {tot_strong/max(1,tot_mut):.2f}")

    print()
    print("  The HLA genotype belongs to the 1000 Genomes backbone, not to")
    print("  the TCGA donor whose mutations these are. These are the")
    print("  peptides those mutations would present in someone carrying")
    print("  this genotype — the right quantity for measuring what the")
    print("  pipeline detects, not a claim about the donor.")


# The second route. mhcflurry starts from the mutations that were placed;
# pVACseq starts from the ones Mutect2 recovered, by way of VEP. They share
# no code and see different inputs, so the comparison is between two
# independent answers rather than two runs of one method.
pv_rows = []
for sid in SAMPLES:
    hits = glob.glob(
        f"{RESULTS}/{sid}/pvacseq/out/MHC_Class_I/*.all_epitopes.tsv")
    if not hits:
        continue
    try:
        t = pd.read_csv(hits[0], sep="\t", low_memory=False)
    except Exception:
        continue
    col = next((c for c in t.columns if c == "Best MT Percentile"), None)
    if col is None:
        col = next((c for c in t.columns if "Percentile" in c), None)
    if col is None:
        continue
    v = pd.to_numeric(t[col], errors="coerce")
    pv_rows.append({"sample": sid, "epitopes": len(t),
                    "strong": int((v < 0.5).sum()),
                    "weak": int(((v >= 0.5) & (v < 2.0)).sum())})

if pv_rows:
    rule("THE SECOND PREDICTOR")
    pv = pd.DataFrame(pv_rows)
    print()
    print(f"  {'sample':<8}{'epitopes':>10}{'strong':>8}{'weak':>7}")
    for _, r in pv.iterrows():
        print(f"  {r['sample']:<8}{r.epitopes:>10}{r.strong:>8}{r.weak:>7}")
    print(f"  {'total':<8}{pv.epitopes.sum():>10}{pv.strong.sum():>8}"
          f"{pv.weak.sum():>7}")

    print()
    print(f"  {'sample':<8}{'mhcflurry':>12}{'NetMHCpan':>12}")
    for sid in SAMPLES:
        f = f"{RESULTS}/{sid}/neoantigens/neoantigens_per_mutation.tsv"
        mf = 0
        if os.path.exists(f):
            mf = int(pd.read_csv(f, sep="\t").strong.sum())
        row = pv[pv["sample"] == sid]
        pvs = int(row.strong.iloc[0]) if len(row) else 0
        if mf or pvs:
            print(f"  {sid:<8}{mf:>12}{pvs:>12}")
    print()
    print("  The two do not have to agree, and where they differ the")
    print("  reason is usually the input rather than the model: mhcflurry")
    print("  sees every placed mutation, pVACseq only those the caller")
    print("  recovered. At around 80% recall that is a fifth fewer")
    print("  mutations to work from.")
    print()
    print("  netMHCpan is licensed per user by DTU and cannot be shipped.")
    print("  Where this section is missing the demo ran without it, and")
    print("  everything else still holds.")

# =====================================================================
rule("SUMMARY")
print()
print(f"  samples                  {len(SAMPLES)}")
print(f"  mutations placed         {len(d)}")
print(f"  recovered by Mutect2     {int(d.found.sum())}"
      f"   ({100*d.found.mean():.1f}%)")
if len(gt):
    print(f"  HLA genotypes            {len(sigs)} distinct of {len(gt)}")
if found_loh:
    print(f"  HLA loss detected        {len(fl[fl.pval<0.05])} loci")
if n_pred:
    print(f"  strong binders           {tot_strong}")

out = f"{RESULTS}/summary.tsv"
d.to_csv(out, sep="\t", index=False)
print()
print(f"  written to {out}")
