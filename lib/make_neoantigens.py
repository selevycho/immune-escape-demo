#!/usr/bin/env python3
"""
Predict neoantigens from the injected mutations using the sample's own
HLA genotype.

Ported from the cohort pipeline. The only change is the column names: the
demo truth file uses lower-case short names, so both spellings are
accepted and the rest of the logic is untouched.

Chain:
  truth file  ->  HGVSp_Short (e.g. p.R130Q)
              ->  the GENCODE isoform whose residue at that position
                  actually matches the wild-type amino acid
              ->  substitute, then slide a window over the mutated
                  position emitting every k-mer that contains it
              ->  mhcflurry, using the alleles OptiType called
              ->  binding table, flagged by strength

Two details that are easy to get wrong:

Isoform choice. TCGA numbers residues against whichever transcript its
annotation pipeline picked, which is not always the longest one — PTEN
p.R130Q lands on Q in the long isoform and on R in the canonical one. So
every isoform of the gene is checked and the first whose residue matches
HGVSp_Short is used. Mutations where no isoform matches are reported and
skipped, never guessed at.

mhcflurry input. When given a file, mhcflurry-predict rejects --alleles;
the allele has to be a column in the file itself, one row per
peptide-allele pair.

One caveat worth stating plainly: the HLA genotype belongs to the 1000
Genomes backbone, while the mutations belong to a TCGA patient. These
predictions are what these mutations would produce in someone with this
genotype — the right quantity for measuring what a pipeline detects, and
not a claim about the donor's own immune landscape.

Usage:
  python make_neoantigens.py <truth_tsv> <optitype_result_tsv> \
                             <gencode_pc_translations_fa_gz> <out_dir> [lengths]
"""
import sys, os, gzip, re, subprocess
import pandas as pd
import numpy as np

TRUTH = sys.argv[1]
OPTI = sys.argv[2]
PROT_FA = sys.argv[3]
OUT_DIR = sys.argv[4]
LENGTHS = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 \
    else [8, 9, 10, 11]

STRONG_PCT = 0.5      # presentation percentile below this = strong binder
WEAK_PCT = 2.0        # below this = weak binder

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

os.makedirs(OUT_DIR, exist_ok=True)

print("[1/7] reading truth set ...", flush=True)
truth = pd.read_csv(TRUTH, sep="\t", low_memory=False)

# the cohort pipeline writes MAF column names; the demo package writes
# short ones. Accepting both means this file works against either.
ALIASES = {
    "hgvsp": "HGVSp_Short", "gene": "Hugo_Symbol",
    "consequence": "Variant_Classification", "type": "Variant_Type",
    "vaf": "VAF", "chrom": "Chromosome_hg38", "pos": "Start_Position_hg38",
    "ref": "Reference_Allele", "alt": "Tumor_Seq_Allele2",
}
truth = truth.rename(columns={k: v for k, v in ALIASES.items()
                              if k in truth.columns and v not in truth.columns})

if "HGVSp_Short" not in truth.columns:
    raise SystemExit("truth set has no HGVSp_Short column")

mis = truth[(truth.Variant_Classification == "Missense_Mutation") &
            truth.HGVSp_Short.notna()].copy()
print("      total mutations    : %d" % len(truth), flush=True)
print("      missense with HGVSp: %d" % len(mis), flush=True)

print("[2/7] reading HLA genotype ...", flush=True)
opti = pd.read_csv(OPTI, sep="\t")
allele_cols = [c for c in ["A1", "A2", "B1", "B2", "C1", "C2"]
               if c in opti.columns]
raw = [str(opti.iloc[0][c]) for c in allele_cols]
alleles = sorted(set("HLA-" + a for a in raw if a and a != "nan"))
print("      called : %s" % ", ".join(raw), flush=True)
print("      unique : %s" % ", ".join(alleles), flush=True)

homozygous = [c for c in ["A", "B", "C"]
              if opti.iloc[0].get(c + "1") == opti.iloc[0].get(c + "2")]
if homozygous:
    print("      homozygous loci: %s" % ", ".join(homozygous), flush=True)

print("[3/7] loading all protein isoforms ...", flush=True)
# GENCODE header: >ENSP|ENST|ENSG|OTTG|OTTT|TRANSCRIPT-NAME|GENE_NAME|LENGTH
isoforms = {}
gene = None
tname = None
buf = []


def store(g, t, seq):
    if g is None or not seq:
        return
    isoforms.setdefault(g, []).append((t, seq))


with gzip.open(PROT_FA, "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            store(gene, tname, "".join(buf))
            parts = line[1:].strip().split("|")
            gene = parts[6] if len(parts) > 6 else None
            tname = parts[5] if len(parts) > 5 else "NA"
            buf = []
        else:
            buf.append(line.strip())
    store(gene, tname, "".join(buf))

for g in isoforms:
    isoforms[g].sort(key=lambda x: -len(x[1]))
n_iso = sum(len(v) for v in isoforms.values())
print("      genes: %d   isoforms: %d" % (len(isoforms), n_iso), flush=True)

print("[4/7] applying substitutions ...", flush=True)
pat1 = re.compile(r"^p\.([A-Z])(\d+)([A-Z])$")
pat3 = re.compile(r"^p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$")

records = []
skipped = []
iso_stats = {"longest": 0, "other": 0}

for _, r in mis.iterrows():
    g = r.Hugo_Symbol
    hg = str(r.HGVSp_Short).strip()

    m = pat1.match(hg)
    if m:
        wt, pos, mt = m.group(1), int(m.group(2)), m.group(3)
    else:
        m3 = pat3.match(hg)
        if not m3:
            skipped.append((g, hg, "unparsable"))
            continue
        wt = AA3TO1.get(m3.group(1).capitalize())
        mt = AA3TO1.get(m3.group(3).capitalize())
        pos = int(m3.group(2))
        if wt is None or mt is None:
            skipped.append((g, hg, "unknown amino acid"))
            continue

    cands = isoforms.get(g)
    if not cands:
        skipped.append((g, hg, "gene not in GENCODE proteins"))
        continue

    chosen = None
    for i, (tn, seq) in enumerate(cands):
        if pos <= len(seq) and seq[pos - 1] == wt:
            chosen = (tn, seq, i)
            break

    if chosen is None:
        longest = cands[0][1]
        got = longest[pos - 1] if pos <= len(longest) else "-"
        skipped.append((g, hg,
                        "no isoform has %s at %d (longest has %s, %d checked)"
                        % (wt, pos, got, len(cands))))
        continue

    tn, seq, idx = chosen
    iso_stats["longest" if idx == 0 else "other"] += 1
    mutseq = seq[:pos - 1] + mt + seq[pos:]

    records.append({"gene": g, "hgvsp": hg, "transcript": tn,
                    "wt": wt, "mt": mt, "pos": pos, "prot_len": len(seq),
                    "wt_seq": seq, "mut_seq": mutseq,
                    "VAF": r.VAF,
                    "chrom": r.Chromosome_hg38,
                    "genomic_pos": r.Start_Position_hg38})

print("      substitutions applied: %d" % len(records), flush=True)
print("        matched longest isoform: %d" % iso_stats["longest"], flush=True)
print("        matched another isoform: %d" % iso_stats["other"], flush=True)
print("      skipped: %d" % len(skipped), flush=True)
for s in skipped[:12]:
    print("        %s %s - %s" % s, flush=True)
if skipped:
    pd.DataFrame(skipped, columns=["gene", "hgvsp", "reason"]).to_csv(
        "%s/skipped_mutations.tsv" % OUT_DIR, sep="\t", index=False)

if not records:
    raise SystemExit("nothing to predict")

print("[5/7] generating peptides ...", flush=True)
peps = []
for rec in records:
    p = rec["pos"]
    ms = rec["mut_seq"]
    ws = rec["wt_seq"]
    for L in LENGTHS:
        lo = max(1, p - L + 1)
        hi = min(p, len(ms) - L + 1)
        for start in range(lo, hi + 1):
            mut_pep = ms[start - 1:start - 1 + L]
            wt_pep = ws[start - 1:start - 1 + L]
            if len(mut_pep) != L:
                continue
            if any(ch not in VALID_AA for ch in mut_pep):
                continue
            peps.append({
                "gene": rec["gene"], "hgvsp": rec["hgvsp"],
                "transcript": rec["transcript"],
                "length": L, "pep_start": start,
                "mut_peptide": mut_pep, "wt_peptide": wt_pep,
                "mut_offset": p - start + 1,
                "VAF": rec["VAF"],
                "chrom": rec["chrom"], "genomic_pos": rec["genomic_pos"],
            })

pdf = pd.DataFrame(peps).drop_duplicates(
    subset=["gene", "hgvsp", "mut_peptide"]).reset_index(drop=True)
uniq_peps = sorted(pdf.mut_peptide.unique())
print("      unique mutant peptides: %d" % len(uniq_peps), flush=True)
print("      by length: %s"
      % pdf.length.value_counts().sort_index().to_dict(), flush=True)

# mhcflurry rejects --alleles when given a file, so the allele goes in
# as a column: one row per peptide-allele pair
grid = pd.DataFrame([(p, a) for p in uniq_peps for a in alleles],
                    columns=["peptide", "allele"])
pep_file = "%s/mhcflurry_input.csv" % OUT_DIR
grid.to_csv(pep_file, index=False)
print("      input rows for mhcflurry: %d" % len(grid), flush=True)

print("[6/7] running mhcflurry ...", flush=True)
pred_file = "%s/mhcflurry_raw.csv" % OUT_DIR
cmd = ["mhcflurry-predict", pep_file,
       "--peptide-column", "peptide",
       "--allele-column", "allele",
       "--no-flanking",
       "--out", pred_file]
print("      %s" % " ".join(cmd), flush=True)
r = subprocess.run(cmd, capture_output=True, text=True)
if r.returncode != 0:
    print(r.stdout[-2000:])
    print(r.stderr[-3000:])
    raise SystemExit("mhcflurry failed")

pred = pd.read_csv(pred_file)
print("      predictions: %d rows" % len(pred), flush=True)
print("      columns: %s" % list(pred.columns), flush=True)

print("[7/7] merging and scoring ...", flush=True)
pcol = "mhcflurry_presentation_percentile"
if pcol not in pred.columns:
    pcol = "mhcflurry_affinity_percentile"
print("      ranking on: %s" % pcol, flush=True)

out = pdf.merge(pred, left_on="mut_peptide", right_on="peptide", how="inner")
out["binder"] = np.where(out[pcol] < STRONG_PCT, "STRONG",
                  np.where(out[pcol] < WEAK_PCT, "WEAK", "NONE"))
out = out.sort_values(pcol)

out.to_csv("%s/neoantigens_all.tsv" % OUT_DIR, sep="\t", index=False)
out[out.binder != "NONE"].to_csv(
    "%s/neoantigens_binders.tsv" % OUT_DIR, sep="\t", index=False)

per_mut = out.groupby(["gene", "hgvsp", "VAF"]).agg(
    peptides=("mut_peptide", "nunique"),
    best_percentile=(pcol, "min"),
    strong=("binder", lambda s: int((s == "STRONG").sum())),
    weak=("binder", lambda s: int((s == "WEAK").sum())),
).reset_index().sort_values("best_percentile")
per_mut.to_csv("%s/neoantigens_per_mutation.tsv" % OUT_DIR,
               sep="\t", index=False)

pd.set_option("display.width", 220)
print()
print("========== NEOANTIGEN PREDICTION ==========")
print("HLA alleles     : %s" % ", ".join(alleles))
print("missense inputs : %d of %d usable" % (len(records), len(mis)))
print("peptides tested : %d" % out.mut_peptide.nunique())
print("predictions     : %d peptide-allele pairs" % len(out))
print("strong binders  : %d  (percentile < %.1f)"
      % (int((out.binder == "STRONG").sum()), STRONG_PCT))
print("weak binders    : %d  (percentile < %.1f)"
      % (int((out.binder == "WEAK").sum()), WEAK_PCT))
n_yield = int((per_mut.strong + per_mut.weak > 0).sum())
print("mutations yielding at least one binder: %d / %d"
      % (n_yield, len(per_mut)))

print()
print("--- TOP 20 PEPTIDES ---")
cols = [c for c in ["gene", "hgvsp", "mut_peptide", "wt_peptide", "allele",
                    "mhcflurry_affinity", pcol, "binder", "VAF"]
        if c in out.columns]
print(out.head(20)[cols].to_string(index=False,
                                   float_format=lambda x: "%.3f" % x))

print()
print("--- BINDERS PER ALLELE ---")
b = out[out.binder != "NONE"]
if len(b):
    print(pd.crosstab(b.allele, b.binder).to_string())
else:
    print("  none")

print()
print("--- TOP 15 MUTATIONS BY BEST PEPTIDE ---")
print(per_mut.head(15).to_string(index=False,
                                 float_format=lambda x: "%.3f" % x))
print()
print("files written to %s" % OUT_DIR)
