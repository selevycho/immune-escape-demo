#!/bin/bash
#
# Assemble the files that go to Zenodo.
#
# The demo needs input the repository cannot carry: five normal BAM slices
# and the mutations to place into them. Those come from a DOI rather than
# from git, both because of size and because data and code have different
# lifecycles — the code will change, the input should not.
#
# What is deliberately not included: the tumour BAMs. The demo builds them
# by injection, and shipping them would let someone verify the analysis
# without ever running the pipeline that produced them, which is the
# opposite of the point.
#
# Usage:
#   ./make_zenodo_package.sh              build under DATA_DIR/zenodo
#   ./make_zenodo_package.sh /some/path   build somewhere else
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

SRC_WS="${SRC_WS:-$(ws_find immune_escape 2>/dev/null)}"
PKG="${1:-${DATA_DIR}/zenodo}"

# the shell running this may be a bare conda base without pandas, so a
# python that can actually read the tables is located rather than assumed
PY=""
for cand in python3 \
    "$(conda info --base 2>/dev/null)/envs/${ENV_GATK}/bin/python" \
    "$(conda info --base 2>/dev/null)/envs/${ENV_MHCFLURRY}/bin/python" \
    "$(conda info --base 2>/dev/null)/envs/cptac_env/bin/python"; do
    if "${cand}" -c "import pandas" >/dev/null 2>&1; then
        PY="${cand}"
        break
    fi
done
if [ -z "${PY}" ]; then
    fail "no python with pandas found; activate an environment that has it"
    exit 1
fi
note "using ${PY}"

if [ -z "${SRC_WS}" ] || [ ! -d "${SRC_WS}" ]; then
    fail "source workspace not found; set SRC_WS"
    exit 1
fi

COHORT="${SRC_WS}/simulation/cohort"
INDELS="${SRC_WS}/simulation/indels"

say "assembling from ${SRC_WS}"
note "into ${PKG}"
mkdir -p "${PKG}"

MISSING=0

for s in ${SAMPLES}; do
    note ""
    note "${s}"

    nb="${COHORT}/${s}/${s}_normal.bam"
    if [ -s "${nb}" ]; then
        cp -L "${nb}" "${PKG}/${s}_normal.bam"
        cp -L "${nb}.bai" "${PKG}/${s}_normal.bam.bai" 2>/dev/null
        note "  normal BAM   $(du -h "${PKG}/${s}_normal.bam" | cut -f1)"
    else
        fail "  normal BAM missing"
        MISSING=$((MISSING+1))
        continue
    fi

    # the truth set carries both kinds of mutation; the indel table is
    # separate upstream because the two were injected in different rounds,
    # but the demo places them together and wants one file
    ts="${COHORT}/${s}/truth_set.tsv"
    it="${INDELS}/indel_truth_all.tsv"
    if [ -s "${ts}" ]; then
        "${PY}" - "${ts}" "${it}" "${s}" "${PKG}/${s}_truth.tsv" << 'PYEND'
import sys, os
import pandas as pd

ts, it, sample, out = sys.argv[1:5]
t = pd.read_csv(ts, sep="\t")

# normalise to the columns the demo reads, whichever names upstream used
ren = {}
for a, b in [("Chromosome_hg38", "chrom"), ("Start_Position_hg38", "pos"),
             ("End_Position_hg38", "end"),
             ("Reference_Allele", "ref"), ("Tumor_Seq_Allele2", "alt"),
             ("Hugo_Symbol", "gene"), ("VAF", "vaf"),
             ("Variant_Type", "type"),
             ("Variant_Classification", "consequence"),
             ("HGVSp_Short", "hgvsp"),
             ("t_ref_count", "t_ref"), ("t_alt_count", "t_alt"),
             ("source_patient", "donor")]:
    if a in t.columns:
        ren[a] = b
t = t.rename(columns=ren)

# HGVSp_Short is what the neoantigen step reads to know which amino acid
# changes; without it a missense mutation cannot be turned into a peptide.
# The read counts and the donor barcode come along because they cost
# nothing and someone will eventually want to know where a mutation was
# seen before it was placed.
keep = [c for c in ["chrom", "pos", "end", "ref", "alt", "gene", "vaf",
                    "type", "consequence", "hgvsp", "t_ref", "t_alt",
                    "donor"] if c in t.columns]
t = t[keep]
t.insert(0, "sample", sample)

if os.path.exists(it):
    i = pd.read_csv(it, sep="\t")
    i = i[i["sample"] == sample] if "sample" in i.columns else i.iloc[0:0]
    if len(i):
        ren2 = {c: c for c in i.columns}
        for a, b in [("Chromosome_hg38", "chrom"),
                     ("Start_Position_hg38", "pos"),
                     ("Hugo_Symbol", "gene"), ("VAF", "vaf")]:
            if a in i.columns:
                ren2[a] = b
        i = i.rename(columns=ren2)
        cols = [c for c in t.columns if c in i.columns]
        # an indel row that duplicates a substitution row would be
        # injected twice, so anything already present is dropped
        key = ["chrom", "pos"]
        if all(k in i.columns for k in key):
            have = set(zip(t.chrom, t.pos))
            i = i[~i.apply(lambda r: (r.chrom, r.pos) in have, axis=1)]
        t = pd.concat([t, i[cols]], ignore_index=True)

t = t.sort_values(["chrom", "pos"]).reset_index(drop=True)
t.to_csv(out, sep="\t", index=False)

n_snv = int((t.get("type", pd.Series(dtype=str)) == "SNP").sum()) \
    if "type" in t.columns else len(t)
print(f"  truth set    {len(t)} rows, {n_snv} substitutions")
PYEND
    else
        fail "  truth set missing"
        MISSING=$((MISSING+1))
    fi
done

say "manifest"
"${PY}" - "${PKG}" "${SAMPLES}" "${SRC_WS}" << 'PYEND'
import sys, os, glob
import pandas as pd

pkg, samples, ws = sys.argv[1], sys.argv[2].split(), sys.argv[3]

man_src = f"{ws}/simulation/cohort/manifest.tsv"
src = pd.read_csv(man_src, sep="\t") if os.path.exists(man_src) else None

rows = []
for s in samples:
    t = f"{pkg}/{s}_truth.tsv"
    b = f"{pkg}/{s}_normal.bam"
    if not (os.path.exists(t) and os.path.exists(b)):
        continue
    d = pd.read_csv(t, sep="\t")
    r = {"sample": s,
         "normal_bam": f"{s}_normal.bam",
         "truth_set": f"{s}_truth.tsv",
         "mutations": len(d),
         "bam_MB": round(os.path.getsize(b) / 1e6, 1)}
    if "type" in d.columns:
        r["substitutions"] = int((d.type == "SNP").sum())
        r["indels"] = int((d.type != "SNP").sum())
    if src is not None:
        h = src[src.sample_id == s]
        if len(h):
            for c in ("cohort", "backbone", "population", "superpopulation"):
                if c in h.columns:
                    r[c] = h.iloc[0][c]
    rows.append(r)

m = pd.DataFrame(rows)
m.to_csv(f"{pkg}/manifest.tsv", sep="\t", index=False)
print(m.to_string(index=False))
PYEND

say "checksums"
( cd "${PKG}" && md5sum ./*.bam ./*.tsv ./*.bai 2>/dev/null > md5sums.txt )
note "$(wc -l < "${PKG}/md5sums.txt") files"

say "readme for the record"
cat > "${PKG}/README.txt" << 'EOF'
Immune Escape Convergence — demo input data
===========================================

Five normal BAM slices and the mutations to be placed into them.

Each sample pairs somatic mutations from a TCGA breast or ovarian tumour
with sequencing reads from a healthy 1000 Genomes individual. The reads
here are the healthy ones, untouched: the pipeline injects the mutations
itself, so that what it later recovers can be checked against what it
placed.

Files
-----
  <sample>_normal.bam        panel slice, 350 genes over 23.81 Mb
  <sample>_normal.bam.bai    index
  <sample>_truth.tsv         mutations to place: position on hg38, the
                             alternate base, and the allele fraction the
                             mutation carried in the donor tumour
  manifest.tsv               the five samples with their backbones
  md5sums.txt                checksums

The tumour BAMs are not included. The pipeline builds them.

Running
-------
  https://github.com/<user>/immune-escape-demo

  git clone ...
  ./setup/install.sh
  ./setup/fetch_data.sh
  ./setup/check_deps.sh
  ./run_demo.sh

Licence
-------
CC BY 4.0. The underlying reads are from the 1000 Genomes Project, which
places no restriction on redistribution; TCGA mutation coordinates are
from the open-access tier.
EOF
note "written"

say "package"
du -sh "${PKG}"
ls -la "${PKG}" | head -20

echo ""
if [ "${MISSING}" -gt 0 ]; then
    echo "  ${MISSING} item(s) missing — the package is incomplete"
    exit 1
fi
echo "  ready to upload to ${ZENODO_DOI}"
echo ""
echo "  size: $(du -sh "${PKG}" | cut -f1)"
echo "  upload through the Zenodo web form, then publish the record"
