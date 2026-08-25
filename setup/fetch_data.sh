#!/bin/bash
#
# Everything the demo needs that does not ship with the code.
#
# Three kinds of thing are missing from the repository and for different
# reasons. The reference genome and the VEP cache are too large to
# version. The IMGT/HLA sequence set is CC BY-NoDerivs, which permits
# copying the data as published but not redistributing a reformatted
# version, so it is fetched and rebuilt here rather than shipped.
# netMHCpan is licensed per user by DTU and cannot be redistributed at
# all; this script only says where to get it.
#
# Roughly 32 GB is downloaded. Everything lands under DATA_DIR from
# config.sh, never inside the repository. Each item is skipped if it is
# already present and complete, so an interrupted run can simply be
# started again.
#
# Usage:
#   ./fetch_data.sh            everything
#   ./fetch_data.sh --no-vep   skip the 26 GB cache, which only step 6 needs
#   ./fetch_data.sh --list     say what would be fetched and stop
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

SKIP_VEP=0
LIST_ONLY=0
for a in "$@"; do
    case "$a" in
        --no-vep) SKIP_VEP=1 ;;
        --list)   LIST_ONLY=1 ;;
        *) echo "unknown option: $a" >&2; exit 1 ;;
    esac
done

# the old genomics-public-data bucket now refuses anonymous reads; the
# same files are served from gcp-public-data--broad-references, which is
# what GATK's own documentation points at
GATK_BUCKET="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
VEP_URL="https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz"
IMGT_RELEASE="${IMGT_RELEASE:-3610}"
IMGT_URL="https://raw.githubusercontent.com/ANHIG/IMGTHLA/${IMGT_RELEASE}/hla_gen.fasta.zip"
IMGT_DAT_URL="https://raw.githubusercontent.com/ANHIG/IMGTHLA/${IMGT_RELEASE}/hla.dat.zip"
HLA_TYPES_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt"
GENCODE_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46"
ZENODO_API="https://zenodo.org/api/records/${ZENODO_RECORD}"

say "where things will go"
note "repository   ${DEMO_ROOT}"
note "data         ${DATA_DIR}"
note "reference    ${REF_DIR}"
note "samples      ${SAMPLE_DIR}"
note "results      ${OUT_DIR}"

if [ "${LIST_ONLY}" = "1" ]; then
    say "what would be fetched"
    note "hg38 reference and index        ~3.2 GB   ${GATK_BUCKET}"
    note "1000 Genomes HLA types          ~200 KB   EBI"
    note "IMGT/HLA sequences              ~140 MB   ANHIG/IMGTHLA"
    [ "${SKIP_VEP}" = "0" ] && \
    note "VEP cache, GRCh38               ~26 GB    Ensembl"
    note "demo samples                    ~2 GB     ${ZENODO_DOI}"
    note ""
    note "netMHCpan is not downloadable: request a licence at"
    note "  https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/"
    note "and unpack it to ${NETMHC_DIR}"
    exit 0
fi

mkdir -p "${REF_DIR}" "${SAMPLE_DIR}" "${OUT_DIR}"

FAILED=0

# a download is worth repeating only if the previous attempt left nothing
# usable, so size is checked rather than mere existence
have() {
    local path="$1" min="${2:-1000}"
    [ -s "${path}" ] && [ "$(stat -Lc%s "${path}" 2>/dev/null || echo 0)" -ge "${min}" ]
}

get() {
    local url="$1" dest="$2" min="${3:-1000}" label="$4"
    if have "${dest}" "${min}"; then
        note "${label}: already here"
        return 0
    fi
    note "${label}: downloading"
    if curl -fL --retry 3 --retry-delay 5 -C - --progress-bar \
            -o "${dest}.part" "${url}"; then
        mv "${dest}.part" "${dest}"
        note "${label}: $(du -Lh "${dest}" | cut -f1)"
        return 0
    fi
    rm -f "${dest}.part"
    fail "${label}: download failed"
    FAILED=$((FAILED+1))
    return 1
}

# ---------------------------------------------------------------- hg38
say "reference genome"
get "${GATK_BUCKET}/Homo_sapiens_assembly38.fasta" \
    "${REF_FASTA}" 3000000000 "hg38 fasta"
get "${GATK_BUCKET}/Homo_sapiens_assembly38.fasta.fai" \
    "${REF_FASTA}.fai" 100000 "hg38 index"
get "${GATK_BUCKET}/Homo_sapiens_assembly38.dict" \
    "${REF_FASTA%.fasta}.dict" 50000 "hg38 dictionary"

# GATK wants a .dict beside the fasta and named for it; the bucket's copy
# is named for the fasta without its extension, which is what GATK expects
if have "${REF_FASTA}" 3000000000 && ! have "${REF_FASTA}.bwt" 100000; then
    note "bwa index: building, this takes about an hour"
    if use_env "${ENV_BAMSURGEON}" 2>/dev/null && command -v bwa >/dev/null; then
        bwa index "${REF_FASTA}" 2>"${REF_DIR}/bwa_index.log" \
            && note "bwa index: done" \
            || { fail "bwa index failed, see ${REF_DIR}/bwa_index.log"
                 FAILED=$((FAILED+1)); }
    else
        fail "bwa not available; run setup/install.sh first, then rerun"
        FAILED=$((FAILED+1))
    fi
else
    have "${REF_FASTA}.bwt" 100000 && note "bwa index: already here"
fi

# ------------------------------------------------------------ HLA types
say "published HLA types"
note "the 1000 Genomes panel typed 2693 individuals by Sanger sequencing;"
note "the demo compares OptiType against them"
get "${HLA_TYPES_URL}" "${HLA_TYPES}" 100000 "HLA types"

# ------------------------------------------------------------ IMGT/HLA
say "IMGT/HLA sequences"
note "CC BY-NoDerivs: fetched as published, reformatted locally, never"
note "redistributed. LOHHLA needs one sequence per allele name."
IMGT_RAW="${REF_DIR}/hla_gen.fasta"
IMGT_ZIP="${REF_DIR}/hla_gen.fasta.zip"
if have "${IMGT_RAW}" 100000000; then
    note "IMGT hla_gen.fasta: already unpacked"
elif get "${IMGT_URL}" "${IMGT_ZIP}" 20000000 "IMGT archive"; then
    note "unpacking"
    unzip -o -q -d "${REF_DIR}" "${IMGT_ZIP}" && rm -f "${IMGT_ZIP}"
fi
if have "${IMGT_RAW}" 100000000; then
    # a fresh clone has no soft/lohhla/data: the directory is created by
    # whichever step writes into it first, and that is this one
    mkdir -p "${LOHHLA_DIR}/data"
    LOHHLA_FA="${LOHHLA_DIR}/data/hla_all_lohhla.fasta"
    if have "${LOHHLA_FA}" 1000000; then
        note "LOHHLA reference: already built"
    else
        note "LOHHLA reference: rebuilding from the IMGT set"
        python3 - "${IMGT_RAW}" "${LOHHLA_FA}" << 'PY'
import sys, re

src, dst = sys.argv[1], sys.argv[2]

# LOHHLA identifies an allele by a name of the form hla_a_02_01 and uses
# exactly one sequence per name regardless of how many subtypes the file
# holds, so only the first occurrence of each two-field name is kept.
seen, out, keep = set(), [], False
name = None
for line in open(src):
    if line.startswith(">"):
        m = re.search(r"(HLA:\S+\s+)?([A-Z0-9]+)\*(\d+):(\d+)", line)
        keep = False
        if m:
            gene, f1, f2 = m.group(2), m.group(3), m.group(4)
            if gene in ("A", "B", "C"):
                name = f"hla_{gene.lower()}_{f1}_{f2}"
                if name not in seen:
                    seen.add(name)
                    keep = True
                    out.append(f">{name}\n")
    elif keep:
        out.append(line)

with open(dst, "w") as fh:
    fh.writelines(out)
print(f"    {len(seen)} allele names written")
PY
        have "${LOHHLA_FA}" 1000000 \
            && note "LOHHLA reference: $(du -h "${LOHHLA_FA}" | cut -f1)" \
            || { fail "rebuild produced nothing usable"; FAILED=$((FAILED+1)); }
    fi
fi

# ------------------------------------------------------------- gencode
say "GENCODE annotation"
note "the binding predictor needs the coding sequence a mutation falls in"
note "and the protein it translates to, neither of which is in the BAM"
get "${GENCODE_BASE}/gencode.v46.annotation.gtf.gz" \
    "${REF_DIR}/gencode.v46.annotation.gtf.gz" 40000000 "GTF"
get "${GENCODE_BASE}/gencode.v46.pc_translations.fa.gz" \
    "${REF_DIR}/gencode.v46.pc_translations.fa.gz" 8000000 "protein FASTA"

say "IMGT exon annotation"
note "LOHHLA reads hla.dat to know where the exons sit; the FASTA alone"
note "gives it sequence but no structure"
mkdir -p "${LOHHLA_DIR}/data"
DAT_OUT="${LOHHLA_DIR}/data/hla.dat"
DAT_ZIP="${REF_DIR}/hla.dat.zip"
if have "${DAT_OUT}" 50000000; then
    note "hla.dat: already here"
elif get "${IMGT_DAT_URL}" "${DAT_ZIP}" 20000000 "hla.dat archive"; then
    note "unpacking"
    unzip -o -q -d "${LOHHLA_DIR}/data" "${DAT_ZIP}" && rm -f "${DAT_ZIP}"
    have "${DAT_OUT}" 50000000 \
        && note "hla.dat: $(du -Lh "${DAT_OUT}" | cut -f1)" \
        || fail "hla.dat did not unpack"
fi

# --------------------------------------------------------------- VEP
if [ "${SKIP_VEP}" = "1" ]; then
    say "VEP cache"
    note "skipped; step 6 will not run without it"
else
    say "VEP cache"
    note "26 GB, needed only by the pVACseq step"
    if [ -d "${VEP_CACHE}/homo_sapiens" ]; then
        note "already here: $(du -sh "${VEP_CACHE}" | cut -f1)"
    else
        mkdir -p "${VEP_CACHE}"
        TARBALL="${VEP_CACHE}/vep_cache.tar.gz"
        if get "${VEP_URL}" "${TARBALL}" 10000000000 "VEP cache"; then
            note "unpacking"
            tar -xzf "${TARBALL}" -C "${VEP_CACHE}" \
                && rm -f "${TARBALL}" \
                && note "unpacked: $(du -sh "${VEP_CACHE}" | cut -f1)" \
                || { fail "unpacking failed"; FAILED=$((FAILED+1)); }
        fi
    fi
fi

# ------------------------------------------------------------- samples
say "demo samples"
note "from ${ZENODO_DOI}"
if command -v python3 >/dev/null; then
    python3 - "${ZENODO_API}" "${SAMPLE_DIR}" << 'PY'
import json, os, sys, urllib.request

api, dest = sys.argv[1], sys.argv[2]
try:
    with urllib.request.urlopen(api, timeout=30) as r:
        rec = json.load(r)
except Exception as e:
    print(f"  ERROR: could not reach Zenodo: {e}")
    print(f"  the record may still be a draft; publish it, or download the")
    print(f"  files by hand into {dest}")
    sys.exit(0)

files = rec.get("files", [])
if not files:
    print(f"  the record has no files yet")
    sys.exit(0)

for f in files:
    name = f.get("key")
    url = f.get("links", {}).get("self")
    size = f.get("size", 0)
    out = os.path.join(dest, name)
    if os.path.exists(out) and os.path.getsize(out) == size:
        print(f"  {name}: already here")
        continue
    print(f"  {name}: downloading {size/1e6:.0f} MB")
    urllib.request.urlretrieve(url, out)
PY
fi

# ------------------------------------------------------------ netMHCpan
say "netMHCpan"
if [ -x "${NETMHC_DIR}/netMHCpan" ]; then
    note "found at ${NETMHC_DIR}"
else
    note "not present, and it cannot be downloaded automatically."
    note ""
    note "  request a licence:"
    note "    https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/"
    note "  unpack the archive to:"
    note "    ${NETMHC_DIR}"
    note "  edit netMHCpan and set NMHOME to that path"
    note ""
    note "the pVACseq step is skipped without it; everything else runs"
fi

say "summary"
for item in "${REF_FASTA}:hg38" "${REF_FASTA}.fai:index" \
            "${REF_DIR}/gencode.v46.annotation.gtf.gz:GENCODE GTF" \
            "${REF_DIR}/gencode.v46.pc_translations.fa.gz:GENCODE proteins" \
            "${HLA_TYPES}:HLA types" \
            "${LOHHLA_DIR}/data/hla_all_lohhla.fasta:LOHHLA reference" \
            "${LOHHLA_DIR}/data/hla.dat:IMGT exon annotation"; do
    p="${item%%:*}"; l="${item##*:}"
    have "${p}" 1000 && note "ok       ${l}" || note "missing  ${l}"
done
[ -d "${VEP_CACHE}/homo_sapiens" ] && note "ok       VEP cache" \
                                   || note "missing  VEP cache"
[ -x "${NETMHC_DIR}/netMHCpan" ] && note "ok       netMHCpan" \
                                 || note "missing  netMHCpan (optional)"
n=$(ls "${SAMPLE_DIR}"/*.bam 2>/dev/null | wc -l)
note "${n} sample BAMs in ${SAMPLE_DIR}"

echo ""
if [ "${FAILED}" -gt 0 ]; then
    echo "  ${FAILED} item(s) failed; rerun to retry only those"
    exit 1
fi
echo "  ready — run setup/check_deps.sh next"
