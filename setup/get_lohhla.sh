#!/bin/bash
#
# Obtain LOHHLA and make it run on this stack.
#
# LOHHLA is not redistributed here. Its licence says plainly that it is
# provided as a personal copy and must not be forwarded to anyone else, so
# this repository carries only the changes — the patches below — and this
# script fetches the original from the authors' own GitHub. Running it
# means accepting their terms, which are in the README it clones.
#
# What the patches do, and why each was necessary:
#
#   novoalign is commercial and available neither through conda nor as a
#   cluster module. bwa mem -a replaces it; the -a flag keeps every
#   alignment, matching novoalign's "-r All 9999". That matters here
#   because a read often matches both alleles of a locus and LOHHLA makes
#   the distinction later, from mismatch positions.
#
#   GATK 3 jars — SortSam, SamToFastq, FilterSamReads — no longer exist in
#   GATK 4. samtools does all three.
#
#   The HLA coordinates in the original are hg19 without chr prefixes.
#   Against an hg38 BAM every extraction returns nothing.
#
#   The seven hg19 alternative haplotype contigs (6_apd_hap1 and the rest)
#   do not exist in hg38. Their hg38 equivalents are not substituted in:
#   reads there carry mapping quality zero and would only add noise.
#
#   LOHHLA declares short option flags like -id and -nBAM. Current
#   optparse rejects anything longer than one letter and the script aborts
#   while building its own option list, before reading any argument.
#
#   samtools fastq drops singletons on a coordinate-sorted file, and
#   LOHHLA extracts from three narrow windows where many mates fall
#   outside. On one sample 6444 of 6458 reads were lost. Sorting by name
#   first keeps the pairs together.
#
# None of the statistics are touched: per-allele coverage, log ratio,
# allelic imbalance and the confidence intervals are the authors' own.
#
# Usage:
#   ./get_lohhla.sh              clone and patch
#   ./get_lohhla.sh --force      discard an existing copy and start again
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

FORCE=0
[ "${1:-}" = "--force" ] && FORCE=1

REPO="https://github.com/mskcc/lohhla.git"
SCRIPT="${LOHHLA_DIR}/LOHHLAscript.R"

step_header "LOHHLA"

if [ -s "${SCRIPT}" ] && [ "${FORCE}" = "0" ]; then
    if grep -q "bwa mem -a" "${SCRIPT}" 2>/dev/null; then
        note "already present and patched"
        note "use --force to fetch a fresh copy"
        exit 0
    fi
    note "present but unpatched — patching in place"
else
    if [ "${FORCE}" = "1" ]; then
        note "removing the existing copy"
        rm -rf "${LOHHLA_DIR}/LOHHLAscript.R" "${LOHHLA_DIR}/.git" \
               "${LOHHLA_DIR}/README.md" "${LOHHLA_DIR}/LOHHLAscript.R.orig"
    fi

    note "cloning from ${REPO}"
    TMP="$(mktemp -d)"
    if ! git clone --depth 1 -q "${REPO}" "${TMP}/lohhla" 2>"${TMP}/err"; then
        fail "clone failed:"
        cat "${TMP}/err" | sed 's/^/     /'
        echo ""
        note "LOHHLA can also be fetched by hand from"
        note "  https://github.com/mskcc/lohhla"
        note "put LOHHLAscript.R at ${SCRIPT} and run this again"
        rm -rf "${TMP}"
        exit 1
    fi

    mkdir -p "${LOHHLA_DIR}"
    cp "${TMP}/lohhla/LOHHLAscript.R" "${SCRIPT}"
    [ -s "${TMP}/lohhla/README.md" ] && cp "${TMP}/lohhla/README.md" "${LOHHLA_DIR}/"
    rm -rf "${TMP}"
    note "cloned"
fi

# the original is kept so the diff can be inspected and so a failed patch
# run can start over without another clone
[ -s "${SCRIPT}.orig" ] || cp "${SCRIPT}" "${SCRIPT}.orig"
cp "${SCRIPT}.orig" "${SCRIPT}"

cd "${LOHHLA_DIR}"

say "patching"

# ------------------------------------------------------------------ 1
note "[1/6] novoindex -> bwa index"
sed -i "s|novoindexCMD <- paste('novoindex ', workDir, '/', full.patient, '.patient.hlaFasta.nix', ' ', workDir, '/', full.patient, '.patient.hlaFasta.fa', sep = '')|novoindexCMD <- paste('bwa index ', workDir, '/', full.patient, '.patient.hlaFasta.fa', sep = '')|" LOHHLAscript.R

# ------------------------------------------------------------------ 2
note "[2/6] novoalign -> bwa mem -a"
python3 - << 'PYEOF'
p = "LOHHLAscript.R"
s = open(p).read()
old_start = "    alignCMD <- paste(NOVODir, '/novoalign -d '"
i = s.find(old_start)
assert i >= 0, "novoalign line not found"
j = s.find("\n", i)
new = ("    alignCMD <- paste('bwa mem -a -t 4 ', workDir, '/', full.patient, "
       "'.patient.hlaFasta.fa ', regionDir, \"/\", BAMid, \".chr6region.1.fastq \", "
       "regionDir, \"/\", BAMid, \".chr6region.2.fastq\", ' 1> ', regionDir, '/', "
       "BAMid, '.chr6region.patient.reference.hlas.sam ', '2> ', regionDir, '/', "
       "BAMid, '_BS_GL.chr6region.patient.reference.hlas.metrics', sep = '')")
s = s[:i] + new + s[j:]
open(p, "w").write(s)
PYEOF

# ------------------------------------------------------------------ 3
note "[3/6] GATK 3 jars -> samtools"
python3 - << 'PYEOF'
p = "LOHHLAscript.R"
s = open(p).read()

# SortSam
old = 'sortBAM <- paste("java -jar ",GATKDir,"/SortSam.jar"'
i = s.find(old)
assert i >= 0, "SortSam line not found"
j = s.find("\n", i)
new = ('    sortBAM <- paste("samtools sort -o ", regionDir, "/", BAMid, '
       '".chr6region.patient.reference.hlas.csorted.bam ", regionDir, "/", BAMid, '
       '".chr6region.patient.reference.hlas.bam", sep="")')
s = s[:i - 4] + new + s[j:]

# SamToFastq inside get.partially.matching.reads
marker = "  cmd <- paste('java -jar ', GATKDir, '/SamToFastq.jar I='"
n = 0
while True:
    i = s.find(marker)
    if i < 0:
        break
    j = s.find("\n", i)
    new = ("  cmd <- paste('samtools fastq -1 ', regionDir, '/fished.1.fastq -2 ', "
           "regionDir, '/fished.2.fastq -0 /dev/null -s /dev/null -n ', "
           "regionDir, '/fished.sam', sep = '')")
    s = s[:i] + new + s[j:]
    n += 1

# the second SamToFastq, in the main loop, has a different argument layout
marker2 = 'samToFastQ <- paste("java -jar ",GATKDir,"/SamToFastq.jar "'
i = s.find(marker2)
assert i >= 0, "main SamToFastq not found"
line_start = s.rfind("\n", 0, i) + 1
line_end = s.find("\n", i)
indent = s[line_start:i]
new = (indent + 'samToFastQ <- paste("samtools fastq -1 ", regionDir, "/", BAMid, '
       '".chr6region.1.fastq -2 ", regionDir, "/", BAMid, '
       '".chr6region.2.fastq -0 /dev/null -s /dev/null -n ", '
       'regionDir, "/", BAMid, ".chr6region.sam", sep="")')
s = s[:line_start] + new + s[line_end:]

# FilterSamReads -> samtools view -N, which arrived in samtools 1.12
marker3 = 'extractCMD <- paste("java -jar ",GATKDir,"/FilterSamReads.jar "'
i = s.find(marker3)
assert i >= 0, "FilterSamReads line not found"
line_start = s.rfind("\n", 0, i) + 1
line_end = s.find("\n", i)
indent = s[line_start:i]
new = (indent + 'extractCMD <- paste("samtools view -b -N ", regionDir, "/", BAMid, '
       '\'.\', allele, ".passed.reads.txt -o ", regionDir, "/", BAMid, '
       '".type.", allele, ".filtered.bam ", regionDir, "/", BAMid, '
       '".type.", allele, ".bam", sep = "")')
s = s[:line_start] + new + s[line_end:]

open(p, "w").write(s)
PYEOF

# ------------------------------------------------------------------ 4
note "[4/6] hg19 coordinates -> hg38 with chr prefixes"
sed -i 's|" 6:29909037-29913661 >> "|" chr6:29941260-29949572 >> "|' LOHHLAscript.R
sed -i 's|" 6:31321649-31324964 >> "|" chr6:31352875-31368067 >> "|' LOHHLAscript.R
sed -i 's|" 6:31236526-31239869 >> "|" chr6:31267749-31273130 >> "|' LOHHLAscript.R

python3 - << 'PYEOF'
p = "LOHHLAscript.R"
lines = open(p).read().split("\n")
haps = ["6_apd_hap1", "6_cox_hap2", "6_dbb_hap3", "6_mann_hap4",
        "6_mcf_hap5", "6_qbl_hap6", "6_ssto_hap7"]
n = 0
for k, ln in enumerate(lines):
    if any(h in ln for h in haps) and ln.strip().startswith("samtoolsCMD"):
        indent = ln[:len(ln) - len(ln.lstrip())]
        lines[k] = indent + "## hg19 haplotype contig, absent in hg38: " + ln.strip()
        for off in (1, 2):
            if k + off < len(lines):
                st = lines[k + off].strip()
                if st.startswith("write.table(samtoolsCMD") or \
                   st.startswith("system(samtoolsCMD"):
                    ind2 = lines[k + off][:len(lines[k + off]) -
                                          len(lines[k + off].lstrip())]
                    lines[k + off] = ind2 + "## " + st
        n += 1
open(p, "w").write("\n".join(lines))
print("        %d contig queries disabled" % n)
PYEOF

# ------------------------------------------------------------------ 5
note "[5/6] removing invalid short option flags"
python3 - << 'PYEOF'
import re
p = "LOHHLAscript.R"
s = open(p).read()
pattern = re.compile(r'make_option\(c\("(-[A-Za-z]+)",\s*("--[A-Za-z]+")\)')
matches = pattern.findall(s)
s = pattern.sub(lambda m: 'make_option(c(%s)' % m.group(2), s)
open(p, "w").write(s)
print("        %d short flags removed" % len(matches))
PYEOF

# ------------------------------------------------------------------ 6
note "[6/6] keeping read pairs through FASTQ conversion"
python3 - << 'PYEOF'
p = "LOHHLAscript.R"
s = open(p).read()
marker = 'samToFastQ <- paste("samtools fastq -1 "'
i = s.find(marker)
assert i >= 0, "main samToFastQ line not found"
line_start = s.rfind("\n", 0, i) + 1
line_end = s.find("\n", i)
indent = s[line_start:i]
new = (indent + 'samToFastQ <- paste("samtools sort -n -o ", regionDir, "/", BAMid, '
       '".chr6region.nsorted.bam ", regionDir, "/", BAMid, '
       '".chr6region.sam && samtools fastq -1 ", regionDir, "/", BAMid, '
       '".chr6region.1.fastq -2 ", regionDir, "/", BAMid, '
       '".chr6region.2.fastq -s ", regionDir, "/", BAMid, '
       '".chr6region.single.fastq -0 /dev/null -n ", '
       'regionDir, "/", BAMid, ".chr6region.nsorted.bam", sep="")')
s = s[:line_start] + new + s[line_end:]
open(p, "w").write(s)
PYEOF

say "verification"

# the word survives in comments, in the help text for an option that is
# now ignored, and in the name of a variable that runs bwa index; what
# matters is that nothing invokes the binary
printf "  %-34s" "no novoalign invoked"
if grep -nE "system\(.*novoalign|paste\(NOVODir" LOHHLAscript.R \
     | grep -qv "^[0-9]*: *#"; then
    echo "FAIL"
else
    echo "ok"
fi

printf "  %-34s" "no GATK 3 jars left"
grep -n "\.jar" LOHHLAscript.R | grep -qv "^[0-9]*: *##" && echo "FAIL" || echo "ok"

printf "  %-34s" "hg38 coordinates present"
grep -q "chr6:299\|chr6:313\|chr6:312" LOHHLAscript.R && echo "ok" || echo "FAIL"

printf "  %-34s" "no multi-letter short flags"
grep -q 'make_option(c("-[A-Za-z][A-Za-z]' LOHHLAscript.R && echo "FAIL" || echo "ok"

printf "  %-34s" "name-sort before fastq"
grep -q "samtools sort -n" LOHHLAscript.R && echo "ok" || echo "FAIL"

printf "  %-34s" "R parses the result"
if use_env "${ENV_LOHHLA}" 2>/dev/null && command -v Rscript >/dev/null; then
    n=$(Rscript -e 'cat(length(parse("LOHHLAscript.R")))' 2>/dev/null)
    [ -n "${n}" ] && echo "ok, ${n} expressions" || echo "FAIL"
else
    echo "skipped, no R"
fi

diff -u LOHHLAscript.R.orig LOHHLAscript.R > lohhla.patch 2>/dev/null
echo ""
note "diff against the original written to ${LOHHLA_DIR}/lohhla.patch"
note "$(wc -l < lohhla.patch) lines changed"
