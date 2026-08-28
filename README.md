# Immune Escape Convergence — demo

A five-sample run of a pipeline that places known mutations into real
sequencing reads, then asks standard tools to find them again. Because
every mutation was put there deliberately, what the tools recover can be
checked against what is actually in the file rather than against a
consensus of other tools.

Each sample pairs somatic mutations from a TCGA breast or ovarian tumour
with panel reads from a healthy 1000 Genomes individual. Neither half is
simulated. What is synthetic is the pairing.

The pipeline runs eight stages:

```
1a  place substitutions          BAMSurgeon addsnv
1b  place indels                 BAMSurgeon addindel
2   call variants                Mutect2
3   type HLA                     OptiType
4   simulate loss of one locus   samtools
5   test for that loss           LOHHLA
6   predict binding              mhcflurry
7   predict binding again        pVACseq with NetMHCpan
8   report                       what was placed against what came back
```

---

## What you need

**Disk** — about 40 GB. Reference data is 30 GB of it and is downloaded
once; the run itself produces around 8 GB.

**Time** — roughly 45 minutes of compute per sample. On a cluster the
five run in parallel and the whole thing takes under an hour; sequentially
it is closer to four.

**conda** — seven environments are built from `setup/environments/*.yml`.
They pin the versions this was run against.

**Two things that cannot be shipped:**

*netMHCpan* is licensed per user by DTU. Request one at
<https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/>, unpack the
archive to `soft/netMHCpan-4.1`, and edit its `netMHCpan` script to point
`NMHOME` at that directory. Without it stage 7 is skipped and everything
else still runs — stage 6 already predicts binding with mhcflurry.

*LOHHLA* is provided by its authors as a personal copy and its licence
forbids passing it on. `setup/get_lohhla.sh` fetches the original from
their GitHub and applies `soft/lohhla/lohhla.patch`, which is ours. What
the patch changes and why is at the top of that script.

---

## Running it

```bash
git clone https://github.com/selevycho/immune-escape-demo.git
cd immune-escape-demo

# where data and results go — pick a filesystem with room
export DATA_DIR=/path/with/40GB

./setup/install.sh          # seven conda environments, about 20 minutes
./setup/get_lohhla.sh       # clone and patch LOHHLA
./setup/fetch_data.sh       # 30 GB, an hour or two
./setup/check_deps.sh       # says what will run and what will be skipped

./run_demo.sh
```

`check_deps.sh` is worth reading before starting. It reports every
requirement and, at the end, which stages will run on this machine and
which will be skipped for want of something optional.

`run_demo.sh` submits to SLURM if `sbatch` is available, chaining the
stages so none begins before its input exists. Without SLURM, or with
`--local`, they run in sequence in the current shell.

Every stage skips work already done, so an interrupted run can simply be
started again:

```bash
./run_demo.sh --from 3      # start at variant calling
./run_demo.sh --only 5      # one stage
./run_demo.sh --dry-run     # say what would happen and stop
```

---

## What you should see

The report at the end compares what was placed with what was recovered.
On this data it says roughly:

```
mutations placed         197
recovered by Mutect2     157   (80%)
HLA genotypes            5 distinct of 5
HLA loss detected        2 loci, both the one that was thinned
strong binders           218
```

Three things in that report are worth more than the headline number.

**Recall follows the allele fraction, not the sample.** Nothing below
about three supporting reads is recoverable — at 34× coverage a mutation
present in 5% of cells sits in under two reads. The per-fraction table
makes that floor visible:

```
5-10%      0%
10-15%    50%
15-20%    71%
20-30%    88%
>30%      97%
```

**The MHC is invisible to the caller.** Not one substitution placed there
is recovered. hg38 carries alternate haplotypes for that region, reads map
ambiguously between them, mapping quality collapses, and Mutect2 declines
to open an active region at all. This constrains the project's own
subject: immune escape happens in the MHC, and standard somatic calling
cannot see mutations there.

**Misses divide into two kinds.** Some variants were called and then
rejected by a filter; others were never called at all. The first is a
threshold and can be changed, the second is a detection limit and cannot.
The report separates them.

---

## What is simulated, and what is not

The reads are real. The mutations are real, taken from TCGA cases with the
allele fraction each carried in its donor tumour. The pairing between them
is the synthetic part.

Two limits are worth stating plainly.

**The HLA loss is a loss of coverage, not of an allele.** Thinning drops
reads across a whole locus, which reduces depth without creating allelic
imbalance — and allelic imbalance is what LOHHLA tests for. Separating the
two alleles of an HLA gene needs more discriminating positions than panel
data at this depth provides; an allele-specific version was attempted and
split reads 152 against 20 where an even split was needed. So stage 5
measures sensitivity to depth, and the report says so.

**The HLA genotype belongs to the backbone, not to the donor.** Peptide
predictions are what these mutations would present in someone carrying
this genotype. That is the right quantity for measuring what a pipeline
detects. It is not a claim about the TCGA patient's immune landscape.

---

## Layout

```
config.sh              every path, in one file — the whole configuration
run_demo.sh            runs the stages in order

setup/
  install.sh           builds the conda environments
  get_lohhla.sh        clones LOHHLA and applies our patch
  fetch_data.sh        downloads reference data and the samples
  check_deps.sh        reports what is present and what will run
  clean.sh             removes scratch, or results, on request
  environments/        one yml per environment

steps/                 one script per stage
lib/                   the two python programs the stages call
panel/                 the 350-gene panel, 23.81 Mb over 5 373 intervals
soft/lohhla/           our patch and notes; the tool itself is fetched
```

Nothing is written inside the repository. `DATA_DIR` holds the reference
data, the samples and every result.

---

## Data

All forty pairs, with the analysis that produced the published numbers,
are at <https://doi.org/10.5281/zenodo.22135219> and
<https://github.com/selevycho/immune_escape_cohort>.

The five samples used here come from
<https://doi.org/10.5281/zenodo.22086838>:
five normal BAM slices and, for each, the mutations to be placed with
their hg38 coordinates and allele fractions. `fetch_data.sh` downloads
them. The tumour BAMs are not included — the pipeline builds them, which
is the point.

Reference data is fetched from its original sources: hg38 from the Broad
public bucket, GENCODE 46 from EBI, IMGT/HLA release 3610 from the ANHIG
repository, the VEP cache from Ensembl, and the published 1000 Genomes HLA
types from EBI. The IMGT release is pinned rather than tracking `Latest`,
both so the download keeps working and so the allele set is the one these
results came from.

---

## Licence

The pipeline code here is MIT.

The tools it calls are not: LOHHLA is under the Francis Crick Institute's
academic terms, netMHCpan under DTU's, IMGT/HLA under CC BY-NoDerivs. None
of the three is redistributed. `fetch_data.sh` and `get_lohhla.sh` obtain
them from their sources, which is where their terms are stated.
