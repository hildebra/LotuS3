# Coarse dereplication

Add `-coarseDerep X` to enable SDM coarse processing at a fractional identity
from **0.95 through 1.0**, inclusive. For example:

```sh
./lotus3 -i Example/ -m Example/miSeqMap.sm.txt -o coarse_run \
  -CL vsearch -coarseDerep 0.975 -mergePreClusterReads 0
```

This uses a 97.5% substitution-identity threshold for internal coarse search and
compression. **The main clustering FASTA and sample map still describe exact
effective R1 dereplicates.** The search/truncation cap continues to come from the
SDM options file. Differences in R2 alone do not split an R1 parent, but remain
separate paired seed candidates. Omitting `-coarseDerep` preserves the existing
workflow; explicit `1.0` still enables quality retention and the variant reader.

## Required SDM interface

The updated bundled SDM **3.53 beta** supports this interface. Earlier development
builds can share version labels but use incompatible output layouts. LotuS checks
`-help_flags` for `-seedSubclusters`, `-derepCoarseClusters`, and the statement
`Standard HQ output: .1.hq.fq; paired additionally .2.hq.fq` when coarse processing
or retained qualities are requested. The same configured binary is used for
preprocessing and seed extension.

## Files and counts

All modes use the same paths under `tmpFiles` (or the selected temporary directory):

| File | Use |
| --- | --- |
| `derep.fas` | Passing exact R1 representatives for clustering/backmapping |
| `derep.fas.rest` | Below-cutoff parents for additional-count recovery |
| `derep.map` | Parent sample counts, including below-cutoff parents |
| `derep.1.hq.fq` | R1 seed candidates, also for single-end input |
| `derep.2.hq.fq` | Synchronized R2 seed candidates for paired input |

With retained qualities, HQ files contain every retained exact sequence or full
pair variant once, with its own abundance and selected observed quality vectors.
They include each representative's own variant and below-cutoff parents. There
is no second parent-only record and no `.subclusters*.fq` export. SDM resolves
`parent.subN` names during seed extension and counts each parent's sample vector
once while evaluating its variants. Paired fragments are counted once.

## Effective options

`-coarseDerep` passes the percentage as `-derepIdentity`, enables
`-derepStoreQuals 1`, and explicitly disables binary differences, optional
subcluster FASTA, reassignment, and coarse-parent output. These correspond to
`-derepStoreDiffs 0 -derepSubclusterFasta 0 -derepReassign 0 -derepCoarseClusters 0`.
The normal main FASTA output is kept; DADA2 retains its existing main FASTQ mode.

Custom SDM options files can enable `derepStoreQuals` without `-coarseDerep`,
including at identity 100. LotuS selects `-seedSubclusters 1` from that effective
quality-retention setting, not from identity or filenames. Without the LotuS
coarse flag, explicitly configured reassignment/coarse-parent output and
optional exports remain controlled by the SDM options file; the recorded output
granularity reflects those settings.

Retained-quality processing requires complete unmerged reads with qualities.
Use `-mergePreClusterReads 0`, including with profiles that enable it by default.
LotuS disables all three native preprocessing merge modes, then enables seed
merging only for paired input. Retained HQ input and output use quality offset
33. Seed extension keeps already processed sequences without repeating cuts.

Savont, taxonomy-only and demultiplex-only runs do not use `-coarseDerep`.
DADA2 retains `-derepPerSR 1`: main dereplicates remain run-specific, while the
map, `.rest`, and standard HQ files accumulate across runs. LotuS does not look
for run-specific HQ files or concatenate the cumulative HQ files repeatedly.

## Completion metadata and older outputs

After successful preprocessing, `primary/sdm_dereplication.json` records the IO
contract, effective identity, output granularity, quality retention, HQ record
layout/paths, sequencing-run mode and binary file signature. The run manifest
also records granularity, quality retention and HQ layout. Missing retained HQ
files fail before preprocessing is recorded as complete.

Older coarse jobs may contain parent-only HQ files, separate subcluster FASTQs,
and coarse-parent FASTA/map rows. Regenerate preprocessing **and** downstream
clustering/assignments together. A normal full LotuS run regenerates these stages;
it does not resume them from cached intermediates. `-exe 1` is not a resume mode,
and taxonomy redo does not rebuild dereplication. File names or old
`derepStoreQuals` command lines alone cannot identify a compatible cached layout.
Binary `.diff` files cannot recover selected qualities.

The cleanup code already recognizes the `derep.*` family, including both old
and current artifacts. See the [supplied SDM integration contract](sdm_dereplication_io.md)
for full semantics and upstream validation requirements.

## Regression checks

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -p coarse_derep.py -v
```

The tests use real bundled SDM preprocessing and seed extension with controlled
clusterer/mapper output. Coverage includes single/paired input, R2-only variants,
custom retention at 100%, later-variant selection, below-cutoff count recovery,
1,000 fragments at 97%/100% with 1/4/12 workers, cumulative DADA2 run IO, output
metadata, and regeneration of older artifacts. The DADA2 test controls its R
stage; it validates IO integration rather than denoising accuracy.
