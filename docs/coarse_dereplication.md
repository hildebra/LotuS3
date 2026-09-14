# Coarse dereplication

`-coarseDerep X` enables SDM coarse dereplication with a fractional identity
threshold from **0.95 through 1.0**, inclusive. For example, add
`-coarseDerep 0.975` for 97.5% identity. Omitting the flag retains the existing
pipeline. Explicit `1.0` still enables exact-variant quality retention.

**This integration awaits a companion SDM seed-reader fix.** It requires SDM
3.52 or newer advertising `-seedSubclusters` in `-help_flags`. The currently
bundled SDM 3.51 and the investigated upstream 3.52 beta cannot run the complete
workflow. LotuS checks this before processing. See the
[SDM worker brief](sdm_coarse_seed_worker_brief.md) for the required fix.

LotuS converts the fraction to SDM's percentage `-derepIdentity` and enables
`-derepStoreQuals 1`. The coarse parent representatives remain the input to
clustering/backmapping. Seed extension instead receives the reconstructed exact
variants and their selected quality vectors:

- Single-end: `tmpFiles/derep.subclusters.fq`.
- Paired: `tmpFiles/derep.subclusters.1.fq` and `derep.subclusters.2.fq`.

These exports include parents below the dereplication abundance cutoff. Sample
counts still come from `derep.map`; the companion seed reader must count each
parent once while considering all of its variants. Selected qualities are
observed vectors, not averages. Native seed merging remains enabled. Binary
`.diff` output is unnecessary for this handoff.

The first integration supports the global SDM-dereplication workflows (UPARSE,
VSEARCH, CD-HIT, SWARM, and UNOISE). DADA2's per-sequencing-run dereplication,
Savont's individual-read input, and taxonomy-only/demultiplex-only runs are not
supported with this flag. Use `-mergePreClusterReads 0`; preprocessing merging
is incompatible with retained-quality variants. This also applies to profiles
that enable preprocessing merging by default.

SDM compares substitutions from the same start across the shorter capped
sequence, independently for both mates. Its identity threshold does not perform
an indel alignment or reverse-complement search. See
[SDM's paired-mode documentation](https://github.com/hildebra/sdm/blob/master/docs/paired-diff-quality.md)
for the matching rules.
