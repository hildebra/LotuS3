# Coarse dereplication as internal storage

`-coarseDerep X` selects SDM's internal coarse grouping at fractional identity
**0.95 through 1.0**, inclusive. For example:

```sh
./lotus3 -i Example/ -m Example/miSeqMap.sm.txt -o coarse_run \
  -CL vsearch -coarseDerep 0.97
```

This option is intended for memory management during preprocessing. It must not
change the biological grouping exported to clustering or the seed candidates.
With otherwise identical settings and input, ordinary and coarse runs must
export the same exact dereplicates, representative sequences and qualities, and
sample abundances. Representative tie/order effects must be distinguished from
changes in candidate membership. Peak memory still needs measurement on the
workload of interest; selecting coarse storage alone is not proof of savings.

## Ordinary output and seed extraction

The same paths are used in both modes:

| File | Contents |
| --- | --- |
| `derep.fas` | Passing exact dereplicates; FASTQ when required by DADA2 |
| `derep.merg.fas` | Passing merged dereplicates when preprocessing merging is enabled |
| `derep.fas.rest` | Below-cutoff dereplicates for additional-count recovery |
| `derep.map` | Exact dereplicate sample counts, including below-cutoff records |
| `derep.1.hq.fq` | One selected R1 representative per output dereplicate, with qualities |
| `derep.2.hq.fq` | Corresponding R2 representatives for paired input |

The ordinary seed extraction command consumes these HQ files and the existing
mapping files. `-coarseDerep` does **not** enable `-seedSubclusters`, force R1-only
search, disable preprocessing merging, or change the main FASTA/FASTQ format.
Search caps, truncation, prefix consolidation, quality filtering, copy thresholds,
and `-mergePreClusterReads` retain their ordinary meaning. Explicit `1.0` also
uses the ordinary output contract.

`derepStoreQuals=0` does not discard the selected representative's qualities.
Those qualities are retained for HQ output. Main dereplicated FASTQ quality
averaging is a separate operation and must also match ordinary processing.

## Optional retained-variant output

Explicitly setting `derepStoreQuals=1` in a custom SDM options file requests a
different seed-candidate policy: every retained exact full sequence/pair variant
is exported with selected observed qualities and a `parent.subN` ID. LotuS
preserves this explicit setting whether or not `-coarseDerep` is supplied, and
selects `-seedSubclusters 1 -i_qual_offset 33` accordingly.

This optional variant policy can produce far more HQ records than ordinary
representative output. It requires unmerged preprocessing and
`-mergePreClusterReads 0`. It is not required to retain ordinary representative
qualities or to use coarse internal storage. R2-only differences can produce
multiple full-pair variants even when the ordinary search key is identical.

## Effective options and compatibility

`-coarseDerep` passes its identity as a percentage via `-derepIdentity`, preserves
the configured `derepStoreQuals` setting (default zero), and forces
`-derepStoreDiffs 0 -derepSubclusterFasta 0 -derepReassign 0 -derepCoarseClusters 0`.
The last two settings keep the output at exact-dereplicate granularity.

Use an SDM build with standard HQ paths and exact output support. LotuS checks
version >=3.52 and the standard HQ / `-derepCoarseClusters` capability markers.
`-seedSubclusters` support is required only for explicit variant retention.
Version/capability labels alone cannot prove output parity; development builds
can share labels. The corrected SDM exact-output path must accumulate ordinary
quality evidence for main FASTQ output, including merged quality evidence.

Savont, taxonomy-only and demultiplex-only runs do not use `-coarseDerep`.
DADA2 keeps its existing main FASTQ format and `derepPerSR` behavior: run-specific
main dereplicates and cumulative map, rest and HQ files.

After preprocessing, `primary/sdm_dereplication.json` records the effective
identity, representative/variant layout, quality-retention and merge/search
settings, output paths, and executable/options hashes. Ordinary HQ qualities
are present even when `quality_retention` (the explicit variant-retention flag)
is zero. The run manifest records these policies too.

Earlier LotuS coarse runs forced variant retention and disabled preprocessing
merging. Regenerate preprocessing and downstream clustering/assignments together
to compare the corrected workflow. A normal full run regenerates these stages;
`-exe 1` is not a resume mode and taxonomy redo does not rebuild dereplication.

## Validation

```sh
python3 -m unittest discover -s tests -p coarse_derep.py -v
python3 tests/audit_derep_counts.py /scratch/run/derep.fas
```

The integration regression compares an ordinary run against coarse identities
97% and 100%, with one and multiple workers. It captures native preprocessing
files before later consumers modify them, then compares HQ sequences/qualities,
seed commands, final seeds, seed statistics and the abundance matrix. It covers
single-end input and paired input with preprocessing merging both off and on.
Other regressions cover explicitly requested variant retention, quality/cut
settings, missing mates, per-run output and count recovery. The clusterer/mapper
fixtures are controlled; SDM preprocessing and seed selection use the real binary.

See [the SDM IO contract](sdm_dereplication_io.md) for details of optional retained
variant output and [the worker brief](sdm_coarse_seed_worker_brief.md) for the
required storage/output separation.
