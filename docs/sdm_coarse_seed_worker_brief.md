# SDM worker brief: use coarse subclusters in seed extension

## Request

Implement native seed extension from SDM's reconstructed exact-subcluster FASTQ
exports. LotuS3 now has a guarded `-coarseDerep X` integration, where
`0.95 <= X <= 1.0`; it passes `-derepIdentity (100 * X)` and
`-derepStoreQuals 1` during preprocessing. Its seed step needs the explicit new
SDM input flag **`-seedSubclusters 1`**, advertised in **`-help_flags`**.

This flag is a proposed interface for this fix, not an existing SDM 3.52 option.
The current LotuS3 integration requires SDM >= 3.52 and checks for that flag before
processing. Keep the default `0` and preserve ordinary HQ-seed behavior.

Investigated upstream commit: `ae8c73cdc2896775ff8999de9d4e64a507096718`
(`sdm 3.52 beta`, checked 2026-09-14).

References:

- [Paired quality retention](https://github.com/hildebra/sdm/blob/master/docs/paired-diff-quality.md)
- [Coarse options](https://github.com/hildebra/sdm/blob/master/docs/options.md#streaming-coarse-dereplication)
- [Seed workflows](https://github.com/hildebra/sdm/blob/master/docs/workflows.md)

## Confirmed failure

A real portable build (`make WITH_HTS=0 STATIC=0 -j4`) was tested with three
200-base single-end reads: two identical sequences and one one-substitution
variant, all assigned to sample `s1`. At 97% identity, SDM exported:

```text
# derep.map
#SMPLS  0:s1
r2;size=3;  0:3

# derep.subclusters.fq headers
@r2.sub1;size=2;
@r2.sub2;size=1;
```

(The map uses tabs.) Supplying that FASTQ to the existing seed extension with
the parent mapping and `-derep_map derep.map` exits with status 226 on Linux
(`exit(738)` in the code), reporting:

```text
r2;size=3; is not r2.sub1;size=2;
```

`UClinks::oneDerepLine` consumes one map row per FASTQ record and checks
`sameHead`. `findSeq2UCinstruction` then matches the entire read identifier
against the mapper query. Neither assumption holds for `.subN` records.
Simply removing `.subN`, or relaxing the header check, does not solve repeated
map consumption, missing candidates, or repeated abundance counting.

## Required input/output contract

Preprocessing remains:

```sh
sdm -i_fastq reads.1.fq,reads.2.fq -paired 2 \
  -map samples.tsv -options options.txt \
  -o_dereplicate derep.fas -min_derep_copies 1 \
  -derepIdentity 97 -derepStoreQuals 1 -derepPerSR 0 \
  -merge_pairs_derep 0 -merge_pairs_filter 0 -merge_pairs_demulti 0
```

Retaining qualities already retains nucleotide differences in memory. Do not
require binary `.diff` output (`-derepStoreDiffs 1`) or FASTA subcluster export.
Binary differences contain neither qualities nor per-variant sample IDs.

After the existing clusterer/backmapper produces parent-to-OTU assignments,
LotuS invokes the equivalent of:

```sh
sdm -i_fastq derep.subclusters.1.fq,derep.subclusters.2.fq -paired 2 \
  -seedSubclusters 1 -merge_pairs_seed 1 \
  -derep_map derep.map -optimalRead2Cluster clusters.uc \
  -uparseVer vsearch -OTU_fallback otu.fna \
  -ucAdditionalCounts additional.uc -ucAdditionalCounts1 rest.uc \
  -options options.txt -o_qual_offset 33 \
  -o_fastq otu_seeds.1.fq,otu_seeds.2.fq -otu_matrix OTU.txt
```

For single-end input, use `derep.subclusters.fq`, `-paired 1`, and normal
single-end seed output. Use each pipeline's existing mapper-format argument;
the example above is VSEARCH. PAF, UPARSE, UNOISE, and CD-HIT assignments also
need to retain their existing interpretation.

Required semantics:

1. Read every exact subcluster, retaining full processed sequence and its
   selected qualities. Resolve the generated final `.subN` suffix to its parent,
   matching the parent's actual map/assignment key, including its total
   abundance where the existing mapper requires it. Do not indiscriminately
   strip periods or size annotations from unrelated identifiers.
2. Apply the parent's existing accepted OTU assignment(s) to its subclusters.
   Consider **all** variants as candidates through existing seed-selection and
   pair-merging code. Do not select only the first variant or silently fall back
   to the original coarse HQ representative. This step does not independently
   remap variants or claim variant-specific sample assignments.
3. Read and apply each parent's sparse sample-count vector **once**. Child
   `;size=N;` values describe variant multiplicity, not additional parent counts.
   Their sum should equal the parent's count. For the example above, the matrix
   must contain 3 counts, not 6; a paired fragment is counted once as well.
4. Keep paired FASTQs synchronized by exact generated ID, count, and order.
   Both qualities must come from the corresponding paired observation. Distinct
   R2 sequences with identical R1 must remain distinct candidate pairs.
5. Preserve ordinary additional-count handling, including `.rest` parents and
   sample-prefixed medium-quality reads, without lost or doubled counts.
   Exported FASTQs include parents below `min_derep_copies`. Preserve existing
   chimera-count splitting and fallback behavior. Unmapped parents must not
   create new OTUs or acquire a made-up mapping.
6. Seed merging stays enabled independently of preprocessing merging. Avoid
   reapplying barcode/primer cuts, quality trimming, orientation, or truncation
   to already processed subcluster sequences. Preserve full-length seed output.
7. Reject missing parents, truncated FASTQ records, mismatched/missing mates,
   repeated or malformed variant IDs, and inconsistent variant/parent counts
   with actionable errors. Reuse existing input validation where possible.

## Implementation guidance

Extend `UClinks` and its helpers in `containers.cpp` / `containers.h`; avoid
reimplementing mapper parsers, quality scoring, or read merging in LotuS.
Relevant functions include:

- `findSeq2UCinstruction`, `oneDerepLine`, `readDerepInfo`, `finishMAPfile`;
- `getMAPPERline`, `uclInOldDNA`, `besterDNA`, `add2OTUmat`;
- `finishUCfile`, `addUCdo`, and the seed finalization path in `IO.cpp`.

Separating candidate evaluation from abundance addition in the existing helper
would let every variant compete while counting its parent once. A grouped
streaming reader or parent-assignment lookup can avoid retaining all expanded
FASTQ sequences at once. Account for mapper order differing from export order.
The old `.map`/HQ behavior and mapper parsers should remain the default path.

Document the new flag in built-in help and workflow documentation. LotuS's
capability probe runs `sdm -help_flags` and searches for `-seedSubclusters`.

## Required validation

Use real SDM preprocessing followed by real seed extension. Include:

- Single-end and paired inputs at 95%, 97.5%, and explicit 100% identity.
- Several samples and coarse parents, repeated exact variants, and an R2-only
  variant. Check exact per-sample matrix totals independently of the output.
- A deliberately better later subcluster that changes the selected seed,
  demonstrating that every candidate is evaluated. Check source qualities and
  full-length paired merging, including a truncation/search cap during derep.
- `.rest` parents, medium-quality extra mappings, unmapped parents, chimera
  counting, reordered parent mappings, parent IDs containing `.sub` themselves,
  and malformed/missing input.
- Legacy HQ input with the new flag omitted: unchanged counts and seed behavior.
- Existing `make test` and paired/diff-quality tests.

The LotuS regression suite can then be run with:

```sh
PYTHONDONTWRITEBYTECODE=1 LOTUS_TEST_SDM=/path/to/patched/sdm \
  python3 -m unittest discover -s tests -p coarse_derep.py -v
```

It already tests real 3.52 preprocessing and captures the seed command. The
native seed/count test is skipped until the executable advertises the new flag;
with the fix, that test must run and pass for both single and paired input.
Report the resulting commit, flag/help contract, and executable path/version to
the LotuS workflow. The bundled LotuS `bin/sdm` is still 3.51 and has not been
replaced; update it or configure the tested executable after the upstream fix.
