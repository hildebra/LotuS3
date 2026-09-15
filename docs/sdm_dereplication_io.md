# LotuS3 integration: standard dereplication IO

This is the integration contract for the SDM build containing the standard HQ output update. Apply these instructions in the LotuS3 repository; this SDM change does not edit LotuS3 itself. The SDM version label alone may not distinguish earlier development builds. Its `-help_flags` output describes retained qualities as `Standard HQ output: .1.hq.fq; paired additionally .2.hq.fq`.

## Storage-only coarse mode (current LotuS policy)

`-coarseDerep` changes internal SDM grouping only. It must preserve ordinary exact
outputs, selected representative qualities, main FASTQ quality averaging and
seed extension. It no longer enables `derepStoreQuals=1` or disables preprocessing
merging. The default uses `derepStoreQuals=0` and the ordinary seed reader.
Representative HQ qualities are still retained with this setting.

The retained-variant commands below describe a separate, explicitly requested
`derepStoreQuals=1` mode. They must not be enabled automatically for memory-saving
coarse storage. See [the current coarse guide](coarse_dereplication.md).

## Files to use

Given `-o_dereplicate /scratch/run/derep.fas`, use the same paths for ordinary and coarse dereplication:

| Path | Pipeline role |
| --- | --- |
| `/scratch/run/derep.fas` | Passing representatives for the existing clustering/mapping stage. Default coarse mode preserves exact R1 dereplicates. |
| `/scratch/run/derep.merg.fas` | Additional passing merged representatives when ordinary dereplication merging is enabled; include these in passing-abundance audits. |
| `/scratch/run/derep.fas.rest` | Below-cutoff parents for the existing additional-count/recovery stage. |
| `/scratch/run/derep.map` | Parent sample abundances, including below-cutoff parents. |
| `/scratch/run/derep.1.hq.fq` | R1 seed candidates, also the filename for single-end input. |
| `/scratch/run/derep.2.hq.fq` | Synchronized R2 seed candidates for paired input. |

SDM derives the stem by removing the last filename extension. Keep the main FASTA output for clustering; if a consumer explicitly needs main FASTQ, supply `-derep_format fq` and an appropriate `-o_dereplicate` filename. A `.fq` extension alone does not select FASTQ. The HQ files are always FASTQ.

By default, coarse groups are internal search and compression references. The main FASTA and `.map` retain exact effective search-key dereplicates, including the configured search/truncation cap and ordinary prefix consolidation when selected. With quality retention, the search source is currently R1. R2 differences do not split that parent; the existing paired better-seed selection can replace both mates. Minimum-copy rules apply to these exact dereplicates. Select `-derepCoarseClusters 1` only when the output should describe coarse clusters instead. `-derepReassign 1` also explicitly selects coarse-cluster output and performs final reassignment. `-derepGlobal 1` controls search scope, not output granularity.

With `-derepStoreQuals 0`, HQ files contain one selected representative per parent, as before. With `-derepStoreQuals 1`, they contain all retained exact variants and their selected quality vectors. The paired records represent full available `(R1,R2)` variants after physical preprocessing cuts: identical R1 with different R2 can share a parent while remaining separate seed candidates. No additional parent-only HQ copy or `.subclusters.fq`/`.subclusters.1.fq`/`.subclusters.2.fq` is written. There is no file-level split between cluster reads and subreads: the representative's own 100% variant, other retained variants, and singleton clusters are all written together. Each exact variant has its own count and appears once; never append a second cluster-count record for its representative.

## Explicit retained-variant preprocessing commands

1. Use the updated SDM executable for both preprocessing and seed extension.
2. Keep `-o_dereplicate`, `-derepPerSR`, copy thresholds, sample delimiter, output quality offset, and normal filtered-output paths. Existing `-suppressOutput` settings do not suppress dereplication files.
3. For coarse processing with retained exact seed candidates, pass `-derepStoreQuals 1` and the intended `-derepIdentity`. Complete unmerged reads with qualities are required. Disable native preprocessing merging with `-merge_pairs_derep 0 -merge_pairs_filter 0 -merge_pairs_demulti 0`. Pair merging can still be enabled during seed extension. This currently changes ordinary merged search/admission; it is not a behavior-preserving substitute for the original merged workflow.
4. Keep `-derepStoreDiffs 0 -derepSubclusterFasta 0 -derepReassign 0 -derepCoarseClusters 0` unless those optional outputs or reassignment are required. Binary differences are not required for retained-quality output or seed extension. Default partitioned coarse processing adds no final cluster-merging pass.
5. Replace any coarse-specific FASTQ path selection with the normal `.1.hq.fq` and `.2.hq.fq` paths. Update expected-output checks, file-size checks, cleanup/archive lists, and resume metadata together. Do not wait for `.subclusters` FASTQ files.

An explicitly requested retained-variant paired command is:

```sh
sdm -i_path INPUT_DIRECTORY -map primary/in.map -options sdm_miSeq.txt \
    -o_dereplicate /scratch/run/derep.fas -derep_format fa \
    -o_fna /scratch/run/demulti.1.fna,/scratch/run/demulti.2.fna \
    -paired 2 -derepIdentity 97 -derepStoreQuals 1 \
    -derepStoreDiffs 0 -derepSubclusterFasta 0 -derepReassign 0 -derepCoarseClusters 0 \
    -merge_pairs_derep 0 -merge_pairs_filter 0 -merge_pairs_demulti 0 \
    -min_derep_copies 8:1,4:2,3:3 -dere_size_fmt 0 \
    -sample_sep ___ -derepPerSR 0 -o_qual_offset 33 \
    -suppressOutput 1 -threads 12
```

Retain LotuS3's existing demultiplexed/filtered output and log arguments as needed; the example only shows the dereplication integration. Search/truncation settings continue to come from the same options and command line.

## Behavioral limits and count audit

Default retained output now uses ordinary R1 admission and ordinary prefix consolidation/recovery. It still cannot combine native merged search/output with paired quality retention. Both ordinary and retained HQ files now preserve full available mates and observed qualities independently of `TruncateSequenceLength` and logical quality trimming. Barcode/primer/adapter and fixed-base cuts remain applied; use `keepBarcodeSeq=0` during preprocessing when seeds must not contain barcode sequence. Reconstructed seed input does not repeat these cuts. Older retained HQ exports may lack tails; regenerate them before seed extension. Pairs with a physically empty R2 remain excluded and are reported. See the [four-mode comparison and required storage/IO changes](derep-behavior-parity.md) before treating retention as equivalent to the original pipeline.

Audit parent abundance as `map = main + merged + rest`, counting each paired fragment once. The new finalized dereplication abundance line reports this partition; the filtering log's low-quality-recovery counter is not a complete count ledger. Run `python3 tests/audit_derep_counts.py /scratch/run/derep.fas` on each fresh complete output set (`derepPerSR=0`).

Record merge flags, `noSearchWithMerge`, search/truncation limits, retention/output policy, and the executable hash in preprocessing cache metadata. Changing these requires regenerating preprocessing and downstream assignments together.

## Seed-extension command changes

Select the reader from the **effective quality-retention setting used to create the files**, not from `derepIdentity` or the filenames:

```text
if preprocessing used derepStoreQuals=1:
    seedSubclusters = 1
else:
    seedSubclusters = 0

R1 = stem + ".1.hq.fq"
R2 = stem + ".2.hq.fq"  # paired input only
```

This rule also applies at `derepIdentity=100` when quality retention is explicitly enabled. Pass the generated parent `.map`, the existing parent-to-OTU assignments, and the original fallback OTU FASTA. For paired retained variants:

```sh
sdm -i_fastq /scratch/run/derep.1.hq.fq,/scratch/run/derep.2.hq.fq \
    -paired 2 -seedSubclusters 1 \
    -derep_map /scratch/run/derep.map \
    -optimalRead2Cluster parent_assignments.uc \
    -OTU_fallback clustered_otus.fna -options sdm_miSeq.txt \
    -sample_sep ___ -i_qual_offset 33 -o_qual_offset 33 \
    -merge_pairs_seed 1 -o_fna OTU.fna -otu_matrix OTU.txt
```

Keep the pipeline's existing assignment format, reference/fallback paths, identity/coverage settings, additional-count arguments, logs, and final-output choices. `-derep_map` is the **input** abundance map for seed extension; preprocessing derives its output map from `-o_dereplicate`.

For single-end retained variants, supply only `derep.1.hq.fq` with `-paired 1 -seedSubclusters 1`, and omit pair merging. For ordinary representative HQ input, use `-seedSubclusters 0` or omit that flag.

Variant headers have the form `parent.subN;size=M;`. The native variant reader resolves the final `.subN` suffix to the parent, evaluates every variant as a seed, and applies the parent's sample vector once. Keep paired files synchronized and each parent's `.sub1`, `.sub2`, ... order intact. Do not strip suffixes, assign every variant the full parent count, independently sort mates, or add R1 and R2 abundances together. Do not filter variant records by the parent minimum-copy threshold: the map and HQ files include below-cutoff parents needed by recovery. Existing seed extension avoids repeating preprocessing cuts, truncation, or quality trimming.

## Sequencing runs and resuming older jobs

With `-derepPerSR 1`, the main output receives a sanitized run label, such as `derep.runA.fas`. The normal `derep.map`, `derep.fas.rest`, `derep.1.hq.fq`, and paired `derep.2.hq.fq` accumulate across completed runs. Do not look for run-suffixed HQ files or concatenate the cumulative HQ files once per run. Optional binary differences and exact FASTA exports retain their run suffixes.

This output correction also changes coarse-parent FASTA/map rows into exact-R1 rows by default. Regenerate preprocessing and downstream clustering/assignments together when resuming outputs from earlier coarse builds.

Record the effective output granularity (`exact R1` or `coarse clusters`), quality-retention setting and the HQ record layout (`representatives` or `exact variants`) in the preprocessing completion/cache metadata. Earlier development builds used the same HQ filenames for parent-only records and separate `.subclusters` FASTQ files for variants. Regenerate preprocessing and dependent outputs when migrating those cached jobs, or explicitly keep the old matching reader/path contract for the entire old job. Never infer that an old HQ file contains variants merely because its preprocessing command included `derepStoreQuals=1`. Generated `.diff` files cannot restore retained qualities.

## Integration checks

- For ordinary exact mode without retention, standard filenames and representative seed behavior remain unchanged.
- Compare 100% and 97% on the same accepted effective R1 keys: main FASTA, `.rest`, `.map`, and HQ records must agree, ignoring record order when seed winners are unambiguous. The all-passing 1,000-pair regression covers this at 1, 4, and 12 workers. The mixed-quality behavioral comparison additionally reconciles failed/alternate-tier recovery and passing/rest totals; it records the remaining empty-mate difference from ordinary R1 processing.
- At 97% with quality retention, every exact full pair appears once in the two HQ files with the same ID/count and the selected observed quality vectors; parent FASTA/map outputs remain available.
- For each parent, HQ variant abundances sum to the parent map abundance. One paired fragment contributes one count, not two. Below-cutoff parents remain available.
- A parent with multiple variants can select a better later variant during seed extension without adding its sample counts again. R2-only variants remain linked to their corresponding R1 qualities.
- Validate single-end input, paired input, `derepStoreQuals=1` at 100%, reassignment if enabled, and cumulative sequencing-run output.
- A fresh run produces no `.subclusters` FASTQ files. Optional `.diff` and `.subclusters*.fna` files only appear when requested.

The SDM regressions for these contracts are `tests/test_diff_quality_pipeline.py`, `tests/test_paired_diff_quality_pipeline.py`, `tests/test_derep_behavior_parity.py`, `tests/test_seed_subclusters.py`, and the output/run cases in `tests/TestCoarseDerep.cpp`.
