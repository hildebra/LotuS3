# Oxford Nanopore amplicons

Use `-p ONT` for the ONT workflow. Its default clusterer is **Savont** (`-CL savont`, also `-CL 9`). SDM prepares the reads, Savont calls ASVs, and LotuS backmaps reads to those ASVs to build per-sample abundances. The existing taxonomy and output stages then continue normally.

## Configure the tools

Install or register the ONT dependencies from this checkout:

```bash
perl helpers/autoInstall.pl --ont-only
```

The installer locates or installs [Savont](https://github.com/bluenote-1577/savont), [minimap2](https://github.com/lh3/minimap2) and [Barbell](https://github.com/rickbeeloo/barbell), then registers their absolute paths in the installation-root `lOTUs.cfg`. New downloads are pinned to Savont 0.7.0, Barbell 0.3.2, and minimap2 2.28. Savont source builds require Rust/Cargo >=1.88, C/C++ compilers, and CMake. When suitable Rust/Cargo is unavailable, the installer tries a pinned Bioconda Savont binary, extracted by a Perl helper using `Archive::Tar`, `IO::Uncompress::Unzip`, and the `zstd` command-line tool; Python and Conda are not required. An unusable fallback aborts with instructions to install Rust and compile Savont. See [installation](installation.md#adding-ont-tools-to-an-existing-installation) for platform details and manual configuration.

The full installer includes ONT tools when you accept the initial all-dependencies choice or answer yes to the ONT question in detailed configuration. `--ont-only` installs the three tools directly. Barbell uses a pinned prebuilt release where available. Its use remains optional: SDM is the default preprocessor (`-ontDemux 0`), and `-ontDemux barbell` enables Barbell preprocessing.

This checkout bundles the validated static Linux x86-64 **SDM 3.51 beta** executable in `bin/sdm`. ONT preprocessing requires SDM >=3.51 or a build advertising `ONT amplicon end matching: enabled` in `-version`. LotuS checks the executable selected by the `sdm` configuration entry before processing, because older SDM versions may silently ignore unknown options. If your configuration points elsewhere, update that executable or point it to this checkout's `bin/sdm`. Other platforms require a compatible SDM build.

The default configuration includes PATH-based `savont` and `barbell` entries. The installer adds missing entries for the selected programs to older configuration files. Barbell is optional; SDM remains responsible for read preparation and counts. The usual LotuS dependencies and reference databases still apply.

The default Barbell command is `barbell kit -k <kit> -i <reads> -o <directory> -t <threads>`. LotuS adds `--maximize` only with `-ontBarbellMaximize 1`. The installer checks for the supported command interface before registering Barbell; upstream development versions may change their command-line options.

## Already demultiplexed reads

Supply single-end FASTQ files, with filenames and primers in a tab-separated map:

```text
#SampleID	fastqFile	ForwardPrimer	ReversePrimer
sample1	sample1.fastq	AGAGTTTGATCMTGGCTCAG	TACGGYTACCTTGTTACGACTT
sample2	sample2.fastq	AGAGTTTGATCMTGGCTCAG	TACGGYTACCTTGTTACGACTT
```

Replace the example primers with those used in your experiment. Fields above are separated by tabs.

```bash
./lotus3 -p ONT -i reads/ -m samples.tsv -o ont_results -t 8
```

This selects `configs/sdm_ONT_SAVONT_opt.txt` and automatically enables `-ontMode 1` in the initial SDM call. Later seed/matrix calls and other sequencing platforms do not receive the ONT matching options. SDM demultiplexes where needed, removes primers/barcodes, and quality-filters the reads. Its primary filtered FASTQ passes the exact ambiguity-percentage check before Savont, retaining individual read copies, sample-prefixed identifiers, and base qualities. SDM dereplication is disabled for Savont; other clusterers retain their existing dereplication workflow.

The new Savont preset uses these starting thresholds:

| SDM setting | Default | Purpose |
| --- | --- | --- |
| `minAvgQuality` | `15` | Reject reads with low average Phred quality. |
| `QualWindowWidth` / `QualWindowThreshhold` | `150` / `14` | Reject reads with a low-quality 150-base window. |
| `maxAmbiguousNT` | `5%` | Allow at most 5% non-ACGT bases in each trimmed read. |
| `minSeqLength` / `maxSeqLength` | `1000` / `2000` | Retain trimmed reads within the requested ONT length range. |
| `RejectSeqWithoutFwdPrim` / `RejectSeqWithoutRevPrim` | `T` / `T` | Require both primers for primer-bearing input, then remove them. |
| Quality-based end trimming, homopolymer-tail trimming and fixed-length truncation | Disabled | Preserve complete amplicon boundaries after primer/barcode removal. |

The Savont preset uses average Q15 and a 150-base/Q14 window, less stringent than the general ONT LSSU preset's average Q27 and window Q18. These defaults require evaluation on your data. Savont still applies its own quality filters, ASV inference, and chimera removal. SDM's average Phred threshold and Savont's percentage-based quality cutoff have different meanings.

The percentage form of `maxAmbiguousNT` is a LotuS extension for Savont. SDM itself accepts only integer counts, so LotuS disables that native count limit in a derived options file saved under `primary/` and checks the exact percentage in a streaming pass over the trimmed FASTQ before Savont and backmapping. The same percentage is applied to optional saved per-sample FASTQ. A 1000-base read may contain up to 50 ambiguous bases; a 2000-base read may contain up to 100. Non-ACGT IUPAC symbols count as ambiguous. Retained records remain unchanged and are not dereplicated. `LotuSLogS/savont_ambiguity_filter.log` reports additional rejections after the SDM summary. An integer value in a custom preset retains SDM's original absolute-count behavior.

LotuS explicitly passes **1000–2000 bp** to Savont as well, overriding Savont 0.7.0's native 1100-base minimum. Changing only the SDM preset's length limits does not change these Savont command-line bounds; the route targets full-length SSU amplicons.

Use `-s <file>` to override the SDM preset. Only reads passing the primary SDM filters and any percentage ambiguity limit enter Savont and its abundance counts; lower-quality secondary output is not rescued, even if a custom preset has more permissive `*` settings. If you select another ONT clusterer explicitly, such as `-CL vsearch`, the default remains `configs/sdm_ONT_LSSU.txt`. The older permissive `sdm_ONT_savont.txt` and `sdm_ONT_primersonly.txt` files remain available for explicit use but are no longer selected automatically.

## Primer state entering SDM

Use `-ontPrimerState present` (the default) when PCR primers remain in the reads after any external or Barbell processing. SDM uses the configured primers and preset to find and remove them.

Use `-ontPrimerState removed` when PCR primers have already been removed. This works with both pre-demultiplexed input and `-ontDemux barbell`. LotuS creates `primary/sdm_input.map` without primer columns and generates `primary/sdm_ONT_removed.txt` with primer searching/rejection disabled. It retains the supplied quality/length settings and still prepares individual filtered FASTQ reads and sample counts for Savont. The original map/preset files are preserved; `primary/in.map` retains all sample metadata and `primary/sdm_original_options.txt` records the original preset. Primer-free reads are not reoriented by primer matching, but barcode handling can still change orientation and SDM replaces original header annotations. The default Savont strand-specific mode remains suitable for this handoff.

For example, for already-trimmed sample FASTQs:

```bash
./lotus3 -p ONT -i reads/ -m samples.tsv -o ont_results -ontPrimerState removed
```

The state describes **PCR primers**, not barcode/adapter trimming. It is explicit and is not inferred from a kit name.

## Raw multiplexed reads

SDM is the default demultiplexer. For pooled reads, supply **barcode sequences**, not kit labels, in `BarcodeSequence`, together with the PCR primers:

```text
#SampleID	BarcodeSequence	ForwardPrimer	ReversePrimer
sample1	ACGTCAGTGCTAGACG	AGAGTTTGATCMTGGCTCAG	TACGGYTACCTTGTTACGACTT
sample2	TGCATCGACAGTTCGA	AGAGTTTGATCMTGGCTCAG	TACGGYTACCTTGTTACGACTT
```

Replace these illustrative barcodes and primers with the actual sequences for your experiment. Use literal tabs between fields.

```bash
./lotus3 -p ONT -i multiplexed.fastq.gz -m samples.tsv -o ont_results \
  -barcodeSearchWindow 200 -ontBarcodeEnds either -t 8
```

The supported layout is `adapter — B — spacer — F — insert — revcomp(R) — spacer — revcomp(B) — adapter`, with optional adapters/spacers and either or both barcode ends. The rear barcode is the reverse complement of the same sample barcode. Both read orientations are supported. This is a bounded matcher, without a kit catalogue or support for arbitrary asymmetric barcode pairs, concatemer splitting, or ligation rescue.

`-barcodeSearchWindow` sets the barcode start offsets examined inward from each read end (`0..window-1`, default 200, allowed 1–10000); a complete match can extend beyond the window. Primer start offsets use the same window after barcode trimming. `-ontBarcodeEnds either` accepts unique evidence at one or both ends. `both` requires separate, consistent hits at both ends. Tied assignments, conflicting sample assignments and ambiguous trim endpoints are rejected. Filename-based assignment, including after Barbell, ignores the barcode-end requirement.

The supplied `sdm_ONT_SAVONT_opt.txt` and `sdm_ONT_LSSU.txt` presets set `maxBarcodeErrs` to **1 edit per end**, including substitutions, insertions and deletions. Override it in an SDM options file passed with `-s`: ONT accepts integers 0–3 and A/C/G/T barcodes longer than that budget. `maxPrimerErrs` still means substitutions with IUPAC matching; **primer indels are not supported**. Existing quality and length settings remain in effect, and the supplied ONT presets disable `TruncateSequenceLength`. In demultiplex-only mode (`-saveDemultiplex 1`), LotuS intentionally omits the filtering preset: SDM's default barcode error budget is zero and primer filtering/removal is disabled.

ONT mode requires single-end FASTQ/FASTQ.gz. Paired inputs, separate MID streams, FASTA and alignment inputs are rejected. Remove `Barcode2ndPair`, `MIDfqFile`, `SampleIDinHead` and `alignmentFile` columns entirely, even if their cells are empty. For already demultiplexed reads, omit barcode columns and give each sample its own FASTQ filename. New SDM assignment counters in the filtering logs distinguish assigned, unassigned, ambiguous and conflicting ends; assignment occurs before primer/quality filtering, so it is not the final retained-read count.

### Optional Barbell preprocessing

For optional Barbell preprocessing, use one FASTQ file with `-ontDemux barbell`, a kit name, and barcode **labels** in the mapping file:

```text
#SampleID	ONTBarcode	ForwardPrimer	ReversePrimer
sample1	BC01	AGAGTTTGATCMTGGCTCAG	TACGGYTACCTTGTTACGACTT
sample2	BC02	AGAGTTTGATCMTGGCTCAG	TACGGYTACCTTGTTACGACTT
```

```bash
./lotus3 -p ONT -i multiplexed.fastq -m barcodes.tsv -o ont_results \
  -ontDemux barbell -ontKit SQK-RBK114-96 -ontMinReads 100 -t 8
```

Barbell trims barcodes first. LotuS keeps only declared barcode labels, copies their reads to per-sample files, and updates `primary/in.map`. Missing barcodes and samples below `-ontMinReads` are removed from the working map and metadata. Empty or `NA` labels are excluded. The original mapping file is preserved. SDM uses the same `configs/sdm_ONT_SAVONT_opt.txt` preset to quality-filter the reads and remove any remaining primers. Barcode columns are absent from the prepared map, so SDM does not trim barcodes again. If Barbell has already removed the PCR primers, add `-ontPrimerState removed`.

**Kit-specific limitations:** Barbell's built-in 16S kits include PCR primers in the trimmed flanks. The new `-ontPrimerState removed` mode handles such primer-free output, but it does not fix Barbell 0.3.2's separate filtering problem with expected two-ended 16S reads. Treat `SQK-16S024` / `SQK-16S114-24` as requiring a verified external fix or custom Barbell workflow. See the [ONT demultiplexing audit](ont_barbell_audit.md).

Barbell uses conservative assignment by default. Add `-ontBarbellMaximize 1` to pass `--maximize` and admit its more permissive barcode patterns; this option requires `-ontDemux barbell`.

Barcode labels are kit-specific: native kits use `NBxx`; rapid kit 114-96 mostly uses `BCxx`, with `RBK` labels at positions 26, 39, 40, 48, 54, and 60. Match the labels emitted by the selected Barbell kit.

Each barcode label must belong to one sample. Use sample IDs containing letters, digits, underscores, or dots, starting with a letter, digit, or underscore. The Barbell map must not contain SDM barcode-sequence or separate barcode-file columns, since the barcodes have already been trimmed.

| Option | Default | Purpose |
| --- | --- | --- |
| `-barcodeSearchWindow` | `200` | SDM barcode/primer start-offset window at each end, 1–10000 bp; requires `-p ONT`. |
| `-ontBarcodeEnds` | `either` | SDM accepts one or both consistent barcode ends; `both` requires two. Ignored for filename assignment. |
| `-ontDemux` | `0` | Set to `barbell` to preprocess one multiplexed FASTQ; requires `-p ONT`. |
| `-ontBarbellMaximize` | `0` | With `1`, pass Barbell `--maximize`; requires Barbell demultiplexing. |
| `-ontPrimerState` | `present` | `present`: normal primer handling; `removed`: skip SDM primer searching/rejection. Applies with or without Barbell. |
| `-ontKit` | unset | Required Barbell kit name when demultiplexing. |
| `-ontBarcodeCol` | `ONTBarcode` | Column containing labels such as `BC01`. |
| `-ontMinReads` | `0` | Minimum reads per declared barcode; zero keeps every nonempty declared barcode. |
| `-ontWriteFastqCol` | `1` | Add a `fastqFile` column. With `0`, that column must already exist, because SDM requires it for directory input. Existing entries are always updated. |
| `-savontSingleStrand` | `1` | Enable Savont's `--single-strand` SNPmer mode for SDM-prepared reads; `0` requires strand evidence from both original orientations. |
| `-savontQualCutoff` | `80` | Pass to `--quality-value-cutoff`. |
| `-savontMinBaseQual` | `15` | Pass to `--minimum-base-quality`. |
| `-savontChimeraErrors` | `0` | Pass to `--chimera-allowable-errors`. |

SDM can orient reads during primer/barcode handling, and its sample-prefix rewriting drops original header annotations. Savont recognizes a trailing `rc` annotation to recover original strand evidence; that information is unavailable after SDM. LotuS therefore now defaults to `-savontSingleStrand 1`. This disables the requirement for SNPmer support from both original strands; it does not discard reverse-complement reads. Use `0` only when the actual SDM output retains both original orientations without relying on lost header annotations.

The other Savont values preserve the ONT branch defaults; they are not a claim of optimal settings for every dataset. For example, override them with `-savontQualCutoff 95 -savontMinBaseQual 20` as appropriate for your experiment.

ONT backmapping defaults to minimap2's `map-ont` preset, with `-backmap_id 0.95` for every ONT clusterer: ONT errors are mostly indels, and many reads passing the average-Q15 preset sit just below 97% identity to their ASV. Set `-backmap_id` explicitly to change it. Minimap2 maps the same primary filtered FASTQ supplied to Savont. `-useMini4map 0` selects the existing VSEARCH/USEARCH mapping route, including contamination searches; LotuS converts the filtered reads to FASTA for compatibility with older mapper versions, preserving every read and identifier. No dereplicated FASTA or dereplication map is produced in Savont mode.

SDM builds the abundance matrix from these mappings and the sample prefixes, with Savont ASVs as its seeds. SDM renames the ASVs and its output keeps no link to the Savont IDs, but its seeds are the Savont sequences themselves, so LotuS pairs each matrix row with its ASV by exact sequence rather than by file order. The run stops if any seed is not one of the Savont ASV sequences. Samples with no accepted mapping hits retain zero-count columns. `CombineSamples` grouping is applied after per-sample counting, so equal original read names in different samples remain distinct through mapping. The original sample metadata is preserved.

Savont mode requires single-end FASTQ input. It accepts `-highmem 0`; SDM dereplication stays disabled for Savont regardless of `-highmem`, and `-derepMin` does not apply. Other clusterers continue to require the existing SDM dereplication path. `-saveDemultiplex 2` optionally saves per-sample output, without the SDM base quota, but Savont receives the primary filtered FASTQ directly. The saved files are written as `demultiplexed/*.fq.gz`, compressed after the ambiguity filter. `-saveDemultiplex 1` retains the existing demultiplex-only early exit and does not run Savont or its quality-filtering workflow.

Savont already performs de novo chimera removal, so LotuS skips its extra de novo check for this clusterer. The existing optional reference chimera checks still apply. Final outputs use the normal LotuS names, including `OTU.fna`, `OTU.txt`, and `OTU.biom`, with ASV identifiers.

## Citations

When used, Barbell is cited as Beeloo et al. (2026), *Barbell reveals and resolves demultiplexing and trimming issues in Nanopore data*, Bioinformatics 42(6):btag349, [doi:10.1093/bioinformatics/btag349](https://doi.org/10.1093/bioinformatics/btag349). Savont is cited as Shaw et al. (2026), *Sensitive long-read amplicon sequence variant recovery with savont*, bioRxiv preprint, [doi:10.64898/2026.05.26.727271](https://doi.org/10.64898/2026.05.26.727271).

LotuS records these in `LotuSLogS/citations.txt` when the corresponding stage runs. A Barbell demultiplex-only run cites Barbell without citing Savont.

## Validation and regression tests

Add `--dry-run` to your intended command to validate configuration, input paths, options, and barcode mappings without running biological processing or modifying the requested output/report paths. This cannot test the installed Savont/Barbell command interfaces or establish ASV accuracy.

Run the wiring regression suite with:

```bash
prove -v tests/ont_integration.t
```

The suite runs real bundled SDM with small synthetic reads and controlled stand-ins for Savont, Barbell, and the mappers. A temporary copy of the pipeline stops after SDM writes the abundance matrix; taxonomy and other unrelated downstream stages are outside these tests. The tests check pooled gzip input with offset/substitution/insertion/deletion barcodes, reverse-complement reads and reversed qualities, search-window and two-end rules, initial-call-only ONT flags, SDM version gating, primer trimming, rejection of low-average-quality/low-window-quality/ambiguous/primer-missing reads, exclusion of secondary-quality reads, empty-input handling, retention of duplicate reads, shared read names across samples, zero-hit samples, grouped counts, VSEARCH FASTA conversion, dereplication for other clusterers, multiple-ASV consensus/count preservation, dropped barcodes, option validation, tool failures, and dry-run preservation. Real ONT data and installed external tools are still needed for biological validation.

`bin/savont2uc.pl` is also provided as a standalone converter for Savont's native read-to-ASV assignments. The normal pipeline uses backmapping and does not call this helper.

## Port provenance

The ONT additions were compared against local `main`/`dev` at `35b27b0` and `ont-amplicon` at `f7978c8`, using their common ancestor `43938ce` to isolate the feature work. The port includes the ONT changes to `lotus3`, the two Savont SDM presets, and the standalone converter. `helpers/autoInstall.pl` has no ONT-specific changes relative to that ancestor; mapping-file creation (then `helpers/autoMap.pl`, now `lotus3 -create_map`) is unchanged between the branches. The local helper implementations and all bundled binaries were retained during the port. The subsequent installer update adds ONT dependency installation and registration to the newer local installer. The SDM integration follow-up replaces `bin/sdm` with the validated static Linux x86-64 3.51 beta build from the local SDM checkout. Main's installer, validation, dry-run, taxonomy, cleanup, and SDM summary-parser improvements remain in place; removed reference-based clustering was not restored.

The follow-up [Perl audit](perl_audit.md) covers installer tests and additional pipeline regressions, including full output generation with real SDM and LCA.
