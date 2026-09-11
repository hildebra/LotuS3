# Oxford Nanopore amplicons

Use `-p ONT` for the ONT workflow. Its default clusterer is **Savont** (`-CL savont`, also `-CL 9`). SDM prepares the reads, Savont calls ASVs, and LotuS backmaps reads to those ASVs to build per-sample abundances. The existing taxonomy and output stages then continue normally.

## Configure the tools

Install or register the ONT dependencies from this checkout:

```bash
perl helpers/autoInstall.pl --ont-only
```

The installer locates or installs [Savont](https://github.com/bluenote-1577/savont), [Barbell](https://github.com/rickbeeloo/barbell), and [minimap2](https://github.com/lh3/minimap2), then registers their absolute paths in the installation-root `lOTUs.cfg`. New downloads are pinned to Savont 0.7.0, Barbell 0.3.2, and minimap2 2.28. Savont source builds require Rust/Cargo >=1.88, C/C++ compilers, and CMake. See [installation](installation.md#adding-ont-tools-to-an-existing-installation) for platform details and manual configuration.

The default configuration includes PATH-based `savont` and `barbell` entries. The installer also adds missing entries to older configuration files. Barbell is optional at runtime when your reads are already demultiplexed. The usual LotuS dependencies and reference databases still apply.

The Barbell command used here is `barbell kit -k <kit> -i <reads> -o <directory> --maximize -t <threads>`. The installer checks for that interface before registering it; upstream development versions may change their command-line options.

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

This selects `configs/sdm_ONT_savont.txt`. SDM demultiplexes where needed and removes primers/barcodes. The preset relaxes SDM quality filtering so Savont receives full, non-dereplicated FASTQ reads and performs its own quality filtering, consensus calling, and chimera removal. It still enforces length and primer rules: the supplied presets allow 500–5000 bp reads and require both primers for the main read set.

Savont 0.7.0 also applies its own **1100–2000 bp** read-length defaults. Consequently, the current Savont route targets full-length SSU amplicons even though the SDM presets permit a wider range. LotuS does not yet expose Savont’s length or rRNA-operon preset options; changing only the SDM length limits does not change Savont’s limits.

Use `-s <file>` to override the SDM preset. If you select another ONT clusterer explicitly, such as `-CL vsearch`, the default is the existing `configs/sdm_ONT_LSSU.txt` quality-filtering preset instead of the permissive Savont preset.

## Raw multiplexed reads

Use one FASTQ file with `-ontDemux barbell`, a kit name, and barcode **labels** in the mapping file:

```text
#SampleID	ONTBarcode	ForwardPrimer	ReversePrimer
sample1	BC01	AGAGTTTGATCMTGGCTCAG	TACGGYTACCTTGTTACGACTT
sample2	BC02	AGAGTTTGATCMTGGCTCAG	TACGGYTACCTTGTTACGACTT
```

```bash
./lotus3 -p ONT -i multiplexed.fastq -m barcodes.tsv -o ont_results \
  -ontDemux barbell -ontKit SQK-RBK114-96 -ontMinReads 100 -t 8
```

Barbell trims barcodes first. LotuS keeps only declared barcode labels, copies their reads to per-sample files, and updates `primary/in.map`. Missing barcodes and samples below `-ontMinReads` are removed from the working map and metadata. Empty or `NA` labels are excluded. The original mapping file is preserved. SDM then uses `configs/sdm_ONT_primersonly.txt` to remove primers before Savont runs.

Each barcode label must belong to one sample. Use sample IDs containing letters, digits, underscores, or dots, starting with a letter, digit, or underscore. The Barbell map must not contain SDM barcode-sequence or separate barcode-file columns, since the barcodes have already been trimmed.

| Option | Default | Purpose |
| --- | --- | --- |
| `-ontDemux` | `0` | Set to `barbell` to preprocess one multiplexed FASTQ; requires `-p ONT`. |
| `-ontKit` | unset | Required Barbell kit name when demultiplexing. |
| `-ontBarcodeCol` | `ONTBarcode` | Column containing labels such as `BC01`. |
| `-ontMinReads` | `0` | Minimum reads per declared barcode; zero keeps every nonempty declared barcode. |
| `-ontWriteFastqCol` | `1` | Add a `fastqFile` column. With `0`, that column must already exist, because SDM requires it for directory input. Existing entries are always updated. |
| `-savontSingleStrand` | `0` | Search both strands; `1` enables Savont's `--single-strand`. |
| `-savontQualCutoff` | `80` | Pass to `--quality-value-cutoff`. |
| `-savontMinBaseQual` | `15` | Pass to `--minimum-base-quality`. |
| `-savontChimeraErrors` | `0` | Pass to `--chimera-allowable-errors`. |

These Savont values preserve the ONT branch defaults; they are not a claim of optimal settings for every dataset. For example, override them with `-savontQualCutoff 95 -savontMinBaseQual 20` as appropriate for your experiment.

ONT backmapping defaults to minimap2's `map-ont` preset, with `-backmap_id 0.97`. `-useMini4map 0` selects the existing VSEARCH/USEARCH mapping route, including contamination searches. SDM retains a separate dereplicated FASTA for backmapping while Savont receives all prepared reads. Savont mode requires single-end FASTQ input and enabled SDM dereplication (`-highmem` nonzero).

Savont already performs de novo chimera removal, so LotuS skips its extra de novo check for this clusterer. The existing optional reference chimera checks still apply. Final outputs use the normal LotuS names, including `OTU.fna`, `OTU.txt`, and `OTU.biom`, with ASV identifiers.

## Validation and regression tests

Add `--dry-run` to your intended command to validate configuration, input paths, options, and barcode mappings without running biological processing or modifying the requested output/report paths. This cannot test the installed Savont/Barbell command interfaces or establish ASV accuracy.

Run the wiring regression suite with:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -p 'ont_integration.py' -v
```

The suite runs real bundled SDM with small synthetic reads and controlled stand-ins for Savont, Barbell, and the mappers. A temporary copy of the pipeline stops after SDM writes the abundance matrix; taxonomy and other unrelated downstream stages are outside these tests. The tests check primer trimming, full-read input, consensus preservation, sample counts, dropped barcodes, option validation, tool failures, and dry-run preservation. Real ONT data and installed external tools are still needed for biological validation.

`bin/savont2uc.pl` is also provided as a standalone converter for Savont's native read-to-ASV assignments. The normal pipeline uses backmapping and does not call this helper.

## Port provenance

The ONT additions were compared against local `main`/`dev` at `35b27b0` and `ont-amplicon` at `f7978c8`, using their common ancestor `43938ce` to isolate the feature work. The port includes the ONT changes to `lotus3`, the two Savont SDM presets, and the standalone converter. `helpers/autoInstall.pl` has no ONT-specific changes relative to that ancestor; `helpers/autoMap.pl` is unchanged between the branches. The local helper implementations and all bundled binaries were retained during the port. The subsequent installer update adds ONT dependency installation and registration to the newer local installer. Main's installer, validation, dry-run, taxonomy, cleanup, and SDM summary-parser improvements remain in place; removed reference-based clustering was not restored.

The follow-up [Perl audit](perl_audit.md) covers installer tests and additional pipeline regressions, including full output generation with real SDM and LCA.
