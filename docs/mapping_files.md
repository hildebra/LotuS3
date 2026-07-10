# Mapping files

LotuS3 uses a mapping file to connect biological sample IDs, input sequence files, primer information and optional metadata. The mapping file is passed with the `-m` flag:

```bash
lotus3 -i reads/ -m map.txt -o myRun
```

The behaviour below reflects the current `readMap()` parser in the `lotus3` Perl pipeline and is cross-checked against the legacy LotuS2 mapping-file documentation: <https://lotus2.earlham.ac.uk/main.php?site=documentation#mappingfile>. The legacy page is still useful because LotuS3 retains much of the same mapping-file vocabulary, but this page documents the stricter LotuS3 parser behaviour where it differs.

## Required format

Mapping files are plain text, tab-delimited files. The mapping file is used by LotuS3 and `sdm` to connect samples to sequence files, barcodes, primers and metadata. The format is close to a QIIME-style mapping file, but LotuS-specific columns are available for sequence-file routing and sample merging.

The first line must be a header line and must start with `#SampleID`:

```text
#SampleID	fastqFile	ForwardPrimer	ReversePrimer	SequencingRun
sampleA	sampleA_R1.fastq.gz,sampleA_R2.fastq.gz	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT	run1
sampleB	sampleB_R1.fastq.gz,sampleB_R2.fastq.gz	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT	run1
```

Important details:

- columns must be separated by **tabs**, not spaces;
- the first header cell must be exactly `#SampleID`;
- the header line must start with `#`;
- empty column headers are not allowed;
- all data rows should have the same number of columns as the header;
- lines starting with `#` after the header are ignored as comments;
- blank or near-empty rows after the header are ignored.

LotuS3 normalises Windows, Unix and classic Mac line endings internally, but tab separation is still required.

## Demultiplexed and multiplexed input

LotuS3 supports two common mapping-file styles.

For **demultiplexed data**, where each sample already has its own FASTQ file or read pair, use `fastqFile` and usually no barcode column is needed. This is the most common layout for current Illumina amplicon projects.

For **multiplexed data**, where several samples are still present in the same FASTQ/FASTA files, include `BarcodeSequence` and primer columns so that `sdm` can identify and demultiplex reads. For dual-indexed data, `Barcode2ndPair` can be used for the barcode on the second read.

Typical older/QIIME-like maps used `SampleID`, `BarcodeSequence` and `LinkerPrimerSequence`. In LotuS3, use `#SampleID` as the first column name, but `LinkerPrimerSequence` remains accepted as a forward-primer column.

## Minimal mapping file with explicit FASTQ files

For demultiplexed FASTQ files, the most common layout is:

```text
#SampleID	fastqFile
sampleA	sampleA_R1.fastq.gz,sampleA_R2.fastq.gz
sampleB	sampleB_R1.fastq.gz,sampleB_R2.fastq.gz
```

Use a comma in the `fastqFile` entry to indicate paired-end reads. For single-end reads, provide one file name only:

```text
#SampleID	fastqFile
sampleA	sampleA.fastq.gz
sampleB	sampleB.fastq.gz
```

When `fastqFile` or `fnaFile` is used, the `-i` argument should point to the directory containing those files, unless the mapping file already contains valid paths.

## Recognised columns

The parser treats the following column names specially.

| Column | Purpose | Notes |
|---|---|---|
| `#SampleID` | Biological sample ID and output table sample name | Required as the first header cell in LotuS3. The older documentation often refers to this field as `SampleID`, but the actual header should be `#SampleID`. |
| `BarcodeSequence` | Barcode/MID assigned to each sample | Used for demultiplexing when reads are not already split by sample. |
| `Barcode2ndPair` | Barcode on the second read | Use for dual-indexed reads where the second read pair carries an additional barcode. |
| `ForwardPrimer` | Forward primer sequence | Can contain IUPAC ambiguous bases. Can be supplied in the map or by command-line `-forwardPrimer`, but not both. |
| `LinkerPrimerSequence` | Legacy/alternative forward-primer column name | Treated as forward-primer information. This column name is common in QIIME-style maps and older LotuS documentation. |
| `ReversePrimer` | Reverse primer sequence | Can contain IUPAC ambiguous bases. Can be supplied in the map or by command-line `-reversePrimer`, but not both. |
| `fastqFile` | FASTQ input file name(s) | Use one file for single-end reads or two comma-separated files for paired-end reads. Paths are interpreted relative to `-i` unless absolute paths are supplied. |
| `fnaFile` | FASTA input file name | Alternative to `fastqFile`; do not define both in the same map. |
| `qualFile` | Quality file corresponding to `fnaFile` | Legacy FASTA+QUAL input mode; rarely used for modern FASTQ workflows. |
| `SampleIDinHead` | Sample identifier already present in FASTA/FASTQ headers | Replaces barcode/MID scanning when the sample ID can be extracted from sequence headers. |
| `MIDfqFile` | Separate FASTQ file containing only MID/barcode reads | Equivalent to the command-line `-barcode`/`-MID` option in older LotuS documentation. Requires paired reads. |
| `SequencingRun` | Sequencing run or batch identifier | Recommended, especially for DADA2 and multi-run projects. Must not contain `NA`. |
| `CombineSamples` | Combine multiple mapping entries into one output sample | Non-empty values are interpreted as the combined sample target. Individual `#SampleID` values must still be unique. |

All other columns are treated as metadata and are carried through the internal mapping table where applicable. The legacy LotuS documentation explicitly notes that the number and names of metadata columns are not limited; only recognised processing-column names activate special pipeline behaviour.

## Sample ID rules

Sample IDs are used as keys throughout the pipeline and should be conservative.

Allowed/recommended characters:

```text
A-Z a-z 0-9 _ .
```

Avoid spaces, dashes and non-ASCII characters.

The current parser applies the following checks:

- empty sample IDs are skipped with a warning;
- duplicate sample IDs abort the run;
- the sample ID `OTU` is not allowed;
- sample IDs with leading or trailing whitespace abort the run;
- sample IDs containing dashes (`-`) abort the run;
- spaces anywhere in the mapping line trigger BIOM-compatibility warnings;
- double quotes (`"`) are removed from mapping values;
- non-ASCII characters trigger BIOM-compatibility warnings.

Recommended examples:

```text
sample_001
patient12_T1
runA.sample03
```

Avoid:

```text
sample-001      # dash is not supported
sample 001      # spaces are problematic
OTU             # reserved name
```

## Primer columns and command-line primers

Primer sequences can be supplied either in the mapping file or on the command line.

Mapping-file example:

```text
#SampleID	fastqFile	ForwardPrimer	ReversePrimer
sampleA	sampleA_R1.fastq.gz,sampleA_R2.fastq.gz	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
```

Command-line example:

```bash
lotus3 -i reads/ -m map.txt -o myRun \
  -forwardPrimer GTGYCAGCMGCCGCGGTAA \
  -reversePrimer GGACTACNVGGGTWTCTAAT
```

Do not provide the same primer type in both places. If `ForwardPrimer` is present in the map and `-forwardPrimer` is also passed, LotuS3 aborts. The same applies to `ReversePrimer` and `-reversePrimer`.

If neither the mapping file nor the command line provides a forward primer, LotuS3 warns that chimera checking may be affected.

## Paired-end and single-end files

In `fastqFile`, paired-end reads are represented by two comma-separated files:

```text
sampleA_R1.fastq.gz,sampleA_R2.fastq.gz
```

Single-end reads are represented by one file:

```text
sampleA.fastq.gz
```

A mapping file should not mix paired-end and single-end entries within the same run. If one row indicates paired-end input and a later row has only one file while paired-end mode is active, LotuS3 aborts with an inconsistent file-number error.

## Barcode and header-based demultiplexing columns

Use `BarcodeSequence` when sample identity is encoded by a barcode/MID sequence. For dual-indexed libraries, `Barcode2ndPair` specifies the barcode on the second read pair.

Use `SampleIDinHead` only when the sample identifier is already present in each FASTA/FASTQ header. In that mode, barcode scanning is replaced by header-based sample assignment.

Use `MIDfqFile` when barcode reads are provided as a separate FASTQ file. This corresponds to the older LotuS `-barcode`/`-MID` option and requires paired reads.

## `SequencingRun`

The `SequencingRun` column is recommended, particularly when using DADA2, because it lets LotuS3 group samples by sequencing run or batch.

Example:

```text
#SampleID	fastqFile	SequencingRun
sampleA	run1/sampleA_R1.fastq.gz,run1/sampleA_R2.fastq.gz	run1
sampleB	run1/sampleB_R1.fastq.gz,run1/sampleB_R2.fastq.gz	run1
sampleC	run2/sampleC_R1.fastq.gz,run2/sampleC_R2.fastq.gz	run2
```

Rules:

- `SequencingRun` values must not be `NA`;
- if the column is present, LotuS3 reports the detected sequencing-run categories;
- if the column is absent, LotuS3 attempts to infer sequencing runs from file or directory structure and writes an updated copy of the map;
- inferred sequencing-run assignments should be checked before relying on them for large or multi-run datasets.

When `SequencingRun` is missing, LotuS3 may infer run groups from:

- the `fastqFile` or `fnaFile` column;
- directory structure in file paths;
- `BarcodeSequence`, if present;
- automatic splitting when no useful run information is available.

## `CombineSamples`

`CombineSamples` can be used when several mapping-file entries should be combined into a single biological sample in the output.

Example:

```text
#SampleID	fastqFile	CombineSamples
sampleA_lane1	sampleA_L001_R1.fastq.gz,sampleA_L001_R2.fastq.gz	sampleA
sampleA_lane2	sampleA_L002_R1.fastq.gz,sampleA_L002_R2.fastq.gz	sampleA
sampleB_lane1	sampleB_L001_R1.fastq.gz,sampleB_L001_R2.fastq.gz	sampleB
```

Rows with an empty `CombineSamples` value are treated as their own sample. Use this column only when the separate rows genuinely represent technical subdivisions of the same biological sample.

## Metadata columns

Additional columns can be included for sample metadata, for example:

```text
#SampleID	fastqFile	Subject	Timepoint	Treatment
sampleA	sampleA_R1.fastq.gz,sampleA_R2.fastq.gz	P01	T0	control
sampleB	sampleB_R1.fastq.gz,sampleB_R2.fastq.gz	P01	T1	treated
```

Metadata values are preserved as provided, except that double quotes are removed. For compatibility with BIOM and downstream tools, avoid spaces, quotes and non-ASCII characters where possible. If spaces are necessary in descriptive metadata, check the downstream output carefully.

## Common templates

### Multiplexed reads with barcode and primer columns

Use this style when reads still need to be demultiplexed by `sdm`:

```text
#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer
S001	ACGTACGT	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
S002	TGCATGCA	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
```

For dual-indexed reads, add `Barcode2ndPair`:

```text
#SampleID	BarcodeSequence	Barcode2ndPair	ForwardPrimer	ReversePrimer
S001	ACGTACGT	AACCGGTT	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
S002	TGCATGCA	CCAATTGG	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
```

### Demultiplexed paired-end Illumina run

```text
#SampleID	fastqFile	ForwardPrimer	ReversePrimer	SequencingRun
sampleA	sampleA_R1.fastq.gz,sampleA_R2.fastq.gz	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT	run1
sampleB	sampleB_R1.fastq.gz,sampleB_R2.fastq.gz	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT	run1
```

### Demultiplexed single-end run

```text
#SampleID	fastqFile	ForwardPrimer	SequencingRun
sampleA	sampleA.fastq.gz	GTGYCAGCMGCCGCGGTAA	run1
sampleB	sampleB.fastq.gz	GTGYCAGCMGCCGCGGTAA	run1
```

### FASTA plus QUAL input 

This format is mainly for legacy data. Prefer FASTQ for new projects.

```text
#SampleID	fnaFile	qualFile	ForwardPrimer	ReversePrimer
S001	S001.fna	S001.qual	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
S002	S002.fna	S002.qual	GTGYCAGCMGCCGCGGTAA	GGACTACNVGGGTWTCTAAT
```

### Multiple sequencing runs

```text
#SampleID	fastqFile	SequencingRun
sampleA	run1/sampleA_R1.fastq.gz,run1/sampleA_R2.fastq.gz	run1
sampleB	run1/sampleB_R1.fastq.gz,run1/sampleB_R2.fastq.gz	run1
sampleC	run2/sampleC_R1.fastq.gz,run2/sampleC_R2.fastq.gz	run2
sampleD	run2/sampleD_R1.fastq.gz,run2/sampleD_R2.fastq.gz	run2
```

### Technical lanes combined into biological samples

```text
#SampleID	fastqFile	CombineSamples	SequencingRun
sampleA_L001	sampleA_L001_R1.fastq.gz,sampleA_L001_R2.fastq.gz	sampleA	run1
sampleA_L002	sampleA_L002_R1.fastq.gz,sampleA_L002_R2.fastq.gz	sampleA	run1
sampleB_L001	sampleB_L001_R1.fastq.gz,sampleB_L001_R2.fastq.gz	sampleB	run1
sampleB_L002	sampleB_L002_R1.fastq.gz,sampleB_L002_R2.fastq.gz	sampleB	run1
```

## Recommended checks before running LotuS3

Before starting a full run, check that:

- the first header cell is exactly `#SampleID`;
- the file is tab-delimited;
- there are no empty headers;
- sample IDs are unique;
- sample IDs do not contain spaces, dashes or non-ASCII characters;
- `fastqFile` and `fnaFile` are not both present;
- `qualFile` is present when using legacy FASTA+QUAL input;
- multiplexed maps contain the required barcode and primer information;
- paired-end entries contain exactly two comma-separated files;
- the run does not mix paired-end and single-end entries;
- all files listed in `fastqFile` or `fnaFile` exist relative to the `-i` directory or as explicit paths;
- `SequencingRun` is present for multi-run datasets;
- `SequencingRun` values are not `NA`;
- primer sequences are provided either in the map or on the command line, not both;
- dual-indexed maps use `Barcode2ndPair` consistently;
- `SampleIDinHead` and barcode-based demultiplexing are not mixed unintentionally.

You can check the installed program and bundled example data with:

```bash
lotus3 --self-test
```

For a new dataset, running LotuS3 first on a small subset of samples is recommended before launching a large full analysis.
