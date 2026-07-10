# Mapping files

LotuS3 uses a mapping file to connect sequencing files to biological sample IDs and sample metadata. The mapping file is passed with the `-m` flag:

```bash
lotus3 -i Example/ -m Example/miSeqMap.sm.txt -o myRun
```

## Purpose

The mapping file tells LotuS3:

- which sample IDs should appear in the final output tables;
- which sequencing files belong to each sample;
- whether reads are single-end, paired-end or associated with a particular sequencing run;
- which optional sample metadata should be retained for downstream interpretation.

## General format

Mapping files are tab-delimited text files. A typical mapping file contains a header row followed by one row per sample or sample-file entry.

Commonly used columns include identifiers such as sample ID, sequencing run and FASTQ file names. The exact required columns can depend on the input layout and workflow. Use the example mapping files distributed with LotuS3 as templates, for example:

```text
Example/miSeqMap.sm.txt
```

## Recommended checks before running LotuS3

Before starting a full run, check that:

- every sample ID is unique unless the workflow explicitly expects repeated entries;
- FASTQ file names in the mapping file match files in the input directory;
- paired-end samples have matching forward and reverse read files;
- sample names do not contain problematic whitespace;
- tab characters, not spaces, separate columns;
- metadata columns do not accidentally duplicate reserved column names.

The built-in self-test checks the bundled example mapping file:

```bash
lotus3 --self-test
```

For user datasets, input validation should fail early when mapping file sample IDs do not match input sample names or when paired-end FASTQ files are inconsistent.
