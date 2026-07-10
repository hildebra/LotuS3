# Troubleshooting

This page lists common LotuS3 installation, self-test and runtime issues.

## Run the self-test first

After installation, run:

```bash
lotus3 --self-test
```

or, from a GitHub checkout:

```bash
./lotus3 --self-test
```

The self-test checks the executable, helper tools, selected databases and example workflows. Some optional database indices are generated during the first relevant full run; missing optional indices should usually be warnings rather than fatal errors.

## Installation troubleshooting

If you install from GitHub, running the autoinstaller should usually be sufficient:

```bash
perl helpers/autoInstall.pl
```

Common installation issues include:

### Warning about RScript

Install or activate a working version of R and RScript.

### Warning about Java or OpenJDK

Java is required for the RDP classifier and some post-filtering steps.

### `Can't locate FindBin.pm in @INC`

Perl cannot find a module that should be in `PERL5LIB`. Check your Perl include paths with:

```bash
perl -e 'print join(":", @INC), "\n"'
```

Then make sure any shell configuration, such as `.bashrc`, does not overwrite the expected Perl paths.

### RScript error in `LULU.R`

R may try to download packages during first use. This can fail without an internet connection.

### DADA2 installation problems

DADA2 can often be installed more reliably through conda/mamba than through the autoinstaller:

```bash
mamba install -c conda-forge -c bioconda bioconductor-dada2
```

The LotuS3 autoinstaller works well inside a conda environment. If the autoinstaller stops at a problem that you then fix manually, re-running it should continue from the previous point.

## Runtime troubleshooting

- The example workflows may require approximately 4 GB of memory. If less memory is available, the run may stop with an out-of-memory error.
- Some optional database indices are generated during the first relevant run. Missing indices should usually be warnings unless the selected workflow requires them immediately.
- If a command works manually but fails inside `--self-test`, compare the executable path, current working directory, database path and number of requested CPU threads.
- If paired-end input fails, check that forward and reverse FASTQ files contain the same number of records and that sample names match the mapping file.
- If reads are unexpectedly removed during PacBio processing, check whether strict primer rejection is active and whether `-forwardPrimer` and `-reversePrimer` were provided.
