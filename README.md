<p align="center">
  <img src="images/lotus3_highres.png" alt="LotuS3 logo" width="350"/>
</p>

<p align="center">
  <a href="https://anaconda.org/bioconda/lotus3"><img src="https://anaconda.org/bioconda/lotus3/badges/downloads.svg" alt="Bioconda downloads"/></a>
  <a href="https://anaconda.org/bioconda/lotus3"><img src="https://anaconda.org/bioconda/lotus3/badges/latest_release_relative_date.svg" alt="Latest Bioconda release"/></a>
</p>

# LotuS3

LotuS3 is a lightweight, fast and configurable pipeline for amplicon sequencing analysis. It supports common marker-gene targets including 16S, 18S and ITS, and can be configured for additional amplicons. LotuS3 integrates several OTU/ASV inference methods, including DADA2, UPARSE, UNOISE3, CD-HIT and VSEARCH, and provides multiple options for taxonomic annotation. Outputs include abundance tables, representative sequences, taxonomy assignments, run logs and citation files, and can be imported into R or other downstream analysis tools.

Full documentation is available in this repository under [`docs/`](docs/). The historical documentation is available at: http://lotus2.earlham.ac.uk/

## Contents

- [Installation](#installation)
- [Check your installation](#check-your-installation)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Configured example run](#configured-example-run)
- [Documentation](#documentation)
- [Updating LotuS3](#updating-lotus3)
- [Citation](#citation)
- [Acknowledgements](#acknowledgements)

## Installation

### Recommended installation: Bioconda

The recommended way to install LotuS3 is through Bioconda:

```bash
conda create -n LotuS3 -c conda-forge -c bioconda --strict-channel-priority LotuS3
conda activate LotuS3
```

If you already have a suitable conda environment, you can also install LotuS3 directly into it:

```bash
conda install -c conda-forge -c bioconda LotuS3
```

LotuS3 requires Perl 5, a C++17-capable compiler, R and Java. These dependencies are normally handled by the conda installation or by the GitHub autoinstaller described below.

For more details, see [`docs/installation.md`](docs/installation.md).

### Developer or pre-release installation: GitHub

The GitHub repository may contain pre-release versions and development updates that are not yet available through Bioconda.

```bash
git clone https://github.com/hildebra/LotuS3.git
cd LotuS3
perl helpers/autoInstall.pl
```

The autoinstaller downloads and installs the required software and databases inside the LotuS3 directory. If you want to install software or databases to custom locations, see [`docs/installation.md`](docs/installation.md).

## Check your installation

After installation, run the built-in self-test:

```bash
lotus3 --self-test
```

If you installed LotuS3 directly from GitHub and are running it from the repository directory, the executable may be invoked as:

```bash
./lotus3 --self-test
```

The self-test checks the LotuS3 executable, helper tools, selected databases and example workflows. Some optional database indices may be generated during the first full pipeline run; missing optional indices should be treated as warnings rather than fatal installation errors.

If the self-test completes successfully, your installation is ready for normal use. See [`docs/troubleshooting.md`](docs/troubleshooting.md) if the self-test fails.

## Inputs

A typical LotuS3 run requires:

- an input directory containing FASTQ files;
- a mapping file linking sample IDs to sequencing files and sample metadata;
- optionally, an `sdm` read-filtering configuration file;
- optionally, forward and reverse primer sequences;
- optionally, a selected reference database and taxonomy assignment method.

Mapping files define how sequencing files are assigned to biological samples. See [`docs/mapping_files.md`](docs/mapping_files.md) for details.

## Outputs

LotuS3 produces the main files required for downstream amplicon analysis, including:

- OTU or ASV abundance tables;
- representative OTU/ASV sequences;
- taxonomy annotations;
- quality-filtering and read-processing summaries;
- run logs;
- method-specific citation files.

The exact output files depend on the selected clustering method, reference database, taxonomy aligner and post-processing options. Each run records relevant settings and citations in the LotuS3 output directory. See [`docs/outputs.md`](docs/outputs.md) for more detail.

## Configured example run

The self-test is the recommended first check after installation. Once the self-test succeeds, you can run LotuS3 on the bundled example data using an explicitly configured 16S/DADA2/SILVA workflow:

```bash
./lotus3 -i Example/ \
  -m Example/miSeqMap.sm.txt \
  -o myTestRun \
  -s configs/sdm_miSeq2.txt \
  -p miSeq \
  -amplicon_type SSU \
  -forwardPrimer GTGYCAGCMGCCGCGGTAA \
  -reversePrimer GGACTACNVGGGTWTCTAAT \
  -CL dada2 \
  -refDB SLV \
  -taxAligner lambda
```

This example explicitly sets read filtering, defines the data as 16S Illumina MiSeq data, removes the supplied PCR primers, uses DADA2 for ASV inference, and annotates ASVs against the SILVA reference database using the lambda aligner.

Building the lambda-formatted SILVA reference database can take a long time the first time this workflow is run. Make sure that the SILVA database was selected during installation if you want to run this example.

LotuS3 has many additional command-line options, but the defaults are chosen to work well for common amplicon workflows. To list available options, run:

```bash
./lotus3
```

For conda installations, the executable may be available on your path as `lotus3`. To locate the executable and associated shared files, use:

```bash
which lotus3
```

## Documentation

The README is intentionally short. Detailed documentation is split across the [`docs/`](docs/) folder:

- [`docs/installation.md`](docs/installation.md) - installation routes, dependencies, autoinstaller use and manual `sdm` compilation;
- [`docs/mapping_files.md`](docs/mapping_files.md) - mapping file purpose, minimal structure and common checks;
- [`docs/outputs.md`](docs/outputs.md) - expected output categories and where to find run-level information;
- [`docs/custom_reference_databases.md`](docs/custom_reference_databases.md) - custom FASTA/taxonomy databases and multi-database annotation;
- [`docs/pacbio.md`](docs/pacbio.md) - PacBio CCS/HiFi amplicon processing;
- [`docs/troubleshooting.md`](docs/troubleshooting.md) - installation, self-test and runtime troubleshooting;
- [`docs/citations.md`](docs/citations.md) - LotuS3, third-party software and database citations.

## Updating LotuS3

If LotuS3 was installed with `git clone`, update the code with:

```bash
git pull
```

LotuS3 also has a built-in update mechanism through the autoinstaller. If LotuS3 was first installed with:

```bash
perl helpers/autoInstall.pl
```

then running the autoinstaller again checks for updates. Previously downloaded proprietary programs and databases do not need to be downloaded again. If no updates are available, the autoinstaller exits without making changes, so it can be run periodically.

More installation and update details are available in [`docs/installation.md`](docs/installation.md).

## Citation

If you use LotuS3, please cite:

**Pipeline** - Özkurt E, Fritscher J, et al. (2022) LotuS2: An ultrafast and highly accurate tool for amplicon sequencing analysis. *Microbiome* 10:176. doi:10.1186/s40168-022-01365-1.

**Off-target removal** - Bedarf JR, Beraza N, Khazneh H, Özkurt E, et al. (2021) Much ado about nothing? Off-target amplification can lead to false-positive bacterial brain microbiome detection in healthy and Parkinson's disease individuals. *Microbiome* 9:75.

LotuS3 writes method-specific citations for each run to:

```text
LotuSLogS/citations.txt
```

Please also cite the clustering, taxonomy and database tools used in your specific run. See [`docs/citations.md`](docs/citations.md) for the full citation list.

## Acknowledgements

LotuS3 was developed at Quadram Institute Bioscience (QIB) and Earlham Institute (EI), Norwich, UK. Various members of the Hildebrand group contributed to the pipeline; see Özkurt et al. (2022).

(c) Falk Hildebrand, Falk.Hildebrand {at} gmail.com
