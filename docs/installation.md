# Installation

LotuS3 can be installed either through Bioconda or directly from GitHub. Bioconda is recommended for most users. The GitHub route is useful for development versions, pre-release features or local testing.

## Recommended installation: Bioconda

Create a dedicated conda environment:

```bash
conda create -n LotuS3 -c conda-forge -c bioconda --strict-channel-priority LotuS3
conda activate LotuS3
```

Alternatively, install LotuS3 into an existing suitable environment:

```bash
conda install -c conda-forge -c bioconda LotuS3
```

Bioconda normally handles the core dependencies, including Perl, R and Java-related packages required by common LotuS3 workflows.

## Developer or pre-release installation: GitHub

The GitHub repository may contain updates that are not yet available through Bioconda:

```bash
git clone https://github.com/hildebra/LotuS3.git
cd LotuS3
perl helpers/autoInstall.pl
```

The autoinstaller downloads and installs required software and databases inside the LotuS3 directory.

## Requirements

LotuS3 requires:

- Perl 5;
- a C++ compiler supporting C++17;
- R and RScript;
- Java or OpenJDK for tools such as the RDP classifier;
- selected third-party tools and databases, depending on the chosen workflow.

These dependencies are normally handled by conda or by the LotuS3 autoinstaller.

## Installing dependencies manually

If you need to prepare an environment manually, the following packages are commonly useful:

```bash
conda install -c bioconda r-base usearch wget perl rdp_classifier
```

For DADA2, conda/mamba installation is often more reliable than ad hoc installation through R:

```bash
mamba install -c conda-forge -c bioconda bioconductor-dada2
```

## Manual `sdm` compilation

LotuS3 includes a statically compiled Linux `sdm` binary, which should work out of the box on most Linux systems. On macOS, or if the bundled binary does not work on your system, compile `sdm` manually:

```bash
cd sdm_src
make
cp sdm ../sdm
```

The autoinstaller can also compile `sdm` when required.

## Updating a GitHub installation

If LotuS3 was installed with `git clone`, update the code with:

```bash
git pull
```

LotuS3 also has a built-in update mechanism through the autoinstaller. If LotuS3 was first installed with:

```bash
perl helpers/autoInstall.pl
```

then running the autoinstaller again checks for updates. Previously downloaded proprietary programs and databases do not need to be downloaded again. If no updates are available, the autoinstaller exits without making changes.

## Checking the installation

After installation, run:

```bash
lotus3 --self-test
```

or, from a GitHub checkout:

```bash
./lotus3 --self-test
```

The self-test checks the executable, helper tools, selected databases and example workflows. Missing optional database indices may be reported as warnings if they can be generated during the first full run.
