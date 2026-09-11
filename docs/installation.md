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

The autoinstaller downloads and installs required software and databases inside the LotuS3 directory. Full installations include minimap2, Savont, and Barbell for ONT processing. Database-only modes do not install programs.

## Adding ONT tools to an existing installation

From this checkout, run:

```bash
perl helpers/autoInstall.pl --ont-only
```

This installs or locates the three ONT tools and writes their absolute executable paths to **`lOTUs.cfg`** in the installation root. Existing compatible executables configured there, in `bin/`, or on `PATH` are reused. Missing configuration entries are added, duplicate entries are collapsed, unrelated settings are retained, and the previous configuration is backed up as `lOTUs.cfg.bak` before the first change in that installer run. This mode does not repeat database or R-package installation.

New downloads use pinned versions: Savont 0.7.0, Barbell 0.3.2, and minimap2 2.28. The installer verifies SHA-256 checksums, builds in temporary directories, and checks executable/version or command-interface compatibility before registering each tool. An installation failure preserves that tool's previous configuration entry; tools installed successfully earlier in the run remain installed.

Savont is built from source and needs **Rust/Cargo 1.88 or newer, C and C++ compilers, and CMake**. Install these build tools before running the installer, or provide a compatible Savont executable on `PATH`. Barbell uses an official binary on Linux x86-64, Linux ARM64, and Apple Silicon; other supported hosts build it with Cargo. Minimap2 uses its Linux x86-64 binary where applicable and is otherwise built with Make, a C compiler, and zlib development headers. Downloads and source builds require network access, and archive extraction needs `tar`, gzip, and bzip2 support. A downloaded binary that cannot execute on the host is rejected; a compatible tool supplied on `PATH` can then be used instead.

For manual installations, add the following entries to `lOTUs.cfg` (or your custom configuration selected with `-c`):

```text
savont /absolute/path/to/savont
barbell /absolute/path/to/barbell
minimap2 /absolute/path/to/minimap2
```

Barbell is needed at runtime only for `-ontDemux barbell`. See [ONT processing](ont.md) for commands and input formats. The ONT-only installer updates the installation-root `lOTUs.cfg`; it does not update custom `-c` configuration files.

## Requirements

LotuS3 requires:

- Perl 5;
- a C++ compiler supporting C++17;
- R and RScript;
- Java or OpenJDK for tools such as the RDP classifier;
- selected third-party tools and databases, depending on the chosen workflow.

Conda normally handles its package dependencies. The source autoinstaller expects system build tools to be installed already; see the ONT build requirements above.

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
