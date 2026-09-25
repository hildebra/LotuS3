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

The autoinstaller downloads and installs selected software and databases inside the LotuS3 directory. It first asks:

```text
Install LotuS3 with all possible dependencies (all databases, ITS, ONT related workflows)?
Simply enter or "1" for yes, "0" for detailed configuration via question.
```

Press **Enter** or enter **1** to select all dependencies supported by the installer: KSGP, SILVA, GreenGenes2, HITdb, PR2 and beeTax reference databases; ITS/UNITE and ITS chimera resources; UTAX databases; both BLAST and Lambda; Savont and Barbell; and the standard programs and R packages. This skips the individual package-selection questions. The existing SILVA license acceptance and any incompatible-R-version question still apply. System prerequisites such as Rscript and compilers must already be available; the installer checks them before any download (see [requirements](#requirements)).

Enter **0** for detailed configuration. It asks about similarity-search programs, reference databases, ITS, UTAX, and then:

```text
-- ONT -- Install ONT tools (Savont and Barbell)?
Simply enter or "1" for yes, "0" for no.
```

Declining ONT skips Savont/Barbell installation and their build checks, and preserves their existing configuration entries. Minimap2 remains part of the core installation because it is also the default mapper for non-ONT workflows. Barbell is included whenever ONT tools are selected; its use during processing remains optional.

An existing installation first offers its refresh/link menu. Refreshing programs then uses the same all-dependencies or detailed choice; refreshing only databases asks only database questions. `--condaDBinstall` remains noninteractive and installs databases only. End-of-input aborts rather than accepting a default.

## Adding ONT tools to an existing installation

From this checkout, run:

```bash
perl helpers/autoInstall.pl --ont-only
```

This installs or locates minimap2, Savont and Barbell and writes their absolute executable paths to **`lOTUs.cfg`** in the installation root. Existing compatible executables configured there, in `bin/`, or on `PATH` are reused. Missing configuration entries are added, duplicate entries are collapsed, unrelated settings are retained, and the previous configuration is backed up as `lOTUs.cfg.bak` before the first change in that installer run. This mode does not repeat database or R-package installation.

The explicit `--ont-only` mode installs all three ONT dependencies without asking the interactive configuration questions. SDM is the primary ONT preprocessing path. Barbell is included using a pinned prebuilt release where available; enable it at runtime with `-ontDemux barbell`. The earlier `--with-barbell` option is still accepted for compatibility, but the ONT tools choice governs installation in detailed mode.

New downloads use pinned versions: Savont 0.7.0, Barbell 0.3.2, and minimap2 2.28. The installer verifies SHA-256 checksums, builds in temporary directories, and checks executable/version or command-interface compatibility before registering each tool. An installation failure preserves that tool's previous configuration entry; tools installed successfully earlier in the run remain installed.

If Savont is not already installed, the installer builds it from source when **Rust/Cargo 1.88 or newer** is available; this also needs C and C++ compilers and CMake. If Cargo or a suitable Rust compiler is unavailable, it instead downloads a pinned **Savont 0.7.0 Bioconda binary** for Linux x86-64/ARM64 or macOS Intel/Apple Silicon. The Perl helper `helpers/extract_conda_executable.pl` extracts the `.conda` package using **`Archive::Tar`, `IO::Uncompress::Unzip`, and the `zstd` command-line tool**. The Perl modules are distributed with Perl; `zstd` must be available on `PATH`. Neither Python nor Conda is required.

The Bioconda fallback checks that the Perl extractor and `zstd` can run, verifies the package SHA-256, and tests `savont --version` and `savont asv --help`, including the options LotuS uses, before registering the executable. These binaries depend on system libraries and may not run on every host. If the download, extraction, or executable check fails, the installer **aborts and tells the user to install Rust >=1.88, Cargo, C/C++ compilers, and CMake to compile Savont**. It preserves the previous Savont configuration entry.

Barbell uses an official binary on Linux x86-64, Linux ARM64, and Apple Silicon; other supported hosts build it with Cargo. Minimap2 uses its Linux x86-64 binary where applicable and is otherwise built with Make, a C compiler, and zlib development headers. Downloads and source builds require network access, and archive extraction needs `tar`, gzip, and bzip2 support. A downloaded binary that cannot execute on the host is rejected; a compatible tool supplied on `PATH` can then be used instead.

For manual installations, add the following entries to `lOTUs.cfg` (or your custom configuration selected with `-c`):

```text
savont /absolute/path/to/savont
barbell /absolute/path/to/barbell
minimap2 /absolute/path/to/minimap2
```

Barbell is needed at runtime only for `-ontDemux barbell`. See [ONT processing](ont.md) for commands and input formats. The ONT-only installer updates the installation-root `lOTUs.cfg`; it does not update custom `-c` configuration files.

## Requirements

LotuS3 requires:

- Perl 5.14 or newer. Only modules that ship with Perl are used, so no CPAN modules need to be installed. Some distributions split the standard modules off the interpreter: on Fedora, RHEL and derivatives install the `perl` package (not only `perl-interpreter`); on Debian and Ubuntu the standard `perl` package is enough;
- a C++ compiler supporting C++17;
- R and RScript;
- Java or OpenJDK for tools such as the RDP classifier;
- selected third-party tools and databases, depending on the chosen workflow.

Conda normally handles its package dependencies. The source autoinstaller expects system build tools to be installed already. Before it downloads any database, a full installation checks for:

- `tar`, `gzip`, `unzip`, `make` and a C compiler (`gcc` or `cc`);
- `xz`, when Lambda is selected on Linux;
- a working `sdm` and `LCA`. The bundled `bin/sdm` and `bin/LCA` are Linux x86-64 builds; on other systems, place builds from [sdm](https://github.com/hildebra/sdm) and [LCA](https://github.com/hildebra/LCA) at those paths first (see [manual sdm compilation](#manual-sdm-compilation)).

A missing item stops the installer with a list of what to install, before anything is downloaded. Java is only reported, since it is needed at run time for RDP classification. Downloads need `wget` or `curl`.

## Download verification

Every file the autoinstaller downloads is checked against a SHA-256 checksum pinned in `helpers/autoInstall.pl`, and all downloads use HTTPS. A file whose checksum differs is deleted and the installation stops, so a changed or substituted upstream file is never installed. To move a tool or database to a new release, download it, check it, and update its URL and checksum together in the `%PINNED_SHA256` table.

Archives already present in `bin/installs/` (IQ-TREE and MAFFT for Linux) are used instead of downloading them again when their checksum matches. The macOS BLAST+ and Clustal Omega files that the installer used to fetch are no longer available upstream; on macOS the installer skips them with a warning. Install `blastn`/`makeblastdb` separately if you need BLAST; alignments use MAFFT.

`bin/vsearch` is a statically linked Linux x86-64 build of VSEARCH 2.32.0 and is registered directly when it runs. It was compiled from the upstream release source (`vsearch-2.32.0.tar.gz`, SHA-256 `99578a8b960a0fb87c1f19dc65aedecddc01cfa91851b697dac8294dd08a6ceb`) with `./configure LDFLAGS=-static`, `make` and `strip`. The resulting executable has SHA-256 `8472f7852e4f320f1cc67e5dc09507e2e31034cd962bd1b8979bec894711c5a3`. On other platforms the installer downloads the pinned VSEARCH 2.32.0 release for Linux ARM64 (static) or macOS.

## Installer options

```text
perl helpers/autoInstall.pl [options]
  (no option)        interactive installation; a rerun offers to refresh databases/programs
  --ont-only         install or register only the ONT tools (minimap2, Savont, Barbell)
  -condaDBinstall    non-interactive download of the standard database set (Bioconda installs)
  -downloadLmbdIdx   download prebuilt Lambda indices instead of building them
  -lambdaIndex       build Lambda indices for the installed databases
  -link_usearch PATH register an existing USEARCH binary and exit
  --no-telemetry     do not send the installation ping (see below)
  --with-barbell     accepted for compatibility; Barbell is part of the ONT tools
  -h, --help         show this help
```

With `-lambdaIndex` on a fresh installation, Lambda indices are built once Lambda itself has been installed, later in the same run.

After a full or database installation, the installer sends one request to the LotuS server with a random installation ID (the `UID` entry in `lOTUs.cfg`) and the LotuS and sdm versions. No file names, paths or data are sent. Use `--no-telemetry` to skip it.

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

LotuS3 includes a statically compiled Linux x86-64 **SDM 3.53 beta** binary at `bin/sdm`. ONT preprocessing requires SDM >=3.51 or a build whose `-version` output contains `ONT amplicon end matching: enabled`; coarse dereplication requires SDM >=3.52. LotuS checks the configured executable before processing. When using an older installation, update the executable selected by the `sdm` entry in `lOTUs.cfg`, even if ONT tools are already installed. The bundled binary should work on most Linux x86-64 systems. On macOS, or if the bundled binary does not work on your system, compile `sdm` from its source repository and place it at `bin/sdm`:

```bash
git clone https://github.com/hildebra/sdm.git sdm_src
make -C sdm_src
cp sdm_src/sdm bin/sdm
```

`LCA` is built the same way from https://github.com/hildebra/LCA into `bin/LCA`. With the sources checked out as `sdm_src/` and `LCA_src/` in the LotuS3 directory, the autoinstaller compiles them itself when the bundled binaries do not run.

## Updating a GitHub installation

If LotuS3 was installed with `git clone`, update the code with:

```bash
git pull
```

Then rerun the autoinstaller:

```bash
perl helpers/autoInstall.pl
```

On an existing installation it offers to (1) refresh databases and reinstall the secondary programs, (2) refresh only the databases, or (3) set the USEARCH path. It does not update LotuS3 itself: the former online updater (`-forceUpdate`) was removed, because it installed an archive it could not verify. sdm, LCA, rtk and the ONT tools are reused when the existing executables still pass their checks; the refresh options download the selected databases and other programs again.

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
