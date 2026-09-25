# LotuS3 regression tests

This directory holds the developer regression suite. It is separate from `lotus3 --self-test`, which checks an installation with the bundled example data (see the [main README](../README.md#check-your-installation)).

## Requirements

- **Linux x86-64.** The tests run the bundled `bin/sdm` and `bin/LCA`, which are Linux builds. On Windows, run the suite inside WSL (see below).
- **Perl 5.14 or newer with its standard modules only** (`Test::More`, `JSON::PP`, `Archive::Tar`, `IO::Compress`, `Digest::SHA`, `Module::CoreList`). `prove` is part of Perl. `perl_dependencies.t` fails if any LotuS3 or test file loads a module that is not part of Perl.
- No Python, R, network access or reference databases are needed.
- **Optional: `zstd`.** Twenty Bioconda `.conda` extraction cases are skipped without it. Put `zstd` on `PATH`, or set `LOTUS_TEST_ZSTD=/path/to/zstd`.

## Running the tests

From the repository root:

```bash
prove tests/
```

A complete run reports `Files=9, Tests=174`. With the bundled SDM 3.53 beta, one case fails: `test_storage_only_preserves_main_fastq_quality_averaging` in `coarse_derep.t`. It runs SDM directly and shows a quality-averaging difference inside SDM, not a LotuS3 defect.

To run one file with per-case output:

```bash
prove -v tests/ont_integration.t
```

To test a different SDM executable instead of the bundled `bin/sdm`:

```bash
LOTUS_TEST_SDM=/path/to/sdm prove tests/
```

### Windows (WSL)

Open a WSL shell (for example Ubuntu), change to the checkout, and run `prove`:

```bash
cd /mnt/c/path/to/LotuS3
prove tests/
```

Running from the Windows filesystem works but is slow (about four minutes). For a faster run (about one minute), copy the checkout into the WSL home directory first:

```bash
mkdir -p ~/lotus3-tests
tar --exclude=.git --exclude=.claude -cf - . | tar -xf - -C ~/lotus3-tests
cd ~/lotus3-tests && chmod +x bin/sdm bin/LCA && prove tests/
```

## What is tested

| File | Cases | Coverage |
| --- | ---: | --- |
| `ont_integration.t` | 55 | ONT preprocessing with SDM, Savont and Barbell wiring, barcode matching, ASV-to-abundance pairing, backmapping identity and mapped-read report, option validation, and the `bin/savont2uc.pl` converter |
| `perl_audit.t` | 31 | Output locking and cleanup guards, input validation, version checks, contamination filters, backmapping read counts, taxonomy parsing, USEARCH version/chimera handling, and complete pipeline runs |
| `audit_fixes.t` | 11 | Command splitting and failure reporting, failing IQ-TREE runs, empty `SequencingRun` values and the copied map header, whitespace and dashes in sample names, SINTAX rank letters, FASTA/taxonomy ID agreement, `-xtalk`, option validation, path-depth guards and `-create_map` pairing |
| `sintax_taxonomy.t` | 5 | SINTAX taxonomy (`-taxAligner sintax`, `-refDB SINTAX`): OTU-first `hiera_BLAST.txt`, higher-level tables, BIOM and phyloseq input, taxonomy-only and ITS runs, early failures |
| `coarse_derep.t` | 28 | Coarse dereplication storage parity, retained-variant seeds, the seed-extension failure abort, and the dereplication count auditor |
| `installer_ont.t` | 25 | ONT tool installation and registration, and Bioconda package extraction |
| `installer_options.t` | 13 | Interactive installer prompts and install modes, the stderr-only Rscript version and the preflight before any download |
| `installer_downloads.t` | 5 | Pinned SHA-256 checksums for every download URL, refusal of unpinned and mismatching downloads, reuse of bundled archives, `--help` and the retired `-forceUpdate` |
| `perl_dependencies.t` | 1 | Every module loaded by `lotus3`, the helpers, `bin/` scripts and the tests is core Perl 5.14 |

Shared code lives in `lib/`:

- `LotusTest.pm`: the pipeline fixture and its stand-in tools.
- `InstallerTest.pm`: the installer fixture.
- `DerepAudit.pm`: the dereplication count audit.

Each case builds a fresh temporary installation. Pipeline cases run the real bundled SDM (and LCA where needed) on small synthetic reads. Savont, Barbell, minimap2, VSEARCH, USEARCH, Lambda, MAFFT, FastTree, Rscript (phyloseq helper) and the installer's download/build tools are replaced by small Perl stand-ins that record how they were called. Most pipeline cases use a copy of `lotus3` that stops once SDM has written the abundance matrix. Nothing is written into the checkout.

The stand-ins check how LotuS3 calls and connects these tools. They do not establish biological accuracy, which needs real data and the real programs.

## Auditing real dereplication output

`audit_derep_counts.pl` checks that an SDM dereplication output set is internally consistent:

```bash
perl tests/audit_derep_counts.pl /path/to/run/derep.fas > count-audit.json
```

It exits with 0 on success, 1 when the audit fails, and 2 on a usage error. See [coarse dereplication](../docs/coarse_dereplication.md) for context.

## Adding tests

- In the relevant `.t` file, add `case test_name => sub { my $t = shift; ... };`. Each case receives a fresh fixture. `tests/lib/LotusTest.pm` provides `run_lotus`, `probe`, `tool_calls` and `check_counts`.
- Write any stand-in programs in Perl (`#!/usr/bin/env perl`). LotuS3 ships no Python, and `*.py` files are ignored by git.
- Keep tests offline and independent of reference databases.
