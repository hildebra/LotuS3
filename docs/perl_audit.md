# Perl audit and ONT installer update

Audited on 2026-09-11 in the local `dev` checkout, after the ONT port from `ont-amplicon` (`f7978c8`) onto local main/dev (`35b27b0`). The changes described here are local source changes; no commit, push, or production installation was performed.

## Scope

The review covered `lotus3` startup and option selection, configuration resolution, output locking and cleanup, ONT demultiplexing and ASV/backmapping interfaces, contamination filtering, taxonomy hierarchy parsing, abundance tables and BIOM, Lambda query/index handling, and tree/output generation. Existing helpers were reused for subprocess capture, configuration writes, input validation, and dry-run isolation. `helpers/autoInstall.pl` was reviewed for the ONT installation and registration paths; existing database installation logic was retained.

This combines static review, failure-path regression tests, and synthetic pipeline runs. It does not establish that every optional workflow or external algorithm is correct.

## Confirmed defects corrected

| Area | Previous behavior | Correction and evidence |
| --- | --- | --- |
| Output locking | Startup could unlink another run's program log before acquiring its lock. Unlocking and then unlinking the lock file also allowed competing processes to lock different inodes. | Logs are changed only by the lock holder. Lock files persist between runs. Tests hold a competing lock and check log preservation and stable lock-file identity. |
| Replacement of previous output | A recognized output directory was reset before detailed configuration and mapping validation. | Before replacing a nonempty previous run, the lock holder runs the existing isolated dry-run validation in a child process. Tests cover a failed preflight preserving results and a successful preflight allowing replacement. Runtime failures after a successful preflight still do not restore the old run. |
| Cleanup/input overlap | An input, map, configuration, or reference inside a recognized output/temporary directory could be removed during cleanup. | A shared directory guard checks input and configured paths, including symlink locations and targets, before deletion. Tests cover maps, FASTQ, configuration, taxonomy files, configured references, and temporary-directory inputs. |
| PhiX and VXtractor filters | PhiX results were returned under `phiX.0` but read under `phiX`, discarding them. VXtractor's returned exclusion set was not passed to matrix filtering. | PhiX sets are combined under the expected key; VXtractor exclusions are forwarded. A pipeline regression checks that PhiX hits reach the final table and contaminant FASTA while retaining the expected read counts. |
| Contaminant alignment acceptance | Minimap2 hits passed when identity **or** coverage passed; coverage used the reference length, which can be a whole genome. | Both >=50% query coverage and >=90% alignment identity are required. Short exact hits and long low-identity hits fail; valid/boundary hits pass. Malformed PAF records fail with file/line context. Hits are counted once per query. |
| Contaminant mapper choice | The test against `PhiX` was made on a name already suffixed with `.0`, so it always selected minimap2. | Contamination searches honor `-useMini4map`. The VSEARCH route is tested with minimap2 unavailable. Its existing 80% query-coverage threshold remains unchanged. |
| Executable checks | A process could print a plausible version, crash, and still pass. Decimal comparisons incorrectly accepted minimap2 2.9 as newer than 2.17. | SDM/LCA/minimap2/VSEARCH probes check exit status, including signals. Version components are compared numerically. Installer program probes also retain signal failures. Tests cover a crashing executable and older minor versions. |
| Unsupported options | Misspelled clusterers silently selected a default. `-highmem 0` reached a known unsupported backmapping path only after processing. | Unknown clusterers and unsupported dereplication/execution/mapping/demultiplex modes fail early. |
| Shell-sensitive paths | The existing path guard missed quotes, escapes, glob characters, and parentheses used unsafely by legacy shell commands. | Those paths now fail explicitly before processing. Tests cover each character class. This extends the existing restriction; it does not add support for paths containing whitespace. |
| Taxonomy hierarchy parsing | Prefix detection inspected the OTU ID for LCA rows, could duplicate rank prefixes, lost trailing empty fields, and accepted duplicate/short records. Duplicate OTU rows could inflate higher-level counts. | Prefixes are applied per rank without altering OTU IDs. Empty ranks and absent reference-hit fields are preserved as unknown. Duplicate OTU IDs and malformed rows fail. LCA and RDP layouts are tested. |
| Empty-table detection | Zero-count OTUs were subtracted twice: one valid plus one zero-count OTU could abort, while an all-zero table could escape the check. | The check uses the actual surviving rows, including contaminants explicitly retained with `-keepOfftargets`. Both zero-count cases and retained PhiX counts are tested. |
| Lambda input preservation | A `.fna` query, including taxonomy-only user input, was renamed in place before indexing/search; a failure left the original path missing. | Lambda receives a temporary `.fa` copy. A simulated index failure preserves the original name and contents. |
| Database index cleanup | Legacy Lambda upgrade cleanup deleted `DB.*`, potentially removing taxonomy sidecars. `-recalcTaxDB` did not remove the current Lambda3 index. UDB cleanup also used a broad glob. | Only exact current `.lba`, `.lba.gz`, or selected UDB index files are cleared when requested. Old Lambda formats and unrelated sidecars remain intact. Preservation and explicit-rebuild behavior are tested. |
| Tree output directory | `-extendedLogs 0` omitted the directory required for MAFFT output. | Tree building creates the alignment directory when needed; a regression exercises this mode. Higher-level taxonomy output opens also now fail explicitly on errors. |

The contamination fixes can change biological results relative to previous main runs. They repair filter logic and propagation; they do not determine suitable scientific thresholds for every dataset. PAF field interpretation follows the [minimap2 format documentation](https://github.com/lh3/minimap2/blob/master/minimap2.1).

## Installer behavior

Both a normal program installation and `perl helpers/autoInstall.pl --ont-only` install/register minimap2, Savont, and Barbell. The latter avoids repeating database and R-package installation. The actual configuration filename is `lOTUs.cfg`; `configs/LotuS.cfg.def` is the template.

New downloads are pinned to [Savont 0.7.0](https://github.com/bluenote-1577/savont/tree/v0.7.0), [Barbell 0.3.2](https://github.com/rickbeeloo/barbell/tree/v0.3.2), and [minimap2 2.28](https://github.com/lh3/minimap2/releases/tag/v2.28). Compatible existing programs are reused. Archive/binary checksums, staging, executable probes, absolute-path registration, configuration backups, old configurations missing ONT keys, duplicate keys, missing build tools, reruns, and checksum failures are covered by installer tests. A failure preserves that tool's previous configuration entry; successful earlier tool installations remain installed.

The pinned Barbell source uses `<barcode>.trimmed.fastq`, matching the port. Its real Linux x86-64 binary passed version/help probes, as did the pinned minimap2 binary. The installer requires Rust/Cargo >=1.88 and C/C++/CMake for Savont source builds. Details are in [installation](installation.md#adding-ont-tools-to-an-existing-installation).

## Validation

Run all 54 regression tests with:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -p '*.py' -v
```

The suites comprise 8 installer tests, 20 ONT integration/converter tests, and 26 Perl-audit tests. All pass in this environment. Perl syntax checks also pass for `lotus3`, `helpers/autoInstall.pl`, `helpers/autoMap.pl`, and `bin/savont2uc.pl`; `git diff --check` passes.

The complete-output tests, with and without Barbell demultiplexing, use real bundled **SDM 3.43** and **LCA 0.29**, with controlled Savont and alignment output. It checks ASV names, consensus preservation, sample counts (4 and 3), taxonomy, BIOM dimensions/data, higher-rank output, manifest creation, and temporary-file cleanup. Other tests stop after real SDM builds the abundance matrix or invoke specific helpers in temporary script copies. Stand-in tools are confined to temporary tests and are never installed into the checkout.

## Remaining limits and findings

- **A real Savont source build was not performed.** Rust/Cargo and CMake are absent in this environment. The download/build/configuration orchestration is tested using fixture archives and a controlled build command, while real upstream archive checksums and Savont's source CLI were checked separately. macOS/ARM builds have not been executed here.
- **The bundled VSEARCH executable crashes on this host.** `bin/vsearch --version` prints `v2.17.1_linux_x86_64` and exits with signal 11. The new version check detects this. A compatible VSEARCH installation and configuration are needed for real workflows using it. No bundled binary was replaced.
- **Savont's own length limits matter.** The pinned release defaults to 1100–2000 bp even though the SDM presets allow 500–5000 bp. LotuS currently exposes neither Savont's length overrides nor its operon preset. This is documented in [ONT processing](ont.md); shorter/longer amplicons need further integration work or another clusterer.
- **Optional ASV post-clustering remains coupled to phyloseq creation.** In `runPhyloObj`, an absent/failed phyloseq setup returns before `-asvPostClust` runs. That pre-existing behavior was identified but not changed here; separating these R workflows needs its own integration validation.
- **Reference index builds are not coordinated across separate runs.** The output lock protects one output directory. Runs sharing an unbuilt reference database can still attempt to build its index simultaneously. Database-index locking is outside this change.
- Full DADA2, USEARCH, CD-HIT, SWARM, ITSx, VXtractor, R/phyloseq, and phylogeny runs with real external programs were not executed. The built-in multi-workflow `--self-test` was not claimed as passing. Real ONT samples are still needed to assess ASV accuracy, demultiplexing quality, and scientific filtering choices.
