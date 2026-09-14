# Worker brief: lightweight SDM support for ONT amplicons

**Status (2026-09-12): completed in SDM 3.51 beta and integrated into LotuS.** This is the original worker brief, retained for provenance. Current options, installation requirements and tested limits are in the [ONT guide](ont.md). LotuS now enables the new SDM matching mode for ONT read preparation and bundles the validated Linux x86-64 executable.

Please implement a small, usable ONT demultiplexing extension in the SDM repository. SDM is the primary preprocessor for LotuS3 ONT reads. Barbell remains an optional external tool and cannot be modified for this work. We do not need feature parity with Barbell or a major SDM rewrite.

Inspect the current SDM checkout and its AGENTS.md first. Reuse and extend the existing barcode search, primer search, reverse-complement, trimming, statistics and output helpers. Preserve unrelated work. Do not modify LotuS3 or replace its bundled SDM executable as part of this task; return the precise integration requirements when done.

## Bounded objective

Handle ordinary single-end ONT FASTQ/FASTQ.gz amplicons with a known barcode near one or both read ends, optional surrounding adapter sequence, and PCR primers. Support forward and reverse-complement reads. Assign samples conservatively, trim the relevant ends, and continue through SDM's existing outputs and dereplication.

Also preserve the existing one-FASTQ-per-sample workflow: this needs primer processing and count tracking, not barcode matching. Do not claim all ONT kits are supported; document the layouts actually tested.

## Small changes wanted

1. **An opt-in ONT matching mode.** Keep current Illumina/PacBio/default behavior unchanged. Look for existing equivalent options before adding these proposed names:

   | Proposed SDM option | Behavior |
   | --- | --- |
   | `-ontMode 0|1` | Default `0`. Enable the small long-read barcode changes only with `1`; initially support single-end FASTQ and reject incompatible paired/MID modes clearly. |
   | `-barcodeSearchWindow <bp>` | Limit barcode start positions to a bounded region at each relevant read end. Suggested ONT starting default: 200 bp; preserve legacy search behavior outside ONT mode. Reject invalid values and handle reads shorter than the window. |
   | `-ontBarcodeEnds either|both` | Default `either` in ONT mode: accept a unique assignment supported at one or both ends; `both` requires consistent evidence at both ends of the same read. If confidently detected ends disagree, reject the read in either mode. |

   Reuse `maxBarcodeErrs` for the barcode error budget. If it currently counts substitutions only, document that ONT mode interprets it as edit distance, including insertions/deletions. Reuse `maxPrimerErrs` and current primer controls; do not add a second set of quality-filter settings.

2. **Bounded, conservative barcode matching.** Check the current implementation first. Add tolerance for a small number of barcode substitutions/insertions/deletions where missing, preferably using an existing alignment helper. Retain an exact-match fast path and limit approximate matching to the end windows/error budget. A unique best supported sample can be assigned; tied/ambiguous sample assignments must be rejected. Evaluate both read orientations without modifying sequence/qualities repeatedly during candidate selection. Determine trim coordinates before applying the final edits, and reverse qualities whenever reversing a read.

3. **Simple end handling.** Reuse `BarcodeSequence` for the sample barcode, recognizing its appropriate orientation at either end. The first version can limit dual-ended support to the same sample barcode at both ends. Preserve the existing meaning of `Barcode2ndPair` for paired reads; do not silently reinterpret it. Asymmetric/custom barcode pairs are outside this task. Use existing technical-adapter and primer machinery for the common layout `adapter — barcode — spacer/flank — primer — amplicon — reverse primer — optional rear barcode/adapter`. Reject unresolved conflicts rather than inventing a rescue rule. No full kit catalogue is required: mapping files provide the actual barcode and primer sequences.

4. **Keep primer handling and outputs consistent.** For primer-bearing reads, detect/trim the configured primers using the existing machinery. For already-trimmed reads, SDM must work with primer columns absent and primer-presence rejection disabled. LotuS already prepares that map/preset, so do not add a redundant primer-state flag unless a demonstrated SDM limitation requires it. Keep length/quality policy under existing options. Preserve full non-dereplicated FASTQ plus `derep.fas`, `derep.map`, HQ representative reads, sample/run grouping, and later abundance-matrix support. Do not change their formats or accidentally reduce the read set sent to Savont.

5. **A few useful counters.** Extend existing summaries with assigned, unassigned, ambiguous and conflicting-end counts. Existing primer-rejection counters can remain. No new reporting framework or per-read annotation database is required.

## Out of scope

No Barbell pattern language, arbitrary concatemer splitting, complex ligation rescue, complete ONT kit registry, new Rust dependency, broad threading/I/O refactor, or changes to clustering/taxonomy. Primer indel support is worth reusing if already available, but do not turn this into a general primer-aligner rewrite. If a necessary step would require a large redesign, deliver the smaller supported subset and explain its limits.

## Focused validation

Use small deterministic synthetic reads to establish:

- Correct assignment, trimming and synchronized sequence/quality lengths in both orientations, including nonzero adapter offsets and native/rapid-like spacing.
- Exact barcodes and a small substitution/insertion/deletion within the configured budget; rejection beyond that budget, for ties, and for conflicting end assignments.
- The difference between `either` and `both`, including a missing end; sensible behavior at window/read boundaries.
- Primer-bearing input and already-trimmed input, with no second primer removal or rejection in the latter case.
- Two samples with known counts, including duplicated sequences shared across samples: all expected full reads and exact per-sample dereplication counts survive. Confirm compatibility with SDM's existing matrix path.
- Existing non-ONT and paired-end regressions still pass. Do a modest runtime sanity check on a repeated fixture, not a large benchmark project.

LotuS currently defaults to `-ontDemux 0`, calls SDM first, sends full per-sample FASTQ to Savont, and uses SDM dereplication maps for backmapping/counts. It exposes `-ontPrimerState present|removed`. The proposed new SDM flags above are **not yet passed by LotuS**: return the final option names, defaults, minimum SDM version, tested example command/map, and any wrapper changes needed. Prefer a compact patch and a clear list of tested capabilities and remaining limits.
