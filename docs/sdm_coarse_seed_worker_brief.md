# SDM coarse seed handoff: superseded

The original SDM worker request has been implemented. Its separate
`.subclusters*.fq` output contract was subsequently replaced by the standard HQ
contract supplied for the updated local SDM 3.53 build.

Use [the current SDM integration contract](sdm_dereplication_io.md) and
[the LotuS coarse-dereplication guide](coarse_dereplication.md). Retained variants
now use `derep.1.hq.fq` and, for pairs, `derep.2.hq.fq`. The main FASTA and map
preserve exact effective R1 dereplicates by default. No companion fix is pending.

The original investigation remains in Git history. Regenerate preprocessing and
downstream clustering together for outputs created under the older contract.
