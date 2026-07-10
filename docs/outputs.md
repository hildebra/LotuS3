# Outputs

LotuS3 writes the main files required for downstream amplicon analysis into the selected output directory.

```bash
lotus3 -i Example/ -m Example/miSeqMap.sm.txt -o myRun
```

In this example, output files are written under `myRun`.

## Main output categories

Depending on the selected workflow, LotuS3 produces:

- OTU or ASV abundance tables;
- representative OTU/ASV sequences;
- taxonomy annotations;
- quality-filtering and read-processing summaries;
- run logs;
- method-specific citation files;
- intermediate files useful for troubleshooting or reproducibility.

## Output variability

The exact output files depend on:

- the clustering or inference method, for example DADA2, UPARSE, UNOISE3, CD-HIT or VSEARCH;
- whether the run produces ASVs or OTUs;
- the selected reference database;
- the selected taxonomy assignment method;
- read filtering, chimera filtering and post-filtering settings;
- whether the run uses Illumina, PacBio CCS/HiFi or another amplicon strategy.

## Logs and citations

Each run records settings and method-specific information in the output directory. LotuS3 also writes citations for the tools and databases used in a run to:

```text
LotuSLogS/citations.txt
```

Use this file to identify which third-party software and databases should be cited for a given analysis.

## Downstream analysis

LotuS3 output can be imported into R or other downstream analysis tools. The abundance tables, taxonomy annotations and representative sequences are the main files usually used for ecological, statistical or visualization workflows.
