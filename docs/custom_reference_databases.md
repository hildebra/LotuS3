# Custom reference databases

LotuS3 can use custom reference databases for taxonomic annotation. A custom database requires both a FASTA-formatted reference sequence file and a tab-delimited taxonomy file.

## Required flags

Use:

- `-refDB` for the FASTA reference database;
- `-tax4refDB` for the taxonomy file.

The format is the same as for databases installed by the LotuS3 autoinstaller. Useful examples include:

```text
DB/SLV_138_SSU.fasta
DB/SLV_138_LSU.tax
```

## Taxonomy format

The taxonomy file uses seven fixed levels:

1. kingdom
2. phylum
3. class
4. order
5. family
6. genus
7. species

Levels are denoted by tags such as `k__`, `p__`, `c__`, and are separated by semicolons. If taxonomy information is missing, use `?`.

Example taxonomy line:

```text
FJ588878	k__Eukaryota; p__Phragmoplastophyta; c__?; o__?; f__?; g__?; s__Osyris wightiana
```

## Example: using SILVA138 as a custom database

```bash
./lotus3 -tax4refDB DB/SLV_138_SSU.tax \
  -refDB DB/SLV_138_SSU.fasta \
  -i Example/ \
  -m Example/miSeqMap.sm.txt \
  -o myTestRun_customDB \
  -forwardPrimer GTGYCAGCMGCCGCGGTAA \
  -reversePrimer GGACTACNVGGGTWTCTAAT \
  -CL uparse \
  -taxAligner vsearch
```

This example uses SILVA138 as a custom database, VSEARCH as the taxonomic search algorithm and UPARSE for OTU clustering.

## Multiple reference databases

Multiple complementary reference databases can be searched by providing comma-separated FASTA and taxonomy files:

```bash
./lotus3 -tax4refDB DB/SLV_138_SSU.tax,DB/HITdb/HITdb_taxonomy.txt \
  -refDB DB/SLV_138_SSU.fasta,DB/HITdb/HITdb_sequences.fna \
  -i Example/ \
  -m Example/miSeqMap.sm.txt \
  -o myTestRun_multiDB \
  -forwardPrimer GTGYCAGCMGCCGCGGTAA \
  -reversePrimer GGACTACNVGGGTWTCTAAT \
  -CL uparse \
  -taxAligner vsearch
```

For built-in databases, the shorter syntax can be used:

```bash
./lotus3 -refDB SLV,HITdb \
  -i Example/ \
  -m Example/miSeqMap.sm.txt \
  -o myTestRun_multiDB \
  -forwardPrimer GTGYCAGCMGCCGCGGTAA \
  -reversePrimer GGACTACNVGGGTWTCTAAT \
  -CL uparse \
  -taxAligner vsearch
```

The order of databases can affect results. For example, `-refDB GG2,SLV` and `-refDB SLV,GG2` can differ because the first database is treated as the primary annotation source.
