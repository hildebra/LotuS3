# PacBio CCS/HiFi amplicon processing

LotuS3 can process PacBio CCS/HiFi amplicon reads. For PacBio data, start by setting:

```bash
-p PacBio
```

This activates default settings suitable for many PacBio amplicon workflows.

## Clustering choice

CD-HIT is often a useful clustering choice for PacBio CCS/HiFi amplicons because it is less dependent on read-quality models and read-length uniformity. Other clustering methods can also be used, depending on the data and analysis goal.

## Read-filtering configuration

Useful starting configuration files are:

```text
configs/sdm_PacBio_ITS.txt
configs/sdm_PacBio_LSSU.txt
```

Copy one of these files and adapt it to the expected amplicon length and primer design of your experiment.

Important parameters include:

```text
minSeqLength	700
maxSeqLength	2000
TruncateSequenceLength	-1
```

`minSeqLength` and `maxSeqLength` should reflect the expected amplicon length distribution. For full-length 16S, a range such as 1300-1600 bp may be more appropriate than a broader default range. For 18S amplicons, 1600-2000 bp may be appropriate depending on the amplified region.

`TruncateSequenceLength -1` disables sequence truncation. This is usually appropriate for PacBio CCS/HiFi reads because the full amplicon is expected to be sequenced and read quality does not show the same length-dependent decay as typical Illumina amplicon reads.

## Primer checks

The PacBio configuration files also use strict primer checks:

```text
ExtensivePrimerChecks	T
RejectSeqWithoutFwdPrim	T
RejectSeqWithoutRevPrim	T
```

These settings help remove faulty amplicons that may arise from long PCRs, suboptimal primers, PCR conditions, sample quality or CCS derivation artefacts. If these settings are active, you must provide the amplicon primers to LotuS3; otherwise, reads lacking primer matches will be removed during quality filtering.

## Example command

```bash
./lotus3 -i /my/PacBioDir/ \
  -s configs/sdm_PacBio_my_copy.txt \
  -m PacBioDir/my_PacBio.map \
  -o /my/PacBio/LotuS3 \
  -p PacBio \
  -id 0.97 \
  -CL cdhit \
  -refDB SLV \
  -forwardPrimer XYZ \
  -reversePrimer XYZ
```
