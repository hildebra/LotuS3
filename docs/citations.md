# Citations

If you use LotuS3, please cite the main pipeline publication and any relevant method-specific tools and databases used in your run.

LotuS3 writes method-specific citations for each run to:

```text
LotuSLogS/citations.txt
```

## Main LotuS3 citations

**Pipeline** - Özkurt E, Fritscher J, et al. (2022) LotuS2: An ultrafast and highly accurate tool for amplicon sequencing analysis. *Microbiome* 10:176. doi:10.1186/s40168-022-01365-1.

**Off-target removal** - Bedarf JR, Beraza N, Khazneh H, Özkurt E, et al. (2021) Much ado about nothing? Off-target amplification can lead to false-positive bacterial brain microbiome detection in healthy and Parkinson's disease individuals. *Microbiome* 9:75.

Related publications:

- LotuS2: https://www.biorxiv.org/content/10.1101/2021.12.24.474111v1
- Off-target removal: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01012-1
- LotuS: http://www.microbiomejournal.com/content/2/1/30

## Software and algorithms used by LotuS3

- **DADA2** - Callahan B, McMurdie PJ, Rosen MJ, et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods* 13:581-583.
- **UPARSE** - Edgar RC. (2013) UPARSE: highly accurate OTU sequences from microbial amplicon reads. *Nature Methods* 10:996-998.
- **VSEARCH** - Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. *PeerJ* 4:e2584.
- **swarm** - Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) Swarm: robust and fast clustering method for amplicon-based studies. *PeerJ* 2:e593.
- **CD-HIT** - Fu L, Niu B, Zhu Z, Wu S, Li W. (2012) CD-HIT: accelerated clustering for next-generation sequencing data. *Bioinformatics* 28:3150-3152.
- **UCHIME** - Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R. (2011) UCHIME improves sensitivity and speed of chimera detection. *Bioinformatics* 27:2194-2200.
- **RDP classifier** - Wang Q, Garrity GM, Tiedje JM, Cole JR. (2007) Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Applied and Environmental Microbiology* 73:5261-5267.
- **lambda aligner** - Hauswedell H, Singer J, Reinert K. (2014) Lambda: the local aligner for massive biological data. *Bioinformatics* 30:i349-i355.
- **BLAST+** - Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. (1990) Basic local alignment search tool. *Journal of Molecular Biology* 215:403-410.
- **Clustal Omega** - Sievers F, Wilm A, Dineen D, Gibson TJ, et al. (2011) Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. *Molecular Systems Biology* 7:539.
- **MAFFT** - Katoh K, Standley DM. (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30:772-780.
- **FastTree 2** - Price MN, Dehal PS, Arkin AP. (2010) FastTree 2: approximately maximum-likelihood trees for large alignments. *PLoS ONE* 5:e9490.
- **IQ-TREE 2** - Nguyen LT, Schmidt HA, von Haeseler A, Minh BQ. (2015) IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. *Molecular Biology and Evolution* 32:268-274.

## Databases

- **SILVA** - Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, et al. (2014) The SILVA and All-species Living Tree Project taxonomic frameworks. *Nucleic Acids Research* 42:D643-D648.
- **Greengenes2** - McDonald D, et al. (2023) Greengenes2 unifies microbial data in a single reference tree. *Nature Biotechnology*. doi:10.1038/s41587-023-01845-1.
- **PR2** - Guillou L, Bachar D, Audic S, et al. (2013) The Protist Ribosomal Reference database: a catalog of unicellular eukaryote small sub-unit rRNA sequences with curated taxonomy. *Nucleic Acids Research* 41:D597-D604.
- **KSGP** - Grant A, Aleidan A, Davies CS, et al. (2023) Improved taxonomic annotation of Archaea communities using LotuS2, the Genome Taxonomy Database and RNAseq data. bioRxiv. https://ksgp.earlham.ac.uk/
- **HITdb** - Ritari J, Salojärvi J, Lahti L, de Vos WM. (2015) Improved taxonomic assignment of human intestinal 16S rRNA sequences by a dedicated reference database. *BMC Genomics* 16:1056.
- **beetax** - Jones JC, Fruciano C, Hildebrand F, et al. (2018) Gut microbiota composition is associated with environmental landscape in honey bees. *Ecology and Evolution* 8:441-451.

## ITS-specific resources

- **UNITE ITS chimera database** - Nilsson et al. (2015) A comprehensive, automatically updated fungal ITS sequence dataset for reference-based chimera control in environmental sequencing efforts. *Microbes and Environments*.
- **UNITE ITS taxonomic reference database** - Kõljalg U, et al. (2013) Towards a unified paradigm for sequence-based identification of fungi. *Molecular Ecology* 22:5271-5277.
- **ITSx** - Bengtsson-Palme J, Ryberg M, Hartmann M, et al. (2013) Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis in environmental sequencing data. *Methods in Ecology and Evolution*.

## Mathematical models and C++ libraries

- Puente-Sánchez F, et al. (2016) A novel conceptual approach to read-filtering in high-throughput amplicon sequencing studies. *Nucleic Acids Research* 44:e40.
- **gzip libraries** - `gzstream.h` and zlib: https://gist.github.com/piti118/1508048 and http://www.zlib.net/
- **robin_hood hash map library** - https://github.com/martinus/robin-hood-hashing
