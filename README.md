# gnomAD Group Ancestry Estimation Shiny App

### Purpose
Estimate the proportion of reference ancestry groups in summary genotype frequency data.

### Data
Our reference panel was created from 1000 Genomes Project (GRCh37/hg19) superpopulations (African, Non-Finish European, East Asian, South Asian) and an Indigenous American population (616,568 SNPs and 43 individuals, GRCh37/hg19). Tri-allelic SNPs and SNPs with missing allele frequency information were removed, leaving 613,298 SNPs across the 22 autosomes.

We estimate the ancestry proportions from gnomAD V2 (GRCh37/hg19). After merging with our reference panel we checked for allele matching and strand flips. Our final dataset had 582,550 genome SNPs and 9,835 exome SNPs across the 22 autosomes.
