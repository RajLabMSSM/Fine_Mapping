1000-genomes
============

Annotation of 1000 Genomes variants

Joe Pickrell

Each file has all of the variants called from ALL.[chr].integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz. These files were downloaded from:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/

For each variant, we list the chromosome, position (hg19 coordinates), and rs number, followed by the annotations. "1" means the variant falls in the annotation, and "0" means it does not. All annotations are described in the Supplementary Material of Pickrell (2014) http://arxiv.org/abs/1311.4843. The "tssdist" annotation is the distance to the nearest transcription start site in the Ensembl gene databases as downloaded from the UCSC genome browser on 3/19/2013.

There are 451 annotations (450 presence/absence annotations and distance to the nearest transcription start site); see Pickrell (2014) for details.

UPDATE (6/17/2014): The updated version fixes a small number of coding SNPs incorrectly marked as non-coding and changes the coordinates of the SNPs to be consistent with 1000 Genomes.
