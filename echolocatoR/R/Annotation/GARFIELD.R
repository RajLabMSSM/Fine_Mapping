
"""
Pre-computed LD (European samples - UK10K sequence data), MAF, TSS distance, p-value files for two example traits (Crohn’s Disease from the IBD Consortium and Height from the GIANT consortium) and annotation files for 1005 GENCODE, ENCODE and Roadmap Epigenomics an- notations can be downloaded from http://www.ebi.ac.uk/birney-srv/GARFIELD/package/garfield- data.tar.gz. Note the data is 5.9Gb in compressed format and needs to be uncompressed prior to analysis (83Gb). Variant genomic position (build 37) is used as an identifier in all data files.

Documentation:
https://bioconductor.org/packages/release/bioc/manuals/garfield/man/garfield.pdf
Turorial:
https://www.ebi.ac.uk/birney-srv/GARFIELD/
"""
# Load in Chimera
system("ml garfield")
# Load in R
library(garfield) # BiocManager::install("garfield")

# garfield.run()
 
# INPUT FILES
 
