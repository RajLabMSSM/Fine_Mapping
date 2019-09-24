## This script can be used to put GWAS summary statistics files into GARFIELD input format
## NOTE: example adta is not supplied and script should be modified so that
chrcol=1 ## the column in GWAS file containing chormosome information
poscol=2 ## the column in GWAS file containing genomic position information
pvalcol=3 ## the column in GWAS file containing GWAS p-value information
TRAITNAME="test" ## name of directory for GWAS trait to be created 
GWASFILENAME="test-GWAS-name" ## name of file containing GWAS summary statistics

# output directory to be used as input for GARFIELD analysis
OUTDIR=../garfield-data/pval/$TRAITNAME
mkdir -p $OUTDIR

for CHR in {1..22} #'X'
do
awk -v chr=$CHR -v chrcol=$chrcol -v poscol=$poscol -v pvalcol=$pvalcol '$chrcol==chr {print $poscol,$pvalcol}' $GWASFILENAME | sort -k1n > $OUTDIR/chr$CHR
echo $chr
done
