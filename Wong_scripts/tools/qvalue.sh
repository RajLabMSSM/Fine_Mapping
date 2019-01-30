# $1 a file with pvalues
# $2 column of pvalues
# $3 TRUE if header
cut -f $2 $1 > $1.pvalue.tmp
if $3; then
	echo "qvalues" > $1.qvalue.tmp
	Rscript ~/tools/qvalue.R $1.pvalue.tmp TRUE >> $1.qvalue.tmp
else
	Rscript ~/tools/qvalue.R $1.pvalue.tmp FALSE > $1.qvalue.tmp
fi
paste $1 $1.qvalue.tmp
rm $1.pvalue.tmp
rm $1.qvalue.tmp

