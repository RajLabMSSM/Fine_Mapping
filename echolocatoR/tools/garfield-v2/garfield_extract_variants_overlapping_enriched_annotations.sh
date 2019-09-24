#!/bin/bash

##################
PRUNETAGSDIR=$1
CLUMPTAGSDIR=$2
ANNOTDIR=$3
PVALDIR=$4
GARFIELD_PREP=$5
GARFIELD_OUTPUT=$6
GARFIELD_significant_annotations=$7
GARFIELD_VARS=$8
PTHRESH=$9
PENRICH=${10}
echo $GARFIELD_significant_annotations $GARFIELD_VARS $PTHRESH $PENRICH
##################

awk -v PTHRESH=$PTHRESH -v PENRICH=$PENRICH '$2==PTHRESH && $4<PENRICH {print $1,$3,$4, $14}' $GARFIELD_OUTPUT > $GARFIELD_significant_annotations

echo -n > $GARFIELD_VARS.tmp
while read -a line
do
  ID=${line[0]}
  awk -v PTHRESH=$PTHRESH -v ID=$ID '$3<PTHRESH && substr($7,ID,1)==1 {print $1,$2, ID}' $GARFIELD_PREP  >> $GARFIELD_VARS.tmp
done < $GARFIELD_significant_annotations

echo -n > $GARFIELD_VARS.tmp.lead.chr
echo -n > $GARFIELD_VARS.tmp.r01.chr
echo -n > $GARFIELD_VARS.tmp.r08.chr
echo -n > $GARFIELD_VARS.tmp.ann.pval.chr

for CHR in {1..22}
do
awk -v CHR=$CHR '$1==CHR {print CHR":"$2, $3}' $GARFIELD_VARS.tmp|sort -n |uniq > $GARFIELD_VARS.tmp.lead.chr.tmp
WC=$(cat $GARFIELD_VARS.tmp.lead.chr.tmp| wc -l)
if [ $WC -gt 0 ]
then
	cat $GARFIELD_VARS.tmp.lead.chr.tmp >> $GARFIELD_VARS.tmp.lead.chr

	awk -v CHR=$CHR 'NR==FNR{a[$1]=$1;next;} {if (CHR":"$1 in a) print $0}' $GARFIELD_VARS.tmp.lead.chr.tmp $PRUNETAGSDIR/chr$CHR | awk -v CHR=$CHR '{a=$2;gsub(",",","CHR":",a); print CHR":"$1, CHR":"a}' > $GARFIELD_VARS.tmp.r01.chr.tmp
	awk -v CHR=$CHR 'NR==FNR{a[$1]=$1;next;} {if (CHR":"$1 in a) print $0}' $GARFIELD_VARS.tmp.lead.chr.tmp $CLUMPTAGSDIR/chr$CHR | awk -v CHR=$CHR '{a=$2;gsub(",",","CHR":",a); print CHR":"$1, CHR":"a}'> $GARFIELD_VARS.tmp.r08.chr.tmp
	cat $GARFIELD_VARS.tmp.r01.chr.tmp >> $GARFIELD_VARS.tmp.r01.chr
	cat $GARFIELD_VARS.tmp.r08.chr.tmp >> $GARFIELD_VARS.tmp.r08.chr

	awk '{print $1}' $GARFIELD_VARS.tmp.lead.chr.tmp | cat - $GARFIELD_VARS.tmp.r01.chr.tmp $GARFIELD_VARS.tmp.r08.chr.tmp | tr " " "\n" | tr "," "\n" |sort| uniq > $GARFIELD_VARS.tmp.all.chr

	awk -v CHR=$CHR 'NR==FNR{arr[$1];next} CHR":"$1 in arr {print CHR":"$0}' $GARFIELD_VARS.tmp.all.chr $ANNOTDIR/chr$CHR | awk -v CHR=$CHR 'NR==FNR{arr[CHR":"$1]=CHR":"$1; all[CHR":"$1]=$2;next} {if ($1 in arr) print $1,all[$1], $2; else print $1,"NA", $2 }' $PVALDIR/chr$CHR - >> $GARFIELD_VARS.tmp.ann.pval.chr

	echo "CHR" $CHR
fi
done  

rm $GARFIELD_VARS.tmp.r08.chr.tmp $GARFIELD_VARS.tmp.r01.chr.tmp $GARFIELD_VARS.tmp.lead.chr.tmp


### then loop over entries and merge information
echo "ID ANNOTATION OR PVAL VAR_INFO INFO_TAGS_USED_TO_ANNOTATE INFO_TAGS_PRUNED_OUT" > $GARFIELD_VARS
while read -a line
do
  ID=${line[0]}
  OR=${line[1]}
  PVAL=${line[2]}
  ANNOT=${line[3]}

  awk -v ID=$ID '$2==ID {print $1}' $GARFIELD_VARS.tmp.lead.chr  > $GARFIELD_VARS.tmp.lead.chr.tmp
  while read -a var
	do
  	VAR=${var[0]}
  	
  	awk -v VAR=$VAR '$1==VAR {print $2}' $GARFIELD_VARS.tmp.r08.chr | tr "," "\n"| awk -v ID=$ID 'NR==FNR{arr[$1];next} $1 in arr {print $1"("$2","substr($3,ID,1)")"}' - $GARFIELD_VARS.tmp.ann.pval.chr > $GARFIELD_VARS.tmp.r08.chr.tmp

  	awk -v VAR=$VAR '$1==VAR {print $2}' $GARFIELD_VARS.tmp.r01.chr | tr "," "\n"| awk -v ID=$ID 'NR==FNR{arr[$1];next} $1 in arr {print $1"("$2","substr($3,ID,1)")"}' - $GARFIELD_VARS.tmp.ann.pval.chr | awk 'NR==FNR {a[$1];next;} {if (!($1 in a)) print $0}' $GARFIELD_VARS.tmp.r08.chr.tmp - > $GARFIELD_VARS.tmp.r01.chr.tmp

  	VARINFO=$(awk -v VAR=$VAR -v ID=$ID '$1==VAR {print $1"("$2","substr($3,ID,1)")"}' $GARFIELD_VARS.tmp.ann.pval.chr)
	TAGSINFO=NA
	if [ $(cat $GARFIELD_VARS.tmp.r08.chr.tmp| wc -l) -gt 0 ]
	then 
		TAGSINFO=$(cat $GARFIELD_VARS.tmp.r08.chr.tmp | tr "\n" "|"| sed 's/.$//')
	fi
	TAGSINFO2=NA
	if [ $(cat $GARFIELD_VARS.tmp.r01.chr.tmp| wc -l) -gt 0 ]
	then 
		TAGSINFO2=$(cat $GARFIELD_VARS.tmp.r01.chr.tmp | tr "\n" "|"| sed 's/.$//')
	fi

	rm $GARFIELD_VARS.tmp.r08.chr.tmp
  	rm $GARFIELD_VARS.tmp.r01.chr.tmp
  	echo $ID $ANNOT $OR $PVAL $VARINFO $TAGSINFO $TAGSINFO2 >> $GARFIELD_VARS
  	#echo "Variant" $VAR

  done < $GARFIELD_VARS.tmp.lead.chr.tmp
  rm $GARFIELD_VARS.tmp.lead.chr.tmp
  echo "Annotation" $ID $ANNOT
done < $GARFIELD_significant_annotations


