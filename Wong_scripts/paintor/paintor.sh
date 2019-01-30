# We want to run paintor on each of the PD loci.
# (eventually we may want to narrow this down to only the significant ones, but whatever.)

# mkdir ~/ad-omics/wongg05/paintor
paintorDataDir=~/ad-omics/wongg05/paintor

###############
# Get the list of events: where the PD GWAS variant has a p < 0.05 in the eQTL or sQTL data.
###############

# The a1/a2 info has columns: chrchr pos chrpos a1 a2.
# The qtl data has columns:   event chrpos dist npval slope se.

# First, we add A1/A2 info to the rosmap sQTL z-scores:
for chr in {1..22}; do
	echo $chr
	python ~/tools/myMerge.py ~/ad-omics/splicing/rosmap/sQTL/vcf/chr$chr.allele.tsv /hpc/users/rajt01/ad-omics/splicing/rosmap/sQTL/ROSMAP.sqtl.nominal.ALL.SE.txt 2 1 > $paintorDataDir/rosmap_sQTL_a1a2_chr$chr.txt
done

# We do the same for the rosmap eQTL z-scores: 
for chr in {1..22}; do
	echo $chr
	python ~/tools/myMerge.py ~/ad-omics/splicing/rosmap/sQTL/vcf/chr$chr.allele.tsv ~/ad-omics/splicing/rosmap/eQTL/eQTL.cis1mb.nominal.ALL.txt 2 1 > $paintorDataDir/rosmap_eQTL_a1a2_chr$chr.txt
done

# We check the length of the labelled files, compared to the sQTL file:
wc -l $paintorDataDir/rosmap_sQTL_a1a2_chr*.txt
wc -l /hpc/users/rajt01/ad-omics/splicing/rosmap/sQTL/ROSMAP.sqtl.nominal.ALL.SE.txt

# Got PD GWAS SNPs from Towfique (email, 9 April 2017) at pd_nalls_snps.txt 
python ~/tools/myMerge.py $paintorDataDir/pd_nalls_snps.txt ~/ad-omics/splicing/rosmap/eQTL/HRC.RSID 0 1 > $paintorDataDir/pd_nalls_snps_chrpos.txt
# We see one of them isn't there:
#while read line
#do
#	echo $line
#	grep $line $paintorDataDir/pd_nalls_snps_chrpos.txt | wc -l
#done < $paintorDataDir/pd_nalls_snps.txt
# It's rs8118008. We get its chr:pos from pdgene.org (chr20:3168166) and add manually to $paintorDataDir/pd_nalls_snps_chrpos.txt

# We find the list of events near our PD GWAS SNPs:
while read line; do
	rsid=`echo $line | awk '{print $1}'`
	chrpos=`echo $line | awk '{print $2}'`
	chr=`echo $chrpos | awk 'split($1, a, ":") {print a[1]}'`
	pos=`echo $chrpos | awk 'split($1, a, ":") {print a[2]}'`
	echo $rsid $chrpos $chr $pos
	awk -v pos=$pos '$2 == pos && $9 < 0.05' $paintorDataDir/rosmap_sQTL_a1a2_chr$chr.txt >> $paintorDataDir/pd_nalls_sQTLs_p0.05.txt
done < $paintorDataDir/pd_nalls_snps_chrpos.txt
python ~/tools/myMerge.py $paintorDataDir/pd_nalls_sQTLs_p0.05.txt $paintorDataDir/pd_nalls_snps_chrpos.txt 2 1 > $paintorDataDir/pd_nalls_snps_chrpos_rsid.txt
# chrpos rsid a1 a2 event chrpos dist npval slope se

# We make a locus file for each event in the new list.
while read line; do
	event=`echo $line | awk '{print $6}'`
	chr=`echo $line | awk 'split($3, a, ":") {print a[1]}'`
	echo $event
	echo -e "chr\tpos\tchrpos\ta1\ta2\tevent\tchrpos2\tdist\tnpval\tslope\tse" > $paintorDataDir/"$event".unp
	rocessed
	awk -v event=$event '$6 == event' $paintorDataDir/rosmap_sQTL_a1a2_chr$chr.txt >> $paintorDataDir/"$event".unprocessed
	echo $line > $paintorDataDir/"$event".info
done < $paintorDataDir/pd_nalls_sQTLs_p0.05.txt

#####################
# We've made locus files. Now we just need corresponding LD files and annotations.
####################

# download 1kg so we can make ld files
for chr in 2; do echo "curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -o ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" | bsub -P acc_ad-omics -W 24:00 -R rusage[mem=4000] -q alloc -o download1kg.o -e download1kg.e; done

#make LD files
while read line; do
        event=`echo $line | awk '{print $6}'`
        chr=`echo $line | awk 'split($3, a, ":") {print a[1]}'`
	echo $event
	ldUtility=~/ad-omics/efthymiou/MS4Aproject/programs/PAINTOR_V3.0/PAINTOR_Utilities/CalcLD_1KG_VCF.py
	python $ldUtility -l $paintorDataDir/"$event".unprocessed -r $paintorDataDir/1kg/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -m $paintorDataDir/1kg/integrated_call_samples_v3.20130502.ALL.panel -i pos -e a1 -a a2 -p EUR -z slope -o $paintorDataDir/"$event"
done < $paintorDataDir/pd_nalls_sQTLs_p0.05.txt

while read line; do
	event=`echo $line | awk '{print $6}'`
	thisEventDir=$paintorDataDir/sqtlLoci/$event
	mv $thisEventDir/sQTL.txt $paintorDataDir/sqtlLoci/"$event".unprocessed
	mv $thisEventDir/event.processed $paintorDataDir/sqtlLoci/"$event".processed
	mv $thisEventDir/eventInfo.txt $paintorDataDir/sqtlLoci/"$event".info
	mv $thisEventDir/event.ld $paintorDataDir/sqtlLoci/"$event".ld
	rmdir $thisEventDir
done < $paintorDataDir/pd_nalls_sQTLs_p0.05.txt

# Split up clip-seq bed files....
bedDir=$paintorDataDir/bed
mkdir -p $bedDir
for experiment in `cut -f 4 /hpc/users/rajt01/ad-omics/tools/qtltools/RBP_sqtl/RBPs/CLIP/clip_db_humanRBP.bed | sort -u`; do
	echo $experiment
	awk -v experiment=$experiment '$4 == experiment' /hpc/users/rajt01/ad-omics/tools/qtltools/RBP_sqtl/RBPs/CLIP/clip_db_humanRBP.bed > $bedDir/clipDbHumanRbp_$experiment.bed
done
# Copy over other ones...
cp ~/ad-omics/data/rosmap/xQTL/bed/macs* $bedDir
cp ~/ad-omics/data/rosmap/xQTL/bed/*.broadPeak $bedDir
# Split up E... mnemonic files...
for file in `ls ~/ad-omics/data/rosmap/xQTL/bed/*coreMarks_mnemonics.bed`; do
	echo $file
	fileEnd=`basename $file`
	for mark in `cut -f 4 $file | sort -u`; do
		if [ "$mark" == "8_ZNF/Rpts" ]; then
			mark="8_ZNF_Rpts"
		fi
		awk -v mark=$mark '$4 == mark' $file > $bedDir/"$fileEnd"_$mark.bed
	done
done

# Get list of bed files...
ls $bedDir/* > $paintorDataDir/bedFileManifest.txt

# Make annotation files.
annUtility=~/ad-omics/efthymiou/MS4Aproject/programs/PAINTOR_V3.0/PAINTOR_Utilities/AnnotateLocus.py
while read line; do
	event=`echo $line | awk '{print $6}'`
	chr=`echo $line | awk 'split($3, a, ":") {print a[1]}'`
	echo $event
	python $annUtility -i $paintorDataDir/bedFileManifest.txt -l $paintorDataDir/sqtlLoci/"$event".processed -o $paintorDataDir/sqtlLoci/"$event".ann -c chr -p pos
done < $paintorDataDir/pd_nalls_sQTLs_p0.05.txt


oneTest=7:23293931:23299599:clu_28219

awk '{if (NR==1) {print $0, "z"} else {print $0, $10/$11}}' $paintorDataDir/sqtlLoci/$oneTest > $paintorDataDir/sqtlLoci/$oneTest.2
mv $paintorDataDir/sqtlLoci/$oneTest.ann $paintorDataDir/sqtlLoci/$oneTest.annotations
mv $paintorDataDir/sqtlLoci/$oneTest.ld $paintorDataDir/sqtlLoci/$oneTest.LD
echo $oneTest > $paintorDataDir/inputFile.txt


# run paintor on one locus
paintor=~/ad-omics/efthymiou/MS4Aproject/programs/PAINTOR_V3.0/PAINTOR
$paintor -input inputFile.txt -Zhead z -LDname ld -in $paintorDataDir/sqtlLoci -out $paintorDataDir/sqtlLoci -enumerate 2 -Gname $oneTest.Enrich.Base -Lname $oneTest.BF.Base

$paintor -input inputFile.txt -Zhead z -LDname ld -in $paintorDataDir/sqtlLoci -out $paintorDataDir/sqtlLoci.Base -enumerate 2 -Gname $oneTest.Enrich.Base -Lname $oneTest.BF.Base


ann=E068_15_coreMarks_mnemonics.bed_15_Quies.bed
$paintor -input inputFile.txt -Zhead  z -LDname ld -in $paintorDataDir/sqtlLoci -out $paintorDataDir/sqtlLoci.$ann -enumerate 2 -annotations $ann  -Gname $paintorDataDir/$oneTest.Enrich.$ann -Lname $paintorDataDir/$oneTest.BF.$ann

$paintor -input inputFile.txt -Zhead  z -LDname ld -in $paintorDataDir/sqtlLoci -out test.txt -enumerate 2 -annotations $ann  -Gname $paintorDataDir/$oneTest.Enrich.$ann -Lname $paintorDataDir/$oneTest.BF.$ann



