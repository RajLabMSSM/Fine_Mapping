

# We run coloc on the PD GWAS and best eQTLs.

while read line; do
	leadSnp=`echo $line | awk '{print $1}'`
	bestGene=`echo $line | awk '{print $10}'`
	echo $leadSnp
	chr=`echo $line | awk 'split($8, a, ":") {print a[1]}'`
	chrpos=`echo $line | awk '{print $8}'`
	#grep $bestGene ../overlapsAndQq/rosmap_genes_$leadSnp.txt | awk '{print $1, "chr"$2, $3, $4, $5}' > rosmap_genes_$leadSnp"_"$bestGene.txt
	#python ~/tools/myMerge.py rosmap_genes_$leadSnp"_"$bestGene.txt /hpc/users/rajt01/ad-omics/data/PD_GWAS/META.PD.NALLS2014.chrPos.rsid.txt 1 0 > rosmap_genes_$leadSnp"_"$bestGene"_"gwas.txt
	#python ~/tools/myMerge.py rosmap_genes_$leadSnp"_"$bestGene"_"gwas.txt  rosmap.chr$chr.toMerge.frq 1 1 > rosmap_genes_$leadSnp"_"$bestGene"_"gwas_maf.txt
	wc -l rosmap_genes_$leadSnp"_"$bestGene"_"gwas.txt rosmap_genes_$leadSnp"_"$bestGene"_"gwas_maf.txt
	
done < ../overlapsAndQq/rosmap_all.txt


plink --vcf /hpc/users/rajt01/rosmap/genotypes/imputation_HRC/vcf/chr1.vcf.gz --freq --out rosmap.chr1
for i in {1..22}; do awk 'NR > 1 {print $1, "chr"$2, $3, $4, $5, $6}' rosmap.chr$i.frq > rosmap.chr$i.toMerge.frq; done

