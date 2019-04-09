set -ue

Rscript test_package.R
for suffix in marginal loglik gamma
do
  for eps in 0 0.2
  do
  	diff ../example/multiLoci/CAVIARBF/priorType0_0.1_exact_eps0_glmnetLASSOMin_cv/00000001_l3.${suffix} ../example/multiLoci/CAVIARBF/priorType0_0.1_exact_eps0_glmnetLASSOMin_cv_gold/00000001_l3.${suffix}
  done
  diff ../example/multiLoci/CAVIARBF/priorType0_0.1_exact_eps0_topK/00000001_l3.${suffix} ../example/multiLoci/CAVIARBF/priorType0_0.1_exact_eps0_topK_gold/00000001_l3.${suffix}
done
echo "PASS"
rm -rf ../example/multiLoci/CAVIARBF/priorType0_0.1_exact_eps0_glmnetLASSOMin_cv
rm -rf ../example/multiLoci/CAVIARBF/priorType0_0.1_exact_eps0.2_glmnetLASSOMin_cv
rm -rf ../example/multiLoci/CAVIARBF/priorType0_0.1_exact_eps0_topK
