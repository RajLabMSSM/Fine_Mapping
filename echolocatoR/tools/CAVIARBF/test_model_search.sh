# This script compares the model_search result against bimbam. 
# It assumes the bimbam is installed in the directory ../bimbam/v1.0
# Assume the current working directory is the caviarbf directory
cd ../bimbam/v1.0
./bimbam -g input/cohort.txt -p input/pheno.txt -pos input/pos.txt -o pref3 -l 3
printf "\n== running model_search ==\n"
../../caviarbf/model_search -p 1 -i output/pref3.multi.txt -m 68 -o output/pref3.multi.txt
printf "\n== bimbam summary output ==\n"
head -n 14 output/pref3.summary.txt
printf "!! Now compare the posterior probability of each model size between model_search and bimbam\n"
printf "They should be similar but not exactly the same because\n"
printf  "1) the null model is included in model_search; 2) BF output precision from bimbam\n"
cd ../../caviarbf
