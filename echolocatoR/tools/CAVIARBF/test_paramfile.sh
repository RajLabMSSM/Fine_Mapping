./model_search -i ./example/pref4.multi.txt -o ./example/pref4.multi.txt.priorFile -m 50 -f ./example/pref4.SNPPrior -sex > temp.console &&
diff ./example/pref4.multi.txt.priorFile.exhaustive ./example/pref4.multi.txt.priorFile.exhaustive.gold &&
diff ./example/pref4.multi.txt.priorFile.exhaustivestepwise ./example/pref4.multi.txt.priorFile.exhaustivestepwise.gold &&
diff ./example/pref4.multi.txt.priorFile.statistics ./example/pref4.multi.txt.priorFile.statistics.gold &&
diff ./example/pref4.multi.txt.priorFile.marginal ./example/pref4.multi.txt.prior0.marginal.gold &&
diff ./example/pref4.multi.txt.priorFile.stepwise ./example/pref4.multi.txt.prior0.stepwise.gold &&
rm ./example/pref4.multi.txt.priorFile.exhaustive &&
rm ./example/pref4.multi.txt.priorFile.exhaustivestepwise &&
rm ./example/pref4.multi.txt.priorFile.statistics &&
rm ./example/pref4.multi.txt.priorFile.marginal &&
rm ./example/pref4.multi.txt.priorFile.stepwise && 
rm temp.console && 
echo "PASS (Parameters from File)"