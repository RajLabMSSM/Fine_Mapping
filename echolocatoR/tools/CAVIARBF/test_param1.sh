./model_search -i ./example/pref4.multi.txt -o ./example/pref4.multi.txt.prior1 -m 50 -p 1 -sex > temp.console &&
diff ./example/pref4.multi.txt.prior1.exhaustive ./example/pref4.multi.txt.prior1.exhaustive.gold &&
diff ./example/pref4.multi.txt.prior1.exhaustivestepwise ./example/pref4.multi.txt.prior1.exhaustivestepwise.gold &&
diff ./example/pref4.multi.txt.prior1.statistics ./example/pref4.multi.txt.prior1.statistics.gold &&
diff ./example/pref4.multi.txt.prior1.marginal ./example/pref4.multi.txt.prior1.marginal.gold &&
diff ./example/pref4.multi.txt.prior1.stepwise ./example/pref4.multi.txt.prior1.stepwise.gold &&
rm ./example/pref4.multi.txt.prior1.exhaustive &&
rm ./example/pref4.multi.txt.prior1.exhaustivestepwise &&
rm ./example/pref4.multi.txt.prior1.statistics &&
rm ./example/pref4.multi.txt.prior1.marginal &&
rm ./example/pref4.multi.txt.prior1.stepwise && 
rm temp.console && 
echo "PASS (Parameter 1)"