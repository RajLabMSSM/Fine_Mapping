./model_search -i ./example/pref4.multi.txt -o ./example/pref4.multi.txt.prior0 -m 50 -p 0 -sex > temp.console &&
diff ./example/pref4.multi.txt.prior0.exhaustive ./example/pref4.multi.txt.prior0.exhaustive.gold &&
diff ./example/pref4.multi.txt.prior0.exhaustivestepwise ./example/pref4.multi.txt.prior0.exhaustivestepwise.gold &&
diff ./example/pref4.multi.txt.prior0.statistics ./example/pref4.multi.txt.prior0.statistics.gold &&
diff ./example/pref4.multi.txt.prior0.marginal ./example/pref4.multi.txt.prior0.marginal.gold &&
diff ./example/pref4.multi.txt.prior0.stepwise ./example/pref4.multi.txt.prior0.stepwise.gold &&
rm ./example/pref4.multi.txt.prior0.exhaustive &&
rm ./example/pref4.multi.txt.prior0.exhaustivestepwise &&
rm ./example/pref4.multi.txt.prior0.statistics &&
rm ./example/pref4.multi.txt.prior0.marginal &&
rm ./example/pref4.multi.txt.prior0.stepwise && 
rm temp.console && 
echo "PASS (Parameter 0)"