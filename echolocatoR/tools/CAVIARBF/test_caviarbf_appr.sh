./caviarbf -z ./example/myfile.Z -r ./example/myfile.LD -t 0 -a 0.1281429 -n 2000 -c 5 -o ./example/myfile.sigma0.1281429.bf --appr &&
diff ./example/myfile.sigma0.1281429.bf ./example/myfile.sigma0.1281429.bf.gold &&
rm ./example/myfile.sigma0.1281429.bf &&
./caviarbf -z ./example/myfile.Z -r ./example/myfile.LD -t 1 -a 0.03 -n 2000 -c 5 -o ./example/myfile.pve0.03.bf --appr &&
diff ./example/myfile.pve0.03.bf ./example/myfile.pve0.03.bf.gold &&
rm ./example/myfile.pve0.03.bf &&
echo "PASS (BF)"
