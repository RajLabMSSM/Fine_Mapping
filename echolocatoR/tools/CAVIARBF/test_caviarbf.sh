set -ue

./caviarbf -z ./example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.z -r ./example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12.LD -t 0 -a 1.6 -n 471 -c 2 -o ./example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12_sigmaa1.6_l2.bf &&
diff ./example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12_sigmaa1.6_l2.bf ./example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12_sigmaa1.6_l2.bf.gold &&
rm ./example/eQTL/region.LOC284581.chr1.205831207.205865215.dosage.p1e-12_sigmaa1.6_l2.bf &&
echo "PASS (BF exact)"
./caviarbf -z ./example/myfile.Z -r ./example/myfile.LD -t 0 -a 0.1281429 -n 2000 -c 5 --appr -o ./example/myfile.sigma0.1281429.bf 
diff ./example/myfile.sigma0.1281429.bf ./example/myfile.sigma0.1281429.bf.gold 
rm ./example/myfile.sigma0.1281429.bf 
./caviarbf -z ./example/myfile.Z -r ./example/myfile.LD -t 1 -a 0.03 -n 2000 -c 5 --appr -o ./example/myfile.pve0.03.bf 
diff ./example/myfile.pve0.03.bf ./example/myfile.pve0.03.bf.gold 
rm ./example/myfile.pve0.03.bf 
echo "PASS (BF approximate)"
# test single causal variant
./caviarbf -z ./example/myfile.Z -r ./example/myfile.LD -t 0 -a 0.1 -n 2000 -c 1 --appr -o ./example/myfile.sigma0.1.c1.bf 
diff ./example/myfile.sigma0.1.c1.bf ./example/myfile.sigma0.1.c1.bf.gold 
# using identity option
./caviarbf -z ./example/myfile.Z -i -t 0 -a 0.1 -n 2000 -c 1 --appr -o ./example/myfile.sigma0.1.i.c1.bf
diff ./example/myfile.sigma0.1.c1.bf ./example/myfile.sigma0.1.c1.bf.gold
rm ./example/myfile.sigma0.1.c1.bf
rm ./example/myfile.sigma0.1.i.c1.bf
echo "PASS (BF identity matrix)"
