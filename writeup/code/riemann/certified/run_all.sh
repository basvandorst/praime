P=${PYTHON:-python3}
D=$(dirname "$0")
ulimit -v 2500000
echo "== existing formulas"; timeout 900 nice $P $D/existing_formulas.py; timeout 60 $P $D/rowland.py
echo "== local criterion, page zeros, eps 1e-5"; EPS=1e-5 timeout 1500 nice $P $D/local_crit.py best
EPS=1e-5 timeout 900 nice $P $D/local_crit.py range 10001 10401 2 8 1.1
EPS=1e-5 timeout 1500 nice $P $D/local_crit.py range 100001 100201 2 9 1.1
echo "== window decoding, page zeros eps 1e-5 against 28-digit zeros eps 1e-12"
for e in "1e-5 z2m.f64" "1e-12 hp2419.npy"; do set -- $e
 for cfg in "12 12 2 1.2 6 1 1001 1401" "60 14 2 1.2 30 1 1001 1401" "120 16 2 1.2 60 1 1001 1401" "150 16 2 1.2 75 1 1001 1401" "80 14 2 1.2 40 1 9001 9301" "100 15 2 1.2 50 1 9001 9301" "130 16 2 1.2 65 1 9001 9301"; do
  echo -n "eps=$1: "; ZFILE=$PWD/$2 EPS=$1 timeout 900 nice $P $D/window_cert.py $cfg truth 2>&1 | tail -1; done; done
echo "== error budget"; timeout 600 nice $P $D/budget.py
echo "== section 14 of RIEMANN.md: certified chain from 3 at height 782, independent check, stretch 9001-9200"
timeout 7200 nice $P $D/lpchain.py 10000 '[[20,1.4,6,2,1.2,1],[100,1.9,6,2,1.2,1],[300,6,10,3,1.2,3],[1000,24,14,3,1.2,12],[3000,84,18,3,1.2,42],[6000,221,24,3,1.2,110],[10000,429,28,3,1.2,214]]' '{"maxbits":40,"limit":1500,"record":[9001,9201],"out":"lpchain785.json"}' | tail -2
timeout 7200 nice $P $D/lp_check.py lpchain785.json 2
timeout 14400 nice $P $D/lpchain.py 10000 '[[20,1.4,6,2,1.2,1],[100,1.9,6,2,1.2,1],[300,6,10,3,1.2,3],[1000,27,14,3,1.2,13],[3000,104,18,3,1.2,52],[6000,276,24,3,1.2,138],[8000,429,28,3,1.2,214],[10000,612,32,3,1.2,306]]' '{"maxbits":40,"limit":4000,"record":[9001,9061],"out":"lpchain628.json"}' | tail -2
timeout 7200 nice $P $D/lp_check.py lpchain628.json 1
for cfg in "357,28,3,1.2,178" "428,28,3,1.2,214" "611,32,3,1.2,305"; do
 timeout 3600 nice $P $D/lpchain.py 9201 "[[9201,$cfg]]" '{"truth_below":9001,"start":9001,"maxbits":40,"limit":1500}' | tail -2; done
echo "== section 15 of RIEMANN.md: the family criterion mod 20011 (zeros of 20,009 L-functions below 13.211), independent check"
timeout 3000 nice $P $D/lfam.py 20011 13.211009174311924 0.005 lfam20011.npz | tail -5
timeout 3000 nice $P $D/famcrit.py lfam20011.npz 10000 1.09 12 4 1.2 1e-9 famcrit20011.txt | tail -6
timeout 5400 nice $P $D/fam_check.py lfam20011.npz famcrit20011.npy 1.09 12 4
echo "== sections 15-17: family mod 395137 (n < 200,000), class peeling mod 30000001, hard decisions of the chain"
timeout 10800 nice $P $D/lfam2.py 395137 15.412844036697248 0.005 lfam395137.npz | tail -3
timeout 7200 nice $P $D/lfam_fix.py lfam395137.npz lfam395137f.npz | tail -9
timeout 7200 nice $P $D/famcrit2.py lfam395137f.npz 200000 1.09 14 6 1.2 1e-9 famcrit395137 | grep -v "grid point"
timeout 7200 nice $P $D/charlog_peel.py 30000001 3
timeout 7200 nice $P $D/lpchain.py 10000 '[[20,1.4,6,2,1.2,1],[100,1.9,6,2,1.2,1],[300,6,10,3,1.2,3],[1000,24,14,3,1.2,12],[3000,84,18,3,1.2,42],[6000,221,24,3,1.2,110],[10000,429,28,3,1.2,214]]' '{"maxbits":40,"limit":1500,"out":"lpchain785b.json"}' | tail -2
timeout 600 nice $P $D/miss_analysis.py lpchain785b_decisions.json
echo "== section 19: the localization law (families mod 101, 1009, 10007, 20011; every n < 10,000 alone)"
for spec in "101 1690 0.03" "1009 177 0.01" "10007 24.5 0.005" "20011 15.5 0.005"; do set -- $spec; timeout 5400 nice $P $D/lfam3.py $1 $2 $3 l3_$1.npz | tail -2; done
timeout 7200 nice $P $D/famcrit3.py l3_101.npz 10000 0.009950164 4 1.2 4e-6 3,4,6,8,10,12,14 fc3_101.json
timeout 7200 nice $P $D/famcrit3.py l3_1009.npz 10000 0.0951758 4 1.2 1e-8 3,4,6,8,10,12 fc3_1009.json
timeout 7200 nice $P $D/famcrit3.py l3_10007.npz 10000 0.6866116 4 1.2 1e-10 4,6,8,10,12,14 fc3_10007.json
timeout 7200 nice $P $D/famcrit3.py l3_20011.npz 10000 1.09 4 1.2 1e-9 3,4,6,8,10,12 fc3_20011.json
$P $D/uncertainty_table.py fc3_101.json fc3_1009.json fc3_10007.json fc3_20011.json
echo "== section 20: the information law (digits sweep, window law, chain from 3 near 3.2 sqrt(x), information limits in mpmath)"
for spec in "3 0.15 0.25 4 3001" "5 0.06 0.08 4 9001" "6 0.035 0.06 4 9001" "8 0.027 0.036 4 9001" "10 0.022 0.03 4 9001"; do timeout 9000 nice $P $D/info_law.py bisect-digits $spec | tail -1; done
timeout 3600 nice $P $D/info_law.py bisect-window 101 40 0.5 0.06 0.1 4 | tail -1
timeout 3600 nice $P $D/info_law.py bisect-window 301 60 0.5 0.033 0.05 4 | tail -1
timeout 900 nice $P $D/info_law.py digits 5 0.07182 9001 1500 amb5.json > /dev/null; timeout 1800 nice $P $D/info_law.py check-ambiguous 5 amb5_ambiguous.json
timeout 3600 nice $P $D/lpchain.py 905 "$($P $D/info_law.py stages 1200 3.2 0.5 40)" '{"maxbits":40,"limit":3000,"record":[860,905],"out":"chain905.json"}' | tail -2
timeout 3600 nice $P $D/lp_check.py chain905.json 2
echo done
