# e3efficiency

cp -r macro BOLO-01-2021-02-02
cd BOLO-01-2021-02-02

# check that in reco2.sh the path
# /data/$1/data/Efficiency/$2/$1-$2-$3.bin
# is correct

source enveeeroot

# if possible reconstruct a run with triple coinc to get the calibrations
eeereco.exe -b -r -c /data/BOLO-01/data/2021-02-01/BOLO-01-2021-02-01-00059.bin

# use something like this to get all commands you need to run

# ch bottom
cat /data/BOLO-01/data/efficiency/summary_C1.txt|grep C23|awk '{print "./mkpoint.sh BOLO-01 2021-02-02",$5,$3,$4,"0"}'

# ch middle
cat /data/BOLO-01/data/efficiency/summary_C1.txt|grep C13|awk '{print "./mkpoint.sh BOLO-01 2021-02-02",$6,$3,$4,"1"}'

# ch top
cat /data/BOLO-01/data/efficiency/summary_C1.txt|grep C12|awk '{print "./mkpoint.sh BOLO-01 2021-02-02",$7,$3,$4,"2"}'

