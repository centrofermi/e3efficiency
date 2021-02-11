#!/bin/bash
# BOLO-01 2016-04-22 15600 29 31 0

rm $1-$2-$3

#reconstruct runs
for i in `seq -f "%05g" $4 $5`;
do
./reco2.sh $1 $2 $i $6
done

ls $1-$2-*_dst.root >> $1-$2-$3

eeeroot.exe -b -q -l DoMerge.C\(\"$1-$2-$3\",\"$1-$2-$3.root\"\)

for i in `seq -f "%05g" $4 $5`;
do
rm $1-$2-$i.root.eee
done

rm *_dst.root


root -b -q -l ProcessRun.C\(\"$1-$2-$3.root\",$6,$3\)

