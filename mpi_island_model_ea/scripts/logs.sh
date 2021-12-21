#!/bin/bash

echo "SOLVER"

for dir in logs/$1/*/*/; do
echo $dir
printf "%1s %1s %16s %3s %4s\r\n" r c avg_sol_fit t chan
tail -20 $dir/stats/meta_sol_${1}_*.csv | cut -d"," -f1,2,4,9,11
done

echo "META"
for dir in logs/$1/*/*/; do
echo $dir
printf "%1s %1s %16s %3s %4s\r\n" r c avg_sol_fit t chan
tail -20 $dir/stats/meta_topo_${1}_*.csv | cut -d"," -f1,2,3,4,5,7
done

echo "RUN $1 VALUE COMPARISON"

grep -E "^1,$2" logs/$1/*/*/*/meta_topo* | tail -10 | cut -d"," -f1,14,16,17

