#!/bin/bash

rg1=( $( seq 128 128 512 | shuf ) )
rg2=( $( seq .08 .02 .32 | shuf ) )
rg3=( $( seq 20 2 40 | shuf ) )
rg4=( $( seq 6 1 16 | shuf ) )
rg5=( $( seq .01 .008 .09 | shuf ) )
rg6=( $( seq .03 .01 .20 | shuf ) )
rg7=( $( seq .10 .1 .32 | shuf ) )

echo "l1=$rg1 s=$rg2 m2=$rg3 l2=$rg4 is=$rg5 ir=$rg6 m=$rg7"

sed -i "s/^ea_1_lambda:.*/ea_1_lambda:$rg1/g" config.txt | grep "ea_1_lambda:"
sed -i "s/^ea_2_mu:.*/ea_2_mu:$rg3/g" config.txt | grep "ea_2_mu:"
sed -i "s/^ea_2_lambda:.*/ea_2_lambda:$rg4/g" config.txt | grep "ea_2_lambda:"
sed -i "s/^sparsity:.*/sparsity:$rg2/g" config.txt | grep "sparsity:"
sed -i "s/^send_cap:.*/send_cap:$rg5/g" config.txt | grep "send_cap:"
sed -i "s/^recv_cap:.*/recv_cap:$rg6/g" config.txt | grep "recv_cap:"
sed -i "s/^ea_1_mutation_rate:.*/ea_1_mutation_rate:$rg6/g" config.txt | grep "ea_1_mutation_rate:"
sed -i "s/^ea_2_mutation_rate:.*/ea_2_mutation_rate:$rg7/g" config.txt | grep "ea_2_mutation_rate:"

grep "ea_1_lambda:" config.txt
grep "ea_2_mu:" config.txt
grep "ea_2_lambda:" config.txt
grep "sparsity:" config.txt
grep "send_cap:" config.txt
grep "recv_cap:" config.txt
grep "ea_1_mutation_rate:" config.txt
grep "ea_2_mutation_rate:" config.txt

#sed -i "s/^ea_1_mu:.*/ea_1_mu:$rg2/g" config.txt | grep "ea_1_mu:"
