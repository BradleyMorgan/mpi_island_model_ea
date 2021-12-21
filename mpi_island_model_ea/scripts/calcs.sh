#!/bin/bash

MSIZE=$(( $1 * $1 ))
PSZ=$(bc <<< "$MSIZE * $2 * 2")
BCHSZ=$(( $1 * 2 ))
SCAP=$(bc <<< "scale=4; $3 / 100")
RCAP=$(bc <<< "scale=4; $4 / 100")
SMAX=$(bc <<< "$1 * $SCAP")
RMAX=$(bc <<< "$1 * $RCAP")
MAXSZ=$(bc <<< "($1 - $SMAX) * ($1 - $RMAX)")
#EXPSZ=$(bc <<< "$MSIZE * $2 * ($SCAP + $RCAP) * 2")
EXPSZ=$(bc <<< "$MAXSZ * $2 * 2")

echo "communicator size: $1"
echo "send cap: $SCAP"
echo "recv cap: $RCAP"
echo "max send channels: $SMAX"
echo "max recv channels: $SMAX"
echo "matrix dimensons: $MSIZE"
echo "modified dimensions: $MAXSZ"
echo "gene probability: $2"
echo "potential density: $PSZ"
echo "expected density: $EXPSZ"
echo "expected bench density: $BCHSZ"

read -p "Write these values to config.txt? " -n 1 -r
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    exit 1
fi

sed -i "s/^sparsity:.*/sparsity:$2/g" config.txt | grep "sparsity:"
sed -i "s/^send_cap:.*/send_cap:$3/g" config.txt | grep "send_cap:"
sed -i "s/^recv_cap:.*/recv_cap:$4/g" config.txt | grep "recv_cap:"
