#!/bin/bash

HOSTFLDS="name,aliases,net,interface,"
CPUFLDS=$(lscpu | cut -d":" -f1 | sed 's/ //g' | awk '{ printf "%s,", $1 }')
RELFLDS="release,build,"
MEMFLDS1=$(free -m -l -w -t | head -1 |  sed 's/^ *//g' | tr -s ' ' '\n')
MEMFLDS2=$(free -m -l -w -t | tail -5 | cut -d":" -f1 | tr '[:upper:]' '[:lower:]')
MEMFLDS=$(for row in $MEMFLDS2; do for col in $MEMFLDS1; do echo -n "${row}_${col},"; done; done)

echo -n "$HOSTFLDS"
echo -n "$CPUFLDS"
echo -n "$RELFLDS"
echo "$MEMFLDS"

HOSTDATA=$(hostname; hostname -A; hostname -a; hostname -I);
CPUDATA=$(lscpu | cut -d":" -f2 | sed 's/^ *//g' | tr ',' ';' | tr -s ' ' ';' | tr '\n' ',')
RELDATA=$(uname -a | tr ' ' ';')
MEMDATA=$(free -m -l -w -t | tail -5 | cut -d":" -f2 | sed 's/^ *//g' | awk '{ printf "%03d,%03d,%03d,%03d,%03d,%03d,%03d\r\n", $1, $2, $3, $4, $5, $6, $7 }')

echo -ne "$HOSTDATA" | awk '{ printf "%s,", $1 }'
echo -ne "$CPUDATA"
echo -ne "$RELDATA"
echo "$MEMDATA"

echo -e "\r\n"


