#!/bin/bash

module load openmpi
cp config.txt config.old
git stash
git pull
mpic++ -g -O3 -std=c++11 main.cpp -o ./island
cp config.txt config.pull
cp config.old config.txt

