#!/bin/bash

module load openmpi
git stash
git pull
mpic++ -g -O3 -std=c++11 ../main.cpp -o ../island

