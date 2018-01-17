#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times

#base="/home/kitung"
base="/Users/sharon/code"
exe_dir="$base/Marginalize-master"

make readin

data="$exe_dir/data"

# default
#./readin -i 20 -j 20 -k 10 -c 1000000000 -b 1.0 -a $data  -o "$data/output.txt" 

#initial try 2
#readin -i 20 -j 20 -k 20 -c 100000000 -b 1.0 -a $data  -o "$data/output.txt" -r 8.0 -m 0.8 -n 70 -p 67 -q 0.475 -s 0.74 -t 0.21 -u 0.101

#initial try 3
# Parameters -r to -u correspond to the initial guess for the MCMC chain
# this set of parameters gives an acceptance ratio of 0.23
./readin -i 80 -j 80 -k 10 -c 10000000000 -b 0.3 -a $data  -o "$data/output.txt" -r 13.0 -m 1.43 -n 107 -p 77 -q 0.493 -s 1.079 -t 0.23 -u 0.217

#initial try 4; ll = -3014
#readin -i 45 -j 45 -k 10 -c 100000000 -b 0.6 -a $data  -o "$data/output.txt" -r 12.489 -m 1.289 -n 85.6 -p 78.9 -q 0.477 -s 0.992 -t 0.36 -u 0.393


#readin -i 20 -j 20 -k 10 -c 100000000 -b 1.0 -a $data  -o "$data/output.txt" -r 10.0 -m 1.289 -n 85.6 -p 78.9 -q 0.477 -s 0.992 -t 0.36 -u 0.393



times
