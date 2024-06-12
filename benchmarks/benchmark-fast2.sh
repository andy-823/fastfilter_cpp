#!/bin/sh
# run the benchmark multiple times with all important algorithms
# for algorithm ids and other parameters, see
# bulk-insert-and-query.cc
# rnd: rando:set ff=unixm number generators to use
for rnd in `seq -1 -1`; do
  # m: number of entries
  for m in 1024000 2560000 6400000 16000000 40000000 100000000; do
    # test: test id
    for test in `seq 1 10`; do
      # sleep 5;
      ./bulk-insert-and-query.exe ${m} 0,2,3,4,40,42,91,92,116,117,118,119,226,227,228,229;
    done;
  done;
done > benchmark-fast2-results.txt 2>&1