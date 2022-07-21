#!/usr/bin/env bash

# find all .txt files in this directory
# cat txt file
# change name into .fasta when mk is done
find -L . \
 -type f \
 -name "*.txt" \
| xargs cat \
| sed "s#\$#.fasta.fetcherror#" \
| xargs mk
