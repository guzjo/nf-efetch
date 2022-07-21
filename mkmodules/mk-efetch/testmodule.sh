#!/bin/bash

rm -rf test/results/
mkdir -p test/results/

bash runmk.sh \
&& mv *.faa \
	*.fna \
	*.fasta \
	*.fetcherror \
	test/results \
&& echo "=== Prueba de modulo exitosa ==="
