#!/bin/bash

# CHUNK 1: In line edit for gene normalized results
sed -i -e '/?/d' *rsem.genes.normalized_results |sed -i -e 's/|[ 0-9 ]//g' *rsem.genes.normalized_results
