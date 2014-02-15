#!/bin/bash

# CHUNK 1: In line edit for gene normalized results
sed -i '' '/?/d' *rsem.genes.normalized_results |sed -i '' 's/|[ 0-9 ]*//g' *rsem.genes.normalized_results
