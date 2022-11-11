#!/bin/bash

# awk 'NR!=1{{print $5} {print $6}}' Bio-Go-links-download.tsv > download_links_all_samples.txt

start=11
end=20
lines=$(sed -n -e ${start},${end}p /localdata/researchdrive/zhongj2/P-limitation/code-phosphate-limitation/organize_metadata/download_links_all_samples.txt)
for line in $lines; do wget $line; donex