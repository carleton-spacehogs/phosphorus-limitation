#!/bin/bash

all_sample=$(cut -f1 Bio-Go-metagenome-metadata.tsv)
for sample_ID in $all_sample; do
	wget https://www.ebi.ac.uk/ena/browser/api/xml/${sample_ID} -O raw_xml/${sample_ID}_metadata.xml
done
