# Author: Zhanghan Tony Ni, December 2022.
# Downloads metadata formatted in JSON
# Work for Professor Rika Anderson for Phosphorus limitation porject.

all_sample=$(cut -f50 matchingMeta_CorrectedDec7.tsv)
for sample_ID in $all_sample; do
	wget https://www.ebi.ac.uk/biosamples/samples/${sample_ID} -O raw_json_v2/${sample_ID}_metadata.json
done