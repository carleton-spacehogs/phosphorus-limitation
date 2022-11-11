# Author: Jimmy (Juntao) Zhong, 06/22/2022.
# Work for Professor Rika Anderson for Phosphorus limitation porject
import xml.etree.ElementTree as ET
import csv

def get_metadata_row(sample_ID):
	sample_tree = ET.parse(f"raw_xml/{sample_ID}_metadata.xml").getroot()
	all_attri = sample_tree.findall('.//SAMPLE_ATTRIBUTE')
	this_sample = [None] * 7
	for attri in all_attri:
		data_type = attri[0].text
		value = attri[1].text
		# I want [collection_date, geo_loc_name, lat_lon, temp, phosphate, nitro, ENA-BASE-COUNT]
		if data_type == "collection_date": this_sample[0] = value
		elif data_type == "geo_loc_name": this_sample[1] = value
		elif data_type == "lat_lon": this_sample[2] = value
		elif data_type == "temp": this_sample[3] = value.replace(" C", "")
		elif data_type == "phosphate": this_sample[4] = value.replace(" umol/L","")
		elif data_type == "nitro": this_sample[5] = value.replace(" umol/L","")
		elif data_type == "ENA-BASE-COUNT": this_sample[6] = value
	return this_sample


with open("Bio-Go-links-download.tsv") as f:
	all_samples = [row.split() for row in f]

all_samples = all_samples[1:] # emit the colnames
out_metadata = []
out_metadata.append(["SRA_Accession","sample_accession","Sample","collection_date",
"geo_loc_name","lat_lon","temperature","phosphate","nitro","ENA-BASE-COUNT"])

existed = set()
for sample in all_samples:
	sample_ID, run_ID, sample_tile = sample[0], sample[1], sample[6]
	
	if sample_ID in existed: # there is a large fastq file and a small one for almost all samples.
		continue
	else:
		existed.add(sample_ID)
	base_info = [sample_ID, run_ID, sample_tile]
	xml_info = get_metadata_row(sample_ID)
	out_metadata.append(base_info + xml_info)

with open("Bio-GO-SHIP_sample_metadata.csv","w") as my_csv:
	csvWriter = csv.writer(my_csv,delimiter=',')
	csvWriter.writerows(out_metadata)
