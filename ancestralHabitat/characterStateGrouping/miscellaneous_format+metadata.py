# Author: Zhanghan Tony Ni, December 2022.
# 12 miscellaneous code blocks that address issues related to format, metadata, etc.
# Work for Professor Rika Anderson for Phosphorus limitation porject
import csv
import json
import sys
import subprocess
import pandas as pd

csv.field_size_limit(sys.maxsize)

# 1. remove empty line "\n"
# with open('rasp_CB_3.csv', 'w') as out:
#     with open('rasp_CB (2).csv','r') as file:
#         for line in file:
#             if not line.isspace():
#                 out.write(line)

# 2. find samples of which metadata is missing
# df1 = pd.read_csv("missingMeta_corrected.csv")
# df2 = pd.read_csv("mergedv2.csv")
# output = pd.merge(df1, df2,on='accession',how='inner')
# print(output)    
# output.to_csv("missingMetaFull_corrected.csv", index=False)

# 3. find difference between two csv files
# ls = []
# with open(' treePractiplabel.csv', 'r') as t2, open('rownameslocationListPrac.csv', 'r') as t1:
#     fileone = t1.readlines()
#     filetwo = t2.readlines()
#     for line in fileone:
#         ls.append(line)
# with open('matchingmeta_accessiononly_processed.csv', 'w') as out:
#     for i in ls:
#         out.write(i)
# with open('missingPrac.csv', 'w') as outFile:
#     for line in filetwo:
#         if line not in fileone:
#             outFile.write(line)

# 4. change accession number to full tip name
# df1 = pd.read_csv("merged.csv")
# df2 = pd.read_csv("output_fullName.csv")
# for i in range(0, len(df1)):
#     cur = df1.loc[i, 'accession']
#     df1.loc[i, 'accession'] = cur[3:]
# df1.to_csv("mergedv2.csv", index=False)

# 5. modify accession number format to match tip name
# df = pd.read_csv("finalInput.csv")
# for i in range(0, len(df)):
#     cur = str(df.loc[i, "accession"])
#     df.loc[i, "accession"] = cur[:3]+".."+cur[4:]
#     print(cur[:3]+".."+cur[4:])
# df.to_csv("finalInput_processed.csv", index=False)

# 6. get genome names from tree file
# with open("ToL_Final_alignment.txt","r") as f:
#     accNum = []
#     for ln in f:
#         if ln.startswith(">"):
#             i = ln[:17]
#             i = i[1:4]+'_'+i[6:]
#             accNum.append(i)
# print(accNum)
# with open('output_fullName.csv', 'w', newline='') as f:
#     write = csv.writer(f)
#     for i in accNum:
#         # print(i)
#         # i = i[1:4]+'_'+i[6:]
#         # print(i)
#         write.writerow([i])

# 7. get the genome names of which meta data is available
# meta = []
# reader = csv.reader(open('merged.csv'))
# for row in reader:
#     cur = str(row[0])
#     print(cur)
#     if cur[3:] in accNum:
#         meta.append(row)
# with open('matchingMeta_CorrectedDec7.csv', 'w', newline = '') as f:
#     write = csv.writer(f)
#     for i in meta:
#         write.writerow(i)

# 8. generate a list of accession number using which metadata can be downloaded
# json = []
# reader = csv.reader(open('jsonToDownload.csv'))
# for row in reader:
#     print(row[0])
#     if row[0] not in json:
#         json.append(row[0])
# print(json)
# with open('jsonToDownload_Distinct.csv', 'w', newline = '') as f:
#     write = csv.writer(f)
#     for i in json:
#         write.writerow([i])

# 9. convert tsv to csv
# tsv_file='merged.tsv'
# csv_table=pd.read_table(tsv_file,sep='\t')
# csv_table.to_csv('merged.csv',index=False)
# COMMAND = "awk -F '	' ' FNR==NR{idx[$1]; next} FNR==1 || $1 in idx' output.csv merged.tsv"  
# subprocess.call(COMMAND, shell=True)

# 10. get matching genomes from the master sheet
# with open('ar53_metadata_r207.tsv') as csvfile:
#     reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
#     result = [row[0].strip() for row in reader ]
# print(result[1])

# 11. create a dictionary in which key is bioproject number and value is accession number
# meta = []
# reader = csv.reader(open('merged.csv'))
# for row in reader:
#     if row[54] in accNum:
#         meta.append(row)
# with open('dict[bioproject number: accession number].csv', 'w', newline = '') as f:
#     write = csv.writer(f)
#     for i in meta:
#         write.writerow(i)

# 12. convert csv to tsv
# with open('matchingMeta_CorrectedDec7.csv','r') as csvin, open('matchingMeta_CorrectedDec7.tsv', 'w') as tsvout:
#     csvin = csv.reader(csvin)
#     tsvout = csv.writer(tsvout, delimiter='\t')
#     for row in csvin:
#         tsvout.writerow(row)