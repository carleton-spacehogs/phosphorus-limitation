# Author: Zhanghan Tony Ni, December 2022.
# Associates character states with tree tips for ancestral state reconstruction. Deals with special cases.
# Work for Professor Rika Anderson for Phosphorus limitation porject.

from doctest import OutputChecker
import pandas as pd

#connect states with ncbi accession number
states = pd.read_csv('locations_grouped_1214_v2.csv')
tips = pd.read_csv('matchingMeta_CorrectedDec7.csv')
merged = pd.merge(states, tips, on='ncbi_biosample', how='inner')

#remove substring for formatting
for i in range(0, len(merged)):
    cur = merged.loc[i, 'accession']
    merged.loc[i, 'accession'] = cur[3:]
    
#find tips of which geographic info is missing
missing = pd.read_csv("missingMeta_corrected.csv")
fullName = pd.read_csv('output_bothNames_v2.csv')
missing_fullName = pd.merge(fullName, missing, on='accession', how='inner')
missing_fullName.to_csv("missingMeta_fullName.csv")

#assign character state for uncategorized tips
for i in range(0, len(missing_fullName)):
    missing_fullName.loc[i, 'location_grouped'] = "ungrouped"

#connect states with tip names
loc_fullName = pd.merge(fullName, merged, on='accession', how='inner')

#change format to match tip names on tree
def changeFormat(df, col):
    for i in range(0, len(df)):
        cur = str(df.loc[i, col])
        df.loc[i, col] = cur[:3]+".."+cur[4:]
changeFormat(loc_fullName, "num")
changeFormat(missing_fullName, "num")

#write to file
output = [loc_fullName[['num', 'location_grouped']],missing_fullName[['num', 'location_grouped']]]
output = pd.concat(output)
output.columns = ['accession', 'location']
output.to_csv("finalInput_1215.csv", index=False)
