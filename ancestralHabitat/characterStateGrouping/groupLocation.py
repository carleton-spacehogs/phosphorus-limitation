# Author: Zhanghan Tony Ni, December 2022.
# Assigns specific geographic locations to general character states.
# Work for Professor Rika Anderson for Phosphorus limitation porject.

from webbrowser import get
import pandas as pd
import ast

#raw location data
df = pd.read_csv("locations_1214.csv")

#lists for categorization
marine_deep = ["Ridge", "sea", "hydrothermal vent", "coral", "seawater", "marine", "ocean", "deferrisoma camini", "crenarchaeota archaeon", "sponge"]
marine_shallow = ["surface","photic zone","Reef","Coal Oil Point","Plankton", "planktonic","bay"]
marine = marine_deep+marine_shallow
terrestrial = ["Salinibacter ruber","almond", "Algeria","Forest", "deep gas storage aquifer", "peat","Kamchatka","Biocathode", "hot spring","hot springs", "lake", "lakes", "well", "stream", "river","aquatic metagenome","groundwater", "freshwater","Siberia", "Gori-gun", "oilfield", "tailings", "Leusden", "Cave", "wastewater","air-conditioning", "mud","heroin","rock", "terrestrial", "coal", "mine", "beach", "facility", "zoo","microbial mat metagenome", "compost", "hydrocarbon", "sludge", "halite", "plant", "effluent", "soil", "fermentation metagenome", "hydrocarbon", "sediment", "bioreactor metagenome", "ammonia-oxidizing enrichment culture", "subsurface", "digester", "ovicells" ]
hostAssociated = ["mouth","oral","intestinal","vagina","blood","tissue", "Tick", "teat", "salmon", "milk", "choanae", "gastrointestinal", "choanae","cattle", "cavia porcellus","mucus", "rumen", "female", "human", "gut", "faeces"]

# additional lists for refined categorization
# marine = ["sea", "hydrothermal vent", "coral", "seawater", "marine", "bay", "ocean", "deferrisoma camini", "crenarchaeota archaeon", "sponge"]
# land = ["well", "lake", "lakes", "rock", "terrestrial", "coal","hot spring","mine","beach", "facility", "zoo", "stream", "microbial mat metagenome", "river", "aquatic metagenome", "compost", "hydrocarbon", "sludge", "halite", "plant", "effluent", "soil", "groundwater", "fermentation metagenome", "hydrocarbon", "sediment", "hot springs", "bioreactor metagenome", "ammonia-oxidizing enrichment culture", "subsurface", "digester", "freshwater", "ovicells"]
# hostAssociated = ["cattle", "cavia porcellus","mucus", "rumen", "female", "human", "gut", "faeces", "symbiont", "ecological metagenome", "biofilm"]
# fresh = ["lake", "lakes", "well", "stream", "river","aquatic metagenome","groundwater", "freshwater", ]
# saltLakes = []
# hotSpring = ["hot spring","hot springs"]

#a special case requiring additional processing
tara = pd.read_csv('Tara_bins_origin.csv')
tara = tara.replace('-', '', regex=True).astype(str)
def loc(depth):
    if depth == "MES":
        return "marine_deep"
    if depth == "DCM" or depth == "SRF":
        return "marine_shallow"
    return "ungrouped"
dict = {}
for i in range(0, len(tara)):
    dict[tara.loc[i,'bin']] = loc(tara.loc[i, 'depth'])

#assign character states
def getGroup(s):
    if s == "[]":
        return "ungrouped"
    s = s.lower()
    listRep = ast.literal_eval(s)
    if 'marine pelagic biome' in s:
        bin = listRep[0].upper()
        if bin in dict.keys():
            return dict[bin]
    for h in hostAssociated:
        if h.lower() in s:
            return "hostAssociated"
    for ms in marine_shallow:
        if ms.lower() in s:
            return "marine_shallow"
    for md in marine_deep:
        if md.lower() in s:

            return "marine_deep"
    for t in terrestrial:
        if t.lower() in s:
            return "terrestrial"
    
    #additional conditions for refined categorization
    # for f in fresh:
    #     if f.lower() in s:
    #         return "fresh"
    # for s in saltLakes:
    #     if s.lower() in s:
    #         return "saltLakes"
    # for h in hotSpring:
    #     if h.lower() in s:
    #         return "hotSpring"
    return "ungrouped"

#write to file
for i in range(0, len(df)):
    location = getGroup((df.loc[i, 'location']))
    df.loc[i, 'location_grouped'] = location
df.to_csv("locations_grouped_1214_v2.csv", index=False)




