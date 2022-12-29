# Author: Zhanghan Tony Ni, December 2022.
# Obtains geographic location from metadata in json format.
# Work for Professor Rika Anderson for Phosphorus limitation porject.
import json
import os
import csv

def checkKey(dic, key):
    if key in dic.keys():
        return True

#extract location info from metadata files 
directory = '/Users/zhanghanni/Desktop/ToL Location/raw_json_v2'
dict = {}
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        file = open(f)
        data = json.load(file)
        value = data['characteristics']
        keys = ["isolate","organism","env_medium", "environmental medium", "isolation site", "Source of Isolate", "geo loc name" ,"isolation source", "isolation-source", "isolation-source", "metagenome source", "env_material", "project name", "env biome"]
        locations = []
        for key in keys:
            if checkKey(value, key):
                locations.append(value[key][0]['text'])
        name = ''
        for i in filename:
            name = name+i
            if i == '_':
                break

        dict[name[:-1]] = locations
        file.close()

with open('locations_1215.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    writer.writerow(["ncbi_biosample", "location"])
    for key, value in dict.items():
       writer.writerow([key, value])