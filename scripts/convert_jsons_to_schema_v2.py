#!/usr/bin/env python

import sys
import json
import os
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '..')
from correctionlib.schemav2 import Binning, Category, Correction, CorrectionSet

# TODO: how do we want to report uncertainties? value and systs separately, or value +/- syst already?
# def build_uncertainties(sf):
#     keys = ["nominal"]
#     keys += ["up_syst", "down_syst"] if "syst" in sf else []
#     keys += ["up_stat", "down_stat"] if "stat" in sf else []
    
#     content = [sf["value"]]
#     content += [sf["value"] + sf["syst"], sf["value"] - sf["syst"]] if "syst" in sf else []
#     content += [sf["value"] + sf["stat"], sf["value"] - sf["stat"]] if "stat" in sf else []
    
#     return Category.parse_obj({
#         "nodetype": "category",
#         "keys": keys,
#         "content": content
#     })

def build_content(sf, binning_array):

    bin_strings = {}
    binning = {}
    for element in binning_array:
        bin_edges = element['binning']
        bin_var = element['variable']
        binning[bin_var] = bin_edges
        bin_strings[bin_var] = []
        for idx in range(len(bin_edges)-1):
            bin_strings[bin_var] += [bin_var + ':[' + str(bin_edges[idx]) + ',' + str(bin_edges[idx+1]) + ']']

    dimensions = len(binning)
    bin_vars = list(binning.keys())

    # iterate through the bin indices
    # this does nested for loops of the N-D binning (e.g. pt, eta)
    # binning starts at 1 (0 is underflow), same as ROOT
    indices = [list(range(1, len(binning[bin_var]))) for bin_var in bin_vars]

    def build_schema_recursively(dim, index):

        # If we reach recursion bottom, build and return the systematics node
        if dim == dimensions + 1:

            index_to_string_mapping = {}
            for d, bin_var in enumerate(bin_vars):
                index_to_string_mapping[bin_var] = bin_strings[bin_var][index[d]-1]

            def get_systematics(sf, keys):
                for key in keys.values():
                    sf = sf[key]
                return sf

            systematics = get_systematics(sf, index_to_string_mapping)

            keys, content = [], []
            for syst, value in systematics.items():
                keys.append(syst)
                syst = syst if syst != "value" else "nominal"
                content.append({"key": syst, "value": value}) 
            return Category.parse_obj({
                "nodetype": "category",
                "input": "scale_factors",
                "content": content
            })

        # If not, build a binning node
        edges = list(map(float, binning[bin_vars[dim-1]]))
        content = [build_schema_recursively(dim+1, tuple(list(index)[0:dim-1]+[i]+list(index)[dim:])) for i in indices[dim-1]]
        return Binning.parse_obj({
            "nodetype": "binning",
            "input": bin_vars[dim-1],
            "edges": edges,
            "flow": "error",
            "content": content,
        })

    content = build_schema_recursively(1, tuple([1] * dimensions))
    return content

if __name__ != "__main__" or len(sys.argv) < 2:
    print(f'Please run this script as {sys.argv[0]} dir_with_json_dirs (output of spark_tnp)')
    sys.exit(1)
else:
    rootdir = sys.argv[1]

all_json_files = []
for root, subdirs, files in os.walk(rootdir):
    for subdir in subdirs:
        for subroot, subsubdirs, subfiles in os.walk(os.path.join(rootdir, subdir)):
            json_files = [os.path.join(subroot, subfile) for subfile in subfiles if subfile.endswith('.json') and 'schemaV1' not in subfile]
            all_json_files += json_files

all_corrections = []

for json_file in all_json_files:
    with open(json_file) as f:

        print('Processing ', json_file)
        
        sf = json.load(f)

        sf_name = list(sf.keys())[0]
        sf_description = sf_name
        sf_vars_string = list(sf[sf_name].keys())[0]

        binning_array = sf[sf_name][sf_vars_string].pop('binning')

        bin_vars = []
        for binning in binning_array:
            bin_vars.append(binning['variable'])

        inputs = [{"name": bin_var, "type": "real", "description": "Probe " + bin_var} for bin_var in bin_vars]
        inputs += [{"name": "scale_factors", "type": "string", "description": "Choose nominal scale factor or one of the uncertainties"}]
        
        data = build_content(sf[sf_name][sf_vars_string], binning_array)

        corr = Correction.parse_obj({
            "version": 1,
            "name": sf_name,
            "description": sf_description,
            "inputs": inputs,
            "output": {"name": "weight", "type": "real", "description": "Output scale factor (nominal) or uncertainty"},
            "data": data
        })

        all_corrections.append(corr)

cset = CorrectionSet.parse_obj({
    "schema_version": 2,
    "corrections": all_corrections
})

# Write out converted json
# with open(os.path.splitext(file_name)[0]+'_schemaV2.json', "w") as fout:
with open('output' + '_schemaV2.json', "w") as fout:
    fout.write(cset.json(exclude_unset=True, indent=4))