#!/usr/bin/env python

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '..')
from correctionlib.correctionlib.schemav1 import Binning, Category, Correction, CorrectionSet

def build_uncertainties(sf):
    keys = ["nominal"]
    keys += ["up_syst", "down_syst"] if "syst" in sf else []
    keys += ["up_stat", "down_stat"] if "stat" in sf else []
    
    content = [sf["value"]]
    content += [sf["value"] + sf["syst"], sf["value"] - sf["syst"]] if "syst" in sf else []
    content += [sf["value"] + sf["stat"], sf["value"] - sf["stat"]] if "stat" in sf else []
    
    return Category.parse_obj({
        "nodetype": "category",
        "keys": keys,
        "content": content
    })


def parse_str(key, prefix="abseta:"):
    if not key.startswith(prefix + "["):
        raise ValueError(f"{key} missing prefix {prefix}")
    lo, hi = map(float, key[len(prefix + "["):-1].split(","))
    return lo, hi


def build_pts(sf):
    # Could happen that higher pt bin edge comes
    # lexicographically before lower one, so sort bins first
    sf_sorted_data = {}
    sf_sorted_hi = {}
    for binstr, data in sf.items():
        if not binstr.startswith("pt:["):
            raise ValueError
        lo, hi = map(float, binstr[len("pt:["):-1].split(","))
        sf_sorted_data[lo] = data
        sf_sorted_hi[lo] = hi

    edges = []
    content = []
    for i in sorted(sf_sorted_data):
        lo = i
        data = sf_sorted_data[i]
        if len(edges) == 0:
            edges.append(lo)
        if edges[-1] != lo:
            raise ValueError
        edges.append(sf_sorted_hi[lo])
        content.append(build_uncertainties(data))
    
    return Binning.parse_obj({
        "nodetype": "binning",
        "edges": edges,
        "content": content,
    })


def build_etas(sf):
    bins = [parse_str(s, "abseta:") for s in sf]
    edges = sorted(set(edge for bin in bins for edge in bin))
    content = [None] * (len(edges) - 1)
    for s, data in sf.items():
        lo, hi = parse_str(s, "abseta:")
        found = False
        for i, bin in enumerate(bins):
            if bin[0] >= lo and bin[1] <= hi:
                content[i] = build_pts(data)
                found = True

        
    return Binning.parse_obj({
        "nodetype": "binning",
        "edges": edges,
        "content": content,
    })


import json, os

if __name__ != "__main__" or len(sys.argv) < 2:
    print(f'Please run this script as {sys.argv[0]} json_to_be_converted')
    sys.exit(1)
else:
    file_name = sys.argv[1]
    
# Load data
if 'https://' in file_name:
    import requests
    sf = requests.get(file_name).json()
else:
    with open(file_name) as f:
        sf = json.load(f)

        
# Run conversion to schema v1
sf_name = list(sf.keys())[0]
sf_type = 'abseta_pt'
sf_description = sf_name
# Remove existing binning description, it is built automatically later
sf[sf_name][sf_type].pop('binning')

corr = Correction.parse_obj(
    {
        "version": 1,
        "name": sf_name,
        "description": sf_description,
        "inputs": [
            {"name": "abseta", "type": "real"},
            {"name": "pt", "type": "real"},
            {"name": "uncertainties", "type": "string"},
        ],
        "output": {"name": "weight", "type": "real"},
        "data": build_etas(sf[sf_name][sf_type]),
    }
)

cset = CorrectionSet.parse_obj({
    "schema_version": 1,
    "corrections": [
        corr
    ]
})

# Write out converted json
with open(os.path.splitext(file_name)[0]+'_schemaV1.json', "w") as fout:
    fout.write(cset.json(exclude_unset=True, indent=4))