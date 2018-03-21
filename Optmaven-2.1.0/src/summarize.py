""" Create a table of summary data for experiments. """

import csv
import os
import sys

import standards

output = sys.argv[1]
if os.path.exists(output):
    raise IOError("File exists: {}".format(output))

fields = ["Experiment", "Chains", "Epitope Residues"]
keys = {"Experiment name": "Experiment",
        "Antigen chains": "Chains",
        "Epitope residues": "Epitope Residues"}

experiments = list()
for experiment in os.listdir(standards.ExperimentsDirectory):
    exp_dir = os.path.join(standards.ExperimentsDirectory, experiment)
    summary_file = os.path.join(exp_dir, "Summary.txt")
    if not os.path.isfile(summary_file):
        continue
    with open(summary_file) as f:
        values = dict()
        for line in f:
            if ":" not in line:
                 continue
            key, value = line.split(":", 1)[0].strip(), "".join(line.split(":", 1)[1:]).strip()
            field = keys.get(key, None)
            if field == "Epitope Residues":
                value = ", ".join([x.strip() for x in value.split(",")])
            if field is not None:
                values[field] = value
        experiments.append(values)

with open(output, "w") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for experiment in sorted(experiments, key=lambda x: x["Experiment"]):
        writer.writerow(experiment)
