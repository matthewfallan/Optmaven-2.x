import csv
import os
import sys

from Bio.PDB import PDBParser
import numpy as np

import standards


KEY_TOTAL = "Total"
KEY_TYPE = "Type"
FIELDS = ["Experiment", "Real", "User", "System", "Drive Usage", "Chains", "Residues", "Atoms", "Positions", "Selected Designs", "Min Relaxed Energy"]

# Get the name of the output file.
try:
    output = sys.argv[1]
    if os.path.exists(output):
        raise IOError("File exists: {}".format(output))
except IndexError:
    raise IOError("You must specify an output file.")


# Get benchmarking results from each experiment.
results = list()
for experiment in os.listdir(standards.ExperimentsDirectory):
    print(experiment)
    edir = os.path.join(standards.ExperimentsDirectory, experiment)
    bmark_file = os.path.join(edir, "Benchmarking.csv")
    summary_file = os.path.join(edir, "Summary.txt")
    results_file = os.path.join(edir, "Results.csv")
    # Skip experiments without benchmarking information.
    if not all([os.path.isfile(f) for f in [bmark_file, summary_file, results_file]]):
        continue
    with open(bmark_file) as f:
        for row in csv.DictReader(f):
            if row[KEY_TYPE] == KEY_TOTAL:
                results.append({k: v for k, v in row.items() if k in FIELDS})
                break
        else:
            raise ValueError("Cannot find {} row in file: {}".format(KEY_TOTAL, bmark_file))
    with open(summary_file) as f:
        positions, selected = None, None
        for line in f:
            if line.startswith("Total positions:"):
                positions = int(line.split(":")[1])
            elif line.startswith("Selected designs:"):
                selected = int(line.split(":")[1])
        if any([x is None for x in (positions, selected)]):
            raise ValueError("Cannot find all information in file: {}".format(summary_file))
    with open(results_file) as f:
        try:
            min_energy = min([float(row["Relaxed energy (kcal/mol)"]) for row in csv.DictReader(f)])
        except ValueError:
            min_energy = np.nan
    # Count the number of atoms and residues in the antigen.
    ag_file = os.path.join(standards.ExperimentsDirectory, experiment, "structures", "antigen_relaxed.pdb")
    structure = PDBParser().get_structure("antigen", ag_file)
    n_chains, n_atoms, n_residues = 0, 0, 0
    for chain in structure.get_chains():
        n_chains += 1
        n_residues += len(list(chain.get_residues()))
        n_atoms += len(list(chain.get_atoms()))
    results[-1].update({"Experiment": experiment, "Chains": n_chains, "Residues": n_residues, "Atoms": n_atoms, "Positions": positions, "Selected Designs": selected, "Min Relaxed Energy": min_energy})

# Output results as a CSV file.
with open(output, "w") as f:
    writer = csv.DictWriter(f, FIELDS)
    writer.writeheader()
    for experiment in sorted(results):
        writer.writerow(experiment)

