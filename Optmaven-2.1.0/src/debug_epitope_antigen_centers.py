import numpy as np

import sys

pdb = sys.argv[1]
epi_file = sys.argv[2]

with open(epi_file) as f:
	epitope = [tuple(line.split()) for line in f]

all_coords = list()
epi_coords = list()
with open(pdb) as f:
	for line in f:
		if not line.startswith(("ATOM", "HETATM")):
			continue
		x = float(line[30: 38])
		y = float(line[38: 46])
		z = float(line[46: 54])
		coord = [x, y, z]
		all_coords.append(coord)
		chain = line[21]
		resnum = line[22: 26].strip()
		if (chain, resnum) in epitope:
			epi_coords.append(coord)
all_coords = np.array(all_coords)
epi_coords = np.array(epi_coords)
all_center = np.mean(all_coords, axis=0)
epi_center = np.mean(epi_coords, axis=0)
print "CENTER", all_center
print "EPITOPE", epi_center
print "DIFF", epi_center - all_center
