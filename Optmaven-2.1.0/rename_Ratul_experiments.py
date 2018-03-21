""" This script renames all of the experiments that are run on ACI-B (instead of Lion-XF) with an "ACI" extension. """

import os

with open("scripts/master2.sh") as f:
	scripts =[line.split("<")[1].strip() for line in f.readlines()]

for script in scripts:
	print(script)
	with open(script) as f:
		lines = f.readlines()
	if not os.path.isfile(script + ".back"):
		with open(script + ".back", "w") as f:
			f.write("".join(lines))
	if "ACI" not in lines[0]:
		lines[0] = lines[0].strip() + "_ACI\n"
	with open(script, "w") as f:
		f.write("".join(lines))
