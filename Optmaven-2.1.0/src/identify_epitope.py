""" This script identifies the epitope residues of an antigen in complex with an antibody. """

import os
import re

import EXPERIMENT
import MOLECULES
import STANDARDS
import USER_INPUT

def residue_compare(residue_code_1, residue_code_2):
    pattern = re.compile("([0-9]+)([A-Z]?)")
    i1, i2 = [(int(groups[0]), groups[1]) for groups in [pattern.match(code).groups() for code in [residue_code_1, residue_code_2]]]
    return cmp(i1, i2)
    

# Make a Molecule for the antigen-antibody complex.
complex_name = "antigen_antibody_complex"
complex_file = USER_INPUT.get_file("Please type the name of the antigen-antibody complex: ", STANDARDS.PDBDirectory)
complex_molecule = MOLECULES.Molecule(complex_name, complex_file)
complex_chains = complex_molecule.get_chains()

# Select the chains.
chains = [chain.get_id() for chain in complex_chains]
antigen_chains = USER_INPUT.select_from_list("Please select the antigen chain(s): ", chains, 1, None)
lAg = len(antigen_chains)
antibody_chains = USER_INPUT.select_from_list("Please select the antibody chain(s): ", [chain_id for chain_id in chains if chain_id not in antigen_chains], 1, None)
lAb = len(antibody_chains)
"""
for index in range(lAg + lAb):
    if index < lAg:
        chain = antigen_chains[index]
    else:
        chain = antibody_chains[index]
    disp("Selecting residues from chain {}.".format(chain.get_id()))
    residues = map(EXPERIMENT.ExperimentResidue, chain)
    # Hetero-residues (i.e. those that are not amino acids) may be excluded.
    hetero_residues = [residue for residue in residues if residue.hetero != " "]
    excluded_hetero_residues = set(USER_INPUT.select_from_list("Please select hetero-residues to exclude from the antigen: ", hetero_residues))
    antigen_input_residues = [residue for residue in residues if residue not in excluded_hetero_residues]
"""
radius = USER_INPUT.get_number("Please enter the contact cutoff radius: ", 0.0, None)
#epitope_file = USER_INPUT.get_file("Please enter the name of the epitope output file: ", os.getcwd(), new_file=True)
epitope_file = os.path.join(os.getcwd(), complex_molecule.get_file_name() + "_epitope.txt")
contacts = complex_molecule.interchain_residue_contacts(antigen_chains, antibody_chains, radius)
epitope_residues = set()
with open(epitope_file, "w") as f:
    text = """PDB\t{}
Antigen\t{}
Antibody\t{}
Cutoff\t{}""".format(complex_molecule.get_file_name().upper(), ", ".join(antigen_chains), ", ".join(antibody_chains), radius)
    for chain in antigen_chains:
        chain_pairs = [pair for pair in contacts.keys() if chain in pair]
        text += "\n" + chain + "\t" + ", ".join(sorted({MOLECULES.residue_code(contact[chain]) for chain_pair in chain_pairs for contact in contacts[chain_pair]}, cmp=residue_compare))
    f.write(text)
