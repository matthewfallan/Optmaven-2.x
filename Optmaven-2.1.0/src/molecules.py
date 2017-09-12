""" This modules manages molecules. """

from collections import defaultdict
import itertools
import os
import warnings

from Bio.PDB import NeighborSearch, MMCIFParser, PDBParser, Residue, Selection, Structure
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
import numpy as np

import charmm
from console import disp
import user_input


_MAX_RANDOM_RESIDUE = 2**30 - 1


class Molecule(object):
    def __init__(self, name, file_, experiment, exclude_hetero_ask=False):
        if name != "".join(name.split()) or os.path.sep in name:
            raise ValueError("Molecule names may not contain whitespace or {}.".format(os.path.sep))
        self.name = name
        if not isinstance(file_, str):
            raise TypeError("Expected str, got {}.".format(type(file_)))
        self.file = file_
        self.exclude_hetero_ask = exclude_hetero_ask
        self.excluded_residue_ids = None
        self.experiment = experiment
        self.dims = ["x", "y", "z"]
        self.ndim = len(self.dims)
    
    def get_structure(self, file_name=None):
        if file_name is None:
            file_name = self.file
        structure = get_parser(file_name).get_structure(self.name, file_name)
        models = Selection.unfold_entities(structure, "M")
        if len(models) != 1:
            raise NotImplementedError("Optmaven cannot currently handle structures with {} models.".format(len(models)))
        if self.exclude_hetero_ask and self.excluded_residue_ids is None:
            self.excluded_residue_ids = dict()
            for chain in Selection.unfold_entities(models[0], "C"):
                hetero_residues = [residue for residue in chain if residue.get_id()[0] != " "]
                hetero_residue_ids = [residue.get_id() for residue in hetero_residues]
                if len(hetero_residue_ids) > 0:
                    disp("Excluding hetero-residues from chain {}.".format(chain.get_id()))
                    self.excluded_residue_ids[chain.get_id()] = set(user_input.select_from_list("Please select hetero-residues to exclude from {}: ".format(self.name), hetero_residue_ids, names=map(residue_code, hetero_residues)))
        if self.excluded_residue_ids is not None:
            for chain in Selection.unfold_entities(models[0], "C"):
                exclude = self.excluded_residue_ids.get(chain.get_id())
                if exclude is not None:
                    [chain.detach_child(res_id) for res_id in exclude]
        return structure

    def get_model(self):
        return Selection.unfold_entities(self.get_structure(), "M")[0]
 
    def get_file_base(self):
        return os.path.basename(self.file)

    def get_file_split(self):
        return os.path.splitext(self.get_file_base())

    def get_file_name(self):
        return self.get_file_split()[0]

    def get_chains(self):
        return Selection.unfold_entities(self.get_model(), "C")

    def file_exists(self):
        return os.path.isfile(self.file)

    def get_renumbered_structure(self, first_residue=1):
        structure = self.get_structure()
        renumber_structure(structure, first_residue)
        return structure

    def write_pdb(self, structure, file_name=None, selector=None, preserve_atom_numbering=False):
        if file_name is None:
            file_name = self.file
        write_pdb(structure, file_name, selector, preserve_atom_numbering=preserve_atom_numbering)
   
    def write_renumbered_pdb(self, file_name, first_residue=1, selector=None):
        self.write_pdb(self.get_renumbered_structure(first_residue), file_name, selector)
        
    def relax(self, relaxed_file=None, in_place=True):
        """ Perform a structural relaxation. """
        self.relax_CHARMM(relaxed_file, in_place)
    
    def relax_CHARMM(self, relaxed_file=None, in_place=True):
        if relaxed_file is None:
            if not in_place:
                raise ValueError("A relaxation must be in place if no alternate relaxed file is given.") 
            relaxed_file = self.file
        with charmm.Relaxation(self.experiment, [self], relaxed_file) as relax:
            pass
        if in_place:
            if relaxed_file != self.file:
                self.file = relaxed_file

    def disassemble_for_CHARMM(self, structure=None, directory=None):
        # Rename all of the residues so that they correspond to the CHARMM naming conventions.
        if structure is None:
            structure = self.get_structure()
        if directory is None:
            directory = self.directory
        renumber_structure(structure)
        his_to_hsd(structure)
        models = structure.get_list()
        if len(models) > 1:
            raise NotImplementedError("Optmaven cannot currently handle Molecules with {} models.".format(len(models)))
        chains = (models[0]).get_list()
        files = list()
        for chain in chains:
            _id = chain.get_id()
            selector = SelectChains([_id])
            base_name = charmm.make_pdb_name(_id, "in")
            file_name = os.path.join(directory, base_name)
            self.write_pdb(structure, file_name, selector, preserve_atom_numbering=True)
            files.append(file_name)
        return chains, files

    def assemble_from_CHARMM(self, chain_ids, files, residue_numbers_by_model):
        # After running CHARMM, assemble the output files into the molecule.
        # Thee unique part (rather than a simple merge) is that the molecules lose their chain ids within charmm, so these ids must be added back. 
        assembled_structure = None
        assembled_model = None
        for _id, _file in zip(chain_ids, files):
            structure = self.get_structure(_file)
            models = structure.get_list()
            if len(models) > 1:
                raise ValueError("A file output from CHARMM should have 1 model, not {}.".format(len(models)))
            model = models[0]
            residue_numbers_by_chain = residue_numbers_by_model[model.get_id()]
            chains = model.get_list()
            if len(chains) > 1:
                raise ValueError("A file output from CHARMM should have 1 chain, not {}.".format(len(chains)))
            chain = chains[0]
            chain.id = _id
            residue_numbers = residue_numbers_by_chain[_id]
            renumber_chain(chain, residue_numbers=residue_numbers)
            if assembled_structure is None:
                assembled_structure = structure
                assembled_model = assembled_structure.get_list()[0]
            else:
                chain.detach_parent()
                assembled_model.add(chain)
        self.write_pdb(assembled_structure)
        return assembled_structure

    def get_coords(self):
        return self.get_atom_array()[self.dims].view(dtype=np.float).reshape(-1, self.ndim)

    def translate_to(self, point=None, file_name=None, in_place=None):
        if point is None:
             point = np.array([0., 0., 0.])
        point = np.array(point)
        coords = self.get_coords()
        center = np.mean(coords, axis=0)
        self.translate(point - center, file_name, in_place)

    def translate(self, vector, file_name=None, in_place=False):
        if file_name is None:
            if in_place is False:
                raise ValueError("A translation file must be given if the translation is to not be in place.")
            file_name = self.file
        coords = self.get_coords()
        vector = np.array(vector).reshape((1, self.ndim))
        trans_coords = coords + vector
        atoms = self.get_atoms()
        # Replace all of the atomic coordinates.
        for coord, atom in zip(trans_coords, atoms):
            atom.set_coord(coord)
        # Write the rotated coordinates.
        structure = atom.get_parent().get_parent().get_parent().get_parent()
        self.write_pdb(structure, file_name)
        if in_place:
            self.file = file_name

    def rotate(self, rot_matrix, file_name=None, in_place=False):
        if file_name is None:
            if in_place is False:
                raise ValueError("A rotation file must be given if the rotation is to not be in place.")
            file_name = self.file
        rot_coords = np.dot(self.get_coords(), rot_matrix)
        atoms = self.get_atoms()
        # Replace all of the atomic coordinates.
        for coord, atom in zip(rot_coords, atoms):
            atom.set_coord(coord)
        # Write the rotated coordinates.
        structure = atom.get_parent().get_parent().get_parent().get_parent()
        self.write_pdb(structure, file_name)
        if in_place:
            self.file = file_name

    def get_atoms(self):
        structure = self.get_structure()
        chains = (structure.get_list()[0]).get_list()
        atoms = [a for c in chains for r in c for a in r]
        return atoms
        
    def get_atom_array(self):
        raw_atoms = [(a.get_parent().get_parent().get_id(), str(a.get_parent().get_id()[1]), a.get_id(), a.coord[0], a.coord[1], a.coord[2]) for a in self.get_atoms()]
        c_format = "S{}".format(max([len(a[0]) for a in raw_atoms]))
        r_format = "S{}".format(max([len(a[1]) for a in raw_atoms]))
        a_format = "S{}".format(max([len(a[2]) for a in raw_atoms]))
        atom_format = [("C", c_format), ("R", r_format), ("A", a_format), (self.dims[0], float), (self.dims[1], float), (self.dims[2], float)]
        atom_array = np.array(raw_atoms, dtype=atom_format)
        return atom_array

    def interchain_residue_contacts(self, chain_ids_1, chain_ids_2, radius):
        """ Generate a list of residue contacts between two chains. """
        all_chains = {chain.get_id(): chain for chain in self.get_chains()}
        selected_chains = {chain_id: chain for chain_id, chain in all_chains.items() if chain_id in chain_ids_1 + chain_ids_2}
        atoms = [atom for chain_id, chain in selected_chains.items() for atom in Selection.unfold_entities(chain, "A")]
        residue_contacts = NeighborSearch(atoms).search_all(radius, "R")
        classified_contacts = defaultdict(list)
        for contact in residue_contacts:
            chain_1, chain_2 = [residue.get_parent().get_id() for residue in contact]
            if chain_1 in chain_ids_1 and chain_2 in chain_ids_2:
                classified_contacts[(chain_1, chain_2)].append({chain_1: contact[0], chain_2: contact[1]})
            elif chain_2 in chain_ids_1 and chain_1 in chain_ids_2:
                classified_contacts[(chain_2, chain_1)].append({chain_1: contact[0], chain_2: contact[1]})
        return classified_contacts


class SelectChains(Select):
    def __init__(self, chain_ids):
        Select.__init__(self)
        if not isinstance(chain_ids, (list, tuple, set)):
            raise TypeError("Chain IDs must be a list, tuple, or set, not {}.".format(type(chain_ids)))
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain):
        return int(chain.get_id() in self.chain_ids)


def write_pdb(structure, file_name, selector=None, preserve_atom_numbering=False):
    writer = PDBIO()
    writer.set_structure(structure)
    if selector is None:
        writer.save(file_name, preserve_atom_numbering=preserve_atom_numbering)
    else:
        writer.save(file_name, selector, preserve_atom_numbering=preserve_atom_numbering)


def his_to_hsd(structure):
    rename_residues(structure, {"HIS": "HSD"})


def hsd_to_his(structure):
    rename_residues(structure, {"HSD": "HIS"})


def rename_residues(structure, names):
    # names is a dict of {old_name: new_name}
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.resname
                try:
                    residue.resname = names[resname]
                except KeyError:
                    pass

"""
def rename_residues(_file, residues):
    # Rename all of the residues of a given name(s) to other names.
    with open(_file) as f:
        contents = f.read()
    for old_name, new_name in residues.items():
        contents = contents.replace(old_name, new_name)
    with open(_file, "w") as f:
        f.write(contents)
"""

def renumber_structure(structure, first_residue=1):
    resnum = first_residue
    for model in structure:
        resnum = renumber_model(model, resnum) + 1


def renumber_model(model, first_residue=1):
    resnum = first_residue
    for chain in model:
        resnum = renumber_chain(chain, resnum) + 1
    return resnum


def renumber_chain(chain, first_residue=1, residue_numbers=None):
    residues = chain.get_list()
    resids = [residue.get_id() for residue in residues]
    resids_set = set(resids)
    if len(resids) != len(resids_set):
        raise ValueError("Repeat residue IDs detected.")
    if residue_numbers is not None and len(residue_numbers) != len(resids):
        raise ValueError("Residue numbers has length {} but residues has length {}.".format(len(residue_numbers), len(resids)))
    for index, residue in enumerate(chain):
        if residue_numbers is not None:
            resnum_new = residue_numbers[index]
        else:
            resnum_new = index + first_residue
        resid_new = get_renumbered_resid(residue, resnum_new)
        if resid_new == resids[index]:
            # Skip this residue if its new number matches the old.
            continue
        while resid_new in resids_set:
            # Biopython will raise an exception if a residue is renumbered such that its number matches the number of any existing residue.
            # If it matches, then temporarily renumber the matching residue to something random.
            matching_residue_index = resids.index(resid_new)
            if matching_residue_index <= index:
                raise ValueError("Residue numbers are out of order.")
            resid_temp = get_renumbered_resid(residue, np.random.randint(0, _MAX_RANDOM_RESIDUE))
            while resid_temp in resids_set:
                resid_temp = get_renumbered_resid(residue, np.random.randint(0, _MAX_RANDOM_RESIDUE))
            renumber_residue(residues[matching_residue_index], resid_temp[1])
            resids_set.remove(resid_new)
            resids_set.add(resid_temp)
            resids[matching_residue_index] = resid_temp
        # Renumber the residue.
        renumber_residue(residue, resnum_new)
        resids_set.remove(resids[index])
        resids_set.add(resid_new)
        resids[index] = resid_new
    return resnum_new


def get_renumbered_resid(residue, number):
    return (residue.get_id()[0], number, residue.get_id()[2])


def renumber_residue(residue, number, increment=False):
    # FIXME: check type of structure.
    if increment:
        residue.id = get_renumbered_resid(residue, residue.get_id()[1] + number)
    else:
        residue.id = get_renumbered_resid(residue, number)


def merge(molecule_list, return_molecule=False, merged_name=None, merged_file=None, write_pdb=False):
    """ Merge all of the molecules in a list into one Molecule. """
    # First, ensure that all of the items in molecules are actually Molecules.
    if not all([isinstance(molecule, Molecule) for molecule in molecule_list]):
        raise TypeError("Expected a list of Molecule objects.")
    # Merge all of the chains into one model in one structure.
    merged_structure = None
    exp = None
    for molecule in molecule_list:
        if exp is None:
            exp = molecule.experiment
        elif molecule.experiment != exp:
            raise ValueError("Cannot merge two Molecules from different Experiments.")
        molecule_structure = molecule.get_structure()
        molecule_models = list(molecule_structure)
        if merged_structure is None:
            merged_structure = molecule_structure
            merged_model = molecule_models[0]
        else:
            for chain in molecule_models[0]:
                chain_id = chain.get_id()
                if chain_id in [merged_chain.get_id() for merged_chain in merged_model]:
                    raise ValueError("Attempted to merge two Molecules that each have a chain {}.".format(chain_id))
                chain.detach_parent()
                merged_model.add(chain)
    if return_molecule:
        if not all([isinstance(arg, str) for arg in [merged_name, merged_file]]):
            raise ValueError("A name and file must be given to make a Molecule.")
        merged_molecule = Molecule(merged_name, merged_file, exp)
        if write_pdb:
            merged_molecule.write_pdb(merged_structure)
        return merged_molecule
    else:
        return merged_structure


def residue_code(residue):
    if isinstance(residue, Residue.Residue):
        return "{}{}".format(residue.get_id()[1], residue.get_id()[2]).strip()
    else:
        raise TypeError("Expected a Residue, got {}.".format(type(residue)))


def get_parser(file_):
    """ Get a parser appropriate to the file format. """
    try:
        file_base, ext = os.path.splitext(file_)
    except ValueError:
        raise ValueError("Cannot obtain extension of file {}".format(file_))
    else:
        try:
            return {
                ".pdb": PDBParser(),
                ".cif": MMCIFParser()
            }[ext]
        except KeyError:
            raise ValueError("Unknown molecular file format: {}".format(ext))
 

if __name__ == "__main__":
    m1 = Molecule("m1", "../data/pdbs/1ce1.pdb")
    m2 = Molecule("m2", "../data/pdbs/1tzi.pdb")
    m1.relax()
