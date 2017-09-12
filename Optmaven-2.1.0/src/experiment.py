""" This module defines the Experiment class of OptMAVEn. """

import cPickle as pkl
import os
import sys
import tempfile
import traceback

from Bio.PDB import Residue, Selection
from Bio.PDB.PDBIO import Select
import numpy as np

from console import disp
import klaus
import maps
import molecules
import standards
import submitter
import user_input


class Experiment(object):
    def __init__(self):
        do = True
        while do:
            name = self.ask("name")
            directory = os.path.join(standards.ExperimentsDirectory, name)
            if os.path.isdir(directory):
                disp("There is already an Experiment named {}.".format(name))
            else:
                do = False
        self.name = name
        self.name_contig = "_".join(name.split())
        self.directory = directory
        os.mkdir(self.directory)
        self.report_directory()
        self.temp = os.path.join(self.directory, ".temp")
        self.get_temp()
        self.file = os.path.join(self.temp, "{}.pickle".format(self.name_contig))
        self.warnings = os.path.join(self.directory, "warnings.txt")
        self.errors = os.path.join(self.directory, "errors.txt")
        #self.walltime = self.ask("walltime (seconds)", number=True)
        self.structure_directory = os.path.join(self.directory, "structures")
        os.mkdir(self.structure_directory)

    def save(self):
        with open(self.file, "w") as f:
            pkl.dump(self, f)

    def submit(self, args=None):
        handle, file_name = tempfile.mkstemp(prefix="{}_".format(self.name_contig), suffix=".sh", dir=self.temp)
        os.close(handle)
        if args is None:
            args = [""]
        print "ARGS:",args
        command = "\n".join(["""python {} {} {}""".format(os.path.realpath(__file__), self.file, arg) for arg in args])
        submitter.submit(file_name, command, self.walltime)

    def save_and_submit(self):
        self.save()
        self.submit()

    def ask(self, attribute, number=False):
        try:
            name = self.name
        except AttributeError:
            name = "this Experiment"
        prompt = "Please enter the {} of {}: ".format(attribute, name)
        if number:
            return user_input.get_number(prompt)
        else:
            return user_input.get(prompt)
    
    def report_directory(self):
        try:
            disp("The results of {} will be located in the following directory:\n{}".format(self.name, self.directory))
        except AttributeError:
            disp("This Experiment has no directory.")

    def get_temp(self):
        """ Return the directory of temporary files. If it does not exist, create it. This is to avoid returning a nonexistent directory. """
        if not os.path.isdir(self.temp):
            os.mkdir(self.temp)
        return self.temp

    def document_error(self, message):
        try:
            with open(self.errors, "a") as f:
                f.write(str(message))
        except OSError as e:
            error_file = "Optmaven_Experiment_errors.txt"
            with open(error_file, "a") as f:
                f.write("{}\n{}".format(message, e.message))


class OptmavenExperiment(Experiment):
    def __init__(self):
        Experiment.__init__(self)
        # Ask whether to customize advanced settings.
        if user_input.get_yn("Would you like to customize Optmaven settings? "):
            raise NotImplementedError("You cannot customize settings at this time.")
            #FIXME
        else:
            self.grid_x = standards.DefaultOptmavenGrid_x
            self.grid_y = standards.DefaultOptmavenGrid_y
            self.grid_z = standards.DefaultOptmavenGrid_z 
            self.grid_zAngle = standards.DefaultOptmavenGrid_zAngle
            self.topology_files = [standards.DefaultTopologyFile]
            self.parameter_files = [standards.DefaultParameterFile]
            self.solvation_files = [standards.DefaultSolvationFile]
            self.charmm_energy_terms = standards.DefaultCharmmEnergyTerms
            self.charmm_iterations = standards.DefaultCharmmIterations
            self.clash_cutoff = standards.DefaultClashCutoff
            self.walltime = standards.DefaultWalltime
            self.batch_size = standards.DefaultBatchSize
        # Define the antigen and epitope.
        entire_input_model = self.get_entire_input_model()
        antigen_input_chains = self.select_antigen_input_chains(entire_input_model)
        #antigen_input_residues = [residue for chain in antigen_input_chains for residue in self.select_antigen_input_residues(chain)]
        self.antigen_chain_ids = [chain.get_id() for chain in antigen_input_chains]
        self.select_epitope_residues(antigen_input_chains)
        self.write_antigen_input_residues(entire_input_model)#, antigen_input_residues)
        self.status = 1
        self.report_directory()
        self.save_and_submit()

    def run(self, args):
        tasks = {
            1: self.relax_antigen,
            2: self.one_maps_energy
        }
        try:
            try:
                task = tasks[self.status]
            except KeyError:
                raise ValueError("Bad Experiment status: {}".format(self.status))
            else:
                task(args)
        except Exception as e:
            tb = traceback.format_exc()
            self.document_error(tb)
        else:
            self.status += 1
            self.save_and_submit()

    def relax_antigen(self, args):
        antigen_molecule = molecules.Molecule(self.antigen_input_name, self.antigen_input_file, self)
        self.antigen_relaxed_name = "relaxed_antigen"
        self.antigen_relaxed_file = os.path.join(self.structure_directory, "antigen_relaxed.pdb")
        antigen_molecule.relax(self.antigen_relaxed_file)
        self.minimize_epitope_z_coordinates(antigen_molecule)

    def minimize_epitope_z_coordinates(self, antigen_molecule):
        # First move the antigen to the origin.
        antigen_molecule.translate_to(in_place=True)
        all_atoms = antigen_molecule.get_atom_array()
        epitope_ids = {(c_id, str(r_id[1])) for c_id, r_ids in self.epitope_residue_ids.items() for r_id in r_ids}
        epitope_selector = np.array([tuple(map(str, _id)) in epitope_ids for _id in all_atoms[["C", "R"]]])
        epitope_atoms = all_atoms[epitope_selector]
        coord_dims = ["x", "y", "z"]
        coord_n = len(coord_dims)
        all_coords = all_atoms[coord_dims].view(dtype=np.float).reshape(-1, coord_n)
        epi_coords = epitope_atoms[coord_dims].view(dtype=np.float).reshape(-1, coord_n)
        all_center = np.mean(all_coords, axis=0)
        epi_center = np.mean(epi_coords, axis=0)
        if np.allclose(all_center, epi_center, atol=0.001):
            # If the antigen and epitope centers are practically concurrent, then just use an null rotation.
            rot_matrix = np.eye(coord_n)
        else:
            epi_vector = epi_center - all_center
            # Make a rotation matrix that rotates the epitope vector to the negative z axis.
            neg_z_axis = np.array([0., 0., -1.])
            rot_axis = np.cross(epi_vector, neg_z_axis)
            if np.allclose(rot_axis, 0, atol=0.001):
                rot_matrix = np.eye(coord_n)
            else:
                rot_matrix = standards.rotate_vi_to_vf(epi_vector, neg_z_axis)
        self.epitope_zmin_file = os.path.join(self.structure_directory, "antigen_epitope_down.pdb")
        antigen_molecule.rotate(rot_matrix, self.epitope_zmin_file, in_place=True)
        # Rotate the antigen so that its z rotation angle is 0.
        #FIXME
        # Position the antigen using a grid search.
        self.position_antigen()
        
    def position_antigen(self):
        self.positions_file = os.path.join(self.temp, "positions.dat")
        with klaus.PositionAntigen(self) as x:
            pass
        self.status += 1
        self.all_maps_energies(None)

    def get_maps_part_energy_directory(self, part):
        return os.path.join(self.maps_energies_directory, part)

    def get_maps_part_energy_file_temp(self, part):
        return os.path.join(self.get_maps_part_energy_directory(part), "{}_energies_temp.dat".format(part))
    
    def get_maps_part_energy_file_finished(self, part):
        return os.path.join(self.maps_energies_directory, "{}_energies.dat".format(part))

    def all_maps_energies(self, args):
        """ Calculate the interacton energy between the antigen and all MAPs parts. """
        self.maps_energies_directory = os.path.join(self.temp, "maps_energies")
        if not os.path.isdir(self.maps_energies_directory):
            os.mkdir(self.maps_energies_directory)
        # Check which parts have not finished.
        unfinished_parts = [part for part in maps.parts if not os.path.isfile(self.get_maps_part_energy_file_finished(part))]
        unstarted_parts = [part for part in maps.parts if not os.path.isfile(self.get_maps_part_energy_file_temp(part))]
        #FIXME: better way to tell which parts are unstarted.
        collected_parts = list()
        for part in unstarted_parts:
            collected_parts.append(part)
            if len(collected_parts) == self.batch_size:
                self.submit(collected_parts)
                collected_parts = list()
        
    def one_maps_energy(self, args):
        """ Calculate the interacton energy between the antigen and one MAPs part. """
        try:
            part = args[2]
        except IndexError:
            #self.all_maps_energies(args)
        else:
            print "CALCULATING ENERGIES FOR", part
        
        
    def get_entire_input_model(self):
        # Select the antigen file.
        self.entire_input_name = "entire_input_file"
        do = True
        while do:
            entire_input_file = user_input.get_file("Please name the file containing the antigen: ", standards.PDBDirectory, fetch_pdb=True)
            try:
                entire_input_model = molecules.Molecule(self.entire_input_name, entire_input_file, self, exclude_hetero_ask=True).get_model()
            except Exception as error:
                disp("There was an error with the PDB import:\n{}".format(error.message))
            else:
                do = False
        return entire_input_model

    def select_antigen_input_chains(self, entire_input_model):
        # Select the antigen chains.
        chains = Selection.unfold_entities(entire_input_model, "C")
        chain_ids = [chain.get_id() for chain in chains]
        return user_input.select_from_list("Please select the antigen chains: ", chains, 1, None, names=chain_ids)
    
    def select_epitope_residues(self, chains):
        self.epitope_residue_ids = dict()
        for chain in chains:
            self.epitope_residue_ids[chain.get_id()] = user_input.select_from_list("Please select the epitope residues from chain {}: ".format(chain.get_id()), [residue.get_id() for residue in chain], 1, None, map(molecules.residue_code, chain))
        #return antigen_input_residues
    
    def write_antigen_input_residues(self, entire_input_structure):#, antigen_input_residues):
        self.antigen_input_name = "antigen_input_file"
        self.antigen_input_file = os.path.join(self.structure_directory, os.path.basename("antigen_input.pdb"))
        self.antigen_input_molecule = molecules.Molecule(self.antigen_input_name, self.antigen_input_file, self)
        selector = molecules.SelectChains(self.antigen_chain_ids)
        self.antigen_input_molecule.write_pdb(entire_input_structure, selector=selector)
    
    def get_molecules_list(self):
        status = self.status
        if status == 1:
            molecules_list = [self.antigen_input_molecule]
        else:
            raise ValueError("Invalid status: {}".format(status))
        missing = [molecule.file for molecule in molecules_list if not molecule.file_exists()]
        if len(missing) > 0:
            raise OSError("Missing molecule files: {}".format(", ".join(missing)))
        return molecules_list

"""
class ExperimentResidue(object):
    def __init__(self, biopython_residue):
        if not isinstance(biopython_residue, Residue.Residue):
            raise TypeError("Expected Bio.PDB.Residue.Residue, got {}".format(type(biopython_residue)))
        self.hetero, self.number, self.insertion = biopython_residue.get_id()
        self.type = biopython_residue.get_resname()
        self.chain_id = biopython_residue.get_parent().get_id()
        self.model_id = biopython_residue.get_parent().get_parent().get_id()

    def get_id(self):
        return (self.model_id, self.chain_id, self.number, self.insertion, self.type, self.hetero)

    def __str__(self):
        return "{}{}".format(self.number, self.insertion.strip())
    
    def __repr__(self):
        return str(self)


class SelectedResidueSelect(Select):
    def __init__(self, selected_residues):
        Select.__init__(self)
        self.selected_residues = {residue.get_id() for residue in selected_residues}

    def accept_residue(self, residue):
        return int(ExperimentResidue(residue).get_id() in self.selected_residues)
"""

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise OSError("Usage: python experiments.py /path/to/experiment.pickle")
    experiment_file = sys.argv[1]
    with open(experiment_file) as f:
        experiment = pkl.load(f)
    experiment.run(sys.argv)
