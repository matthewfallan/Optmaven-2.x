""" This module defines the Experiment class of OptMAVEn. """

from collections import defaultdict, OrderedDict
import cPickle as pkl
#import multiprocessing
import os
import shutil
import sys
import tempfile
import time
import traceback

from Bio.PDB import Residue, Selection
from Bio.PDB.PDBIO import Select
from Bio import SeqIO
import numpy as np

import charmm
from console import clear, disp
import klaus
import kmeans
import maps
import molecules
import standards
import submitter
import user_input


class Experiment(object):
    def __init__(self):
        clear()
        do = True
        while do:
            name = self.ask("name", valid_path=True)
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

    def submit(self, args=None, jobs=None, options=None, queue=True):
        if args is None:
            args = [""]
        if jobs is None:
            handle, file_name = tempfile.mkstemp(prefix="{}_".format(self.name_contig), suffix=".sh", dir=self.temp)
            os.close(handle)
            command = "\n".join(["{} {} {} {}".format(standards.PythonCommand, os.path.realpath(__file__), self.file, arg) for arg in args])
            submitter.submit(file_name, command, self.walltime, options=options, queue=queue)
        else:
            s = submitter.PbsBatchSubmitter(self)
            print "SUBMITTING" #FIXME
            s.submit(standards.PythonCommand, [os.path.realpath(__file__), self.file], jobs)

    def save_and_submit(self, queue=True):
        self.save()
        self.submit(queue=queue)

    def ask(self, attribute, number=False, valid_path=False):
        try:
            name = self.name
        except AttributeError:
            name = "this Experiment"
        prompt = "Please enter the {} of {}: ".format(attribute, name)
        do = True
        while do:
            if number:
                answer = user_input.get_number(prompt)
            else:
                answer = user_input.get(prompt)
            if not valid_path or standards.is_path_component(answer):
                do = False
            else:
                disp("The {} of {} must be a valid component of a path.".format(attribute, name))
        return answer
    
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

    def run_next(self, args=None):
        self.change_status()
        self.run(args)


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
            self.number_of_designs = standards.DefaultNumberOfDesigns
        # Define the antigen and epitope.
        entire_input_model = self.get_entire_input_model()
        antigen_input_chains = self.select_antigen_input_chains(entire_input_model)
        #antigen_input_residues = [residue for chain in antigen_input_chains for residue in self.select_antigen_input_residues(chain)]
        self.antigen_chain_ids = [chain.get_id() for chain in antigen_input_chains]
        self.select_epitope_residues(antigen_input_chains)
        self.write_antigen_input_residues(entire_input_model)#, antigen_input_residues)
        self.status = 0
        self.report_directory()
        self.save_and_submit()

    def change_status(self, new_status=None):
        if new_status is None:
            new_status = self.status + 1
        self.status = new_status
        self.save()

    def run(self, args=None):
        tasks = [
            self.relax_antigen,
            self.maps_energy_batch,
            self.collect_maps_energies,
            self.select_parts_batch,
            self.select_antibodies,
            self.relax_complexes_batch,
            self.create_report,
            self.completed
        ]
        try:
            try:
                task = tasks[self.status]
            except (IndexError, TypeError):
                raise ValueError("Bad Experiment status: {}".format(self.status))
            else:
                task(args)
        except Exception as e:
            tb = traceback.format_exc()
            self.document_error(tb)

    def create_report(self, args):
        report_info = list()
        report_fields = [
            ("Total positions", len(self.positions)),
            ("Selected designs", self.select_number)
        ]
        report_info.extend([info_format(field, value) for field, value in report_fields])
        for index in range(self.select_number):
            with open(self.results_pickle_file(index)) as f:
                result = pkl.load(f)
            result_info = result.output()
            report_info.append(result_info)
        report_text = "\n".join(report_info)
        self.report_file = os.path.join(self.directory, "report.txt")
        with open(self.report_file, "w") as f:
            f.write(report_text)
        self.change_status()

    def completed(self, args):
        disp("{} has finished running. Please view the results in {}".format(self.name, self.directory))

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
        epi_vector = epi_center - all_center
        # Make a rotation matrix that rotates the epitope vector to the negative z axis.
        rot_matrix = standards.rotate_vi_to_vf(epi_vector, -standards.zAxis)
        self.epitope_zmin_file = os.path.join(self.structure_directory, "antigen_epitope_down.pdb")
        antigen_molecule.rotate(rot_matrix, self.epitope_zmin_file, in_place=True)
        # Rotate the antigen so that its z rotation angle is 0.
        #FIXME
        # Translate the antigen along the z axis to center the epitope at the origin.
        antigen_molecule.translate(standards.zAxis * np.linalg.norm(epi_vector),in_place=True)
        # Position the antigen using a grid search.
        self.position_antigen()
        
    def position_antigen(self):
        self.positions_file = os.path.join(self.temp, "positions.dat")
        with klaus.PositionAntigen(self) as x:
            pass
        self.change_status()
        self.maps_energies_all(None)

    def get_maps_part_energy_directory(self, part):
        return os.path.join(self.maps_energies_directory, part)

    def get_maps_part_energy_file_temp(self, part):
        return os.path.join(self.get_maps_part_energy_directory(part), "{}_energies_temp.dat".format(part))
    
    def get_maps_part_energy_file_finished(self, part):
        return os.path.join(self.get_maps_part_energy_directory(part), "{}_energies.dat".format(part))

    def maps_energies_all(self, args):
        """ Calculate the interacton energy between the antigen and all MAPs parts. """
        self.maps_energies_directory = os.path.join(self.temp, "maps_energies")
        if not os.path.isdir(self.maps_energies_directory):
            os.mkdir(self.maps_energies_directory)
        self.save()
        jobs = {part: self.get_maps_part_energy_file_finished(part) for part in maps.parts}
        self.submit(jobs=jobs)

    def maps_energy_batch(self, args):
        """ Calculate the interacton energy between the antigen and a batch of MAPs parts. """
        parts_file = args[2]
        with open(parts_file) as f:
            parts = f.read().split()
        for part in parts:
            self.maps_energy_single(part)
   
    def maps_energy_single(self, part):
        with klaus.MapsEnergies(self, part, self.get_maps_part_energy_directory(part), self.get_maps_part_energy_file_temp(part), self.get_maps_part_energy_file_finished(part)) as energies:
            pass
    
    def collect_maps_energies(self, args):
        maps_energies = defaultdict(list)
        for part in maps.parts:
            _file = self.get_maps_part_energy_file_finished(part)
            with open(_file) as f:
                for line in f:
                    try:
                        zAngle, x, y, z, energy = map(float, line.split())
                    except ValueError:
                        raise ValueError("Cannot parse line in MAPs energy file {}:\n{}".format(_file, line))
                    position = (zAngle, x, y, z)
                    maps_energies[position].append([maps.get_part_cdr(part), maps.get_part_number(part), energy])
        self.change_status()
        #FIXME: remove maps energies directory
        #FIXME: remove the writing of this pickle file.
        self.select_parts_all(maps_energies)

    def get_select_parts_directory(self, index):
        return os.path.join(self.select_parts_directory, "position_{}".format(index))

    def get_select_parts_energy_file(self, index):
        return os.path.join(self.get_select_parts_directory(index), "energies.pickle")

    def get_select_parts_file_temp(self, index):
        return os.path.join(self.get_select_parts_directory(index), "parts_temp.dat")

    def get_select_parts_file_finished(self, index):
        return os.path.join(self.get_select_parts_directory(index), "parts.dat")

    def select_parts_all(self, maps_energies):
        self.select_parts_directory = os.path.join(self.temp, "select_parts")
        try:
            os.mkdir(self.select_parts_directory)
        except OSError:
            pass
        self.positions = dict(enumerate(maps_energies))
        # Pickle the MAPs energies so that the select parts scripts can use them.
        for index, position in self.positions.iteritems():
            position_energies = maps_energies[position]
            try:
                os.mkdir(self.get_select_parts_directory(index))
            except OSError:
                pass
            with open(self.get_select_parts_energy_file(index), "w") as f:
                pkl.dump(position_energies, f)
        self.save()
        jobs = {index: self.get_select_parts_file_finished(index) for index in self.positions}
        self.submit(jobs=jobs)

    def select_parts_batch(self, args):
        index_file = args[2]
        with open(index_file) as f:
            indexes = f.read().split()
        for index in indexes:
            self.select_parts_single(index)
        
    def select_parts_single(self, index, solution_cuts=None):
        index = int(index)
        position = self.positions[index]
        # Load the integer cuts
        clash_cuts = maps.get_integer_cuts()
        if solution_cuts is None:
            solution_cuts = list()
        with open(self.get_select_parts_energy_file(index)) as f:
            position_energies = pkl.load(f)
        # Load the MAPS energies
        # Using cplex to get the optimal solutions. 
        # The solution is a dictionary and the key and values are partname and number respectively
        # The energy is the total interaction energy between the selected parts and antigen
        #solution, energy = OptMAVEn_selector(energies, struCuts, solutionCuts)
        solution, energy = maps.select_parts(position_energies, clash_cuts, solution_cuts)
        antibody = molecules.ProtoAntibody(solution, position, energy)
        with open(self.get_select_parts_file_finished(index), "w") as f:
            pkl.dump(antibody, f)

    def select_antibodies(self, args):
        # Cluster the antibodies based on their coordinates.
        antibodies = {chain: list() for chain in maps.light_chains}
        for design in map(self.get_select_parts_file_finished, self.positions):
            with open(design) as f:
                antibody = pkl.load(f)
            antibodies[antibody.light_chain].append(antibody)
        clusters = {chain: [sorted(cluster) for cluster in kmeans.optimal_kmeans(chain_abs)] for chain, chain_abs in antibodies.iteritems()}
        #FIXME
        #clusters = {"L": [sorted(cluster) for cluster in kmeans.optimal_kmeans(antibodies["L"])]}
        # Select the best antibodies from the clusters.
        self.highest_ranked_designs = list()
        cluster_depth = 0
        self.select_number = min(self.number_of_designs, len(self.positions))
        while len(self.highest_ranked_designs) < self.select_number:
            # Collect the best unused design from each cluster that has not been exhausted.
            cluster_heads = list()
            for light_chain, light_chain_clusters in clusters.iteritems():
                for cluster in light_chain_clusters:
                    if len(cluster) > cluster_depth:
                        cluster_heads.append(cluster[cluster_depth])
            # Add these designs to the list of best designs, in order of increasing energy.
            cluster_heads.sort()
            while len(self.highest_ranked_designs) < self.select_number and len(cluster_heads) > 0:
                self.highest_ranked_designs.append(cluster_heads.pop(0))
        self.unrelaxed_complex_directory = os.path.join(self.get_temp(), "unrelaxed_complexes")
        try:
            os.mkdir(self.unrelaxed_complex_directory)
        except OSError:
            pass
        self.unrelaxed_complex_files = [os.path.join(self.unrelaxed_complex_directory, "complex_{}.pdb".format(i)) for i in range(self.select_number)]
        [ab.to_molecule("unrelaxed", _file, self) for ab, _file in zip(self.highest_ranked_designs, self.unrelaxed_complex_files)]
        self.change_status() 
        self.relax_complexes_all()

    def relaxed_complex_directory(self):
        return os.path.join(self.directory, "antigen-antibody-complexes")

    def relaxed_complex_file(self, index):
        return os.path.join(self.relaxed_complex_directory(), "complex_{}.pdb".format(index))

    def results_directory(self, index):
        return os.path.join(self.relaxed_complex_directory(), "Result_{}".format(index))

    def results_pickle_directory(self):
        return os.path.join(self.get_temp(), "results")

    def result_name(self, index):
        return "Result_{}".format(index)

    def results_pickle_file(self, index):
        return os.path.join(self.results_pickle_directory(), "result_{}.pickle".format(index))

    def relax_complexes_all(self):
        try:
            os.mkdir(self.relaxed_complex_directory())
        except OSError:
            pass
        try:
            os.mkdir(self.results_pickle_directory())
        except OSError:
            pass
        jobs = {index: self.results_pickle_file(index) for index in range(self.select_number)}
        self.save()
        self.submit(jobs=jobs)

    def relax_complexes_batch(self, args):
        index_file = args[2]
        with open(index_file) as f:
            indexes = map(int, f.read().split())
        for index in indexes:
            self.relax_complex_single(index)

    def relax_complex_single(self, index):
        garbage = list()
        _complex = molecules.AntibodyAntigenComplex("relaxed", self.unrelaxed_complex_files[index], self)
        # Calculate interaction energy before relaxation.
        antigen, antibody = _complex.disassemble_into_chains([_complex.antigen_chains, _complex.antibody_chains])
        garbage.extend([antigen.file, antibody.file])
        with charmm.InteractionEnergy(self, [antigen, antibody]) as e:
            unrelaxed_energy = e.energy
        # Relax the complex.
        _complex.relax(relaxed_file=self.relaxed_complex_file(index), in_place=False)
        # Calculate interaction energy after relaxation.
        antigen, antibody = _complex.disassemble_into_chains([_complex.antigen_chains, _complex.antibody_chains])
        garbage.extend([antigen.file, antibody.file])
        with charmm.InteractionEnergy(self, [antigen, antibody]) as e:
            relaxed_energy = e.energy
        result = OptmavenResult(self.result_name(index), self.results_directory(index), _complex, self.highest_ranked_designs[index], unrelaxed_energy, relaxed_energy) 
        with open(self.results_pickle_file(index), "w") as f:
            pkl.dump(result, f)
        
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
        do = True
        while do:
            selected_chains = user_input.select_from_list("Please select the antigen chains: ", chains, 1, None, names=chain_ids)
            if any(chain in standards.MapsChains for chain in selected_chains):
                disp("The following are not valid names for antigen chains: {}".format(", ".join(standards.MapsChains)))
            else:
                do = False
        return selected_chains

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


class OptmavenResult(object):
    """ Stores information about one result from an Optmaven experiment. """
    def __init__(self, name, directory, molecule, proto_antibody, unrelaxed_energy, relaxed_energy):
        self.name = name
        self.directory = directory
        self.molecule = molecule
        self.proto_antibody = proto_antibody
        self.unrelaxed_energy = unrelaxed_energy
        self.relaxed_energy = relaxed_energy

    def output(self):
        try:
            os.mkdir(self.directory)
        except OSError:
            pass
        info = list()
        info.append(info_format("Name", self.name))
        # Write the molecule.
        self.molecule_file = os.path.join(self.directory, "{}.pdb".format(self.name))
        shutil.copyfile(self.molecule.file, self.molecule_file)
        info.append(info_format("PDB file", self.molecule_file)) 
        # Get the sequence and output it as a fasta.
        self.fasta_file = os.path.join(self.directory, "{}.fasta".format(self.name))
        SeqIO.write(SeqIO.parse(self.molecule_file, "pdb-atom"), self.fasta_file, "fasta")
        info.append("FASTA sequences")
        with open(self.fasta_file) as f:
            info.append(f.read())
        info.append("MAPs parts")
        [info.append(info_format(field, value)) for field, value in self.proto_antibody.get_namesake_parts().items()]
        info.append("Position")
        [info.append(info_format(field, value)) for field, value in self.proto_antibody.get_labeled_position().items()]
        info.append(info_format("Unrelaxed energy", self.unrelaxed_energy))
        info.append(info_format("Relaxed energy", self.relaxed_energy))
        return "\n".join(info)


def info_format(field, value):
    return "{}\t{}".format(field, value)
        

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise OSError("Usage: python experiments.py /path/to/experiment.pickle")
    experiment_file = sys.argv[1]
    with open(experiment_file) as f:
        experiment = pkl.load(f)
    experiment.run(sys.argv)
