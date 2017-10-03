""" This module defines the Experiment class of OptMAVEn. """

from collections import defaultdict, OrderedDict
import cPickle as pkl
import csv
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

import benchmarking
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
        self.purpose = "Initialization"
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
        self.structure_directory = os.path.join(self.directory, "structures")
        os.mkdir(self.structure_directory)
        self.benchmarking = user_input.get_yn("Would you like to turn on benchmarking? ")
        if self.benchmarking:
            self.benchmark_directory = os.path.join(self.get_temp(), "benchmarking")
            os.mkdir(self.benchmark_directory)
            self.add_drive_usage()

    def save(self):
        with open(self.file, "w") as f:
            pkl.dump(self, f)

    def submit(self, args=None, jobs=None, options=None, queue=True):
        if args is None:
            args = [""]
        if jobs is None:
            handle, file_name = tempfile.mkstemp(prefix="{}_".format(self.name_contig), suffix=".sh", dir=self.temp)
            os.close(handle)
            command = "{} {} {} {}".format(standards.PythonCommand, os.path.realpath(__file__), self.file, " ".join(map(str, args)))
            if self.benchmarking:
                time_file = self.make_time_file()
                self.add_time_file(time_file)
            else:
                time_file = None
            submitter.submit(file_name, command, self.walltime, options=options, queue=queue, purpose=self.purpose, time_file=time_file)
        else:
            s = submitter.PbsBatchSubmitter(self)
            s.submit(standards.PythonCommand, [os.path.realpath(__file__), self.file], jobs)

    def save_and_submit(self, queue=True):
        self.save()
        self.submit(queue=queue)

    def add_benchmark(self, task):
        handle, file_name = tempfile.mkstemp(dir=self.benchmark_directory, suffix=".pickle")
        os.close(handle)
        with open(file_name, "w") as f:
            pkl.dump(task, f)

    def add_time_file(self, _file, status_offset=0):
        task, purpose = self.get_task(self.status + status_offset)
        task = benchmarking.TimeFile(_file, purpose)
        self.add_benchmark(task)

    def add_drive_usage(self):
        du = benchmarking.DriveUsage(self)
        self.add_benchmark(du)
        return du

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

    def make_time_file(self):
        handle, name = tempfile.mkstemp(dir=self.get_temp(), prefix="time_", suffix=".txt")
        os.close(handle)
        return name

    def document_error(self, message):
        try:
            with open(self.errors, "a") as f:
                f.write(str(message))
        except OSError as e:
            error_file = "Optmaven_Experiment_errors.txt"
            with open(error_file, "a") as f:
                f.write("{}\n{}".format(message, e.message))

    def run(self, args=None):
        try:
            try:
                task, self.purpose = self.get_task(self.status)
            except (IndexError, TypeError):
                raise ValueError("Bad Experiment status: {}".format(self.status))
            else:
                task(args)
        except Exception as e:
            tb = traceback.format_exc()
            self.document_error(tb)

    def run_next(self, args=None):
        self.change_status()
        self.run(args)

    def get_molecule(self, name, prompt):
        do = True
        while do:
            molecule_file = user_input.get_file(prompt, standards.PDBDirectory, fetch_pdb=True)
            try:
                molecule = molecules.Molecule(name, molecule_file, self, exclude_hetero_ask=True)
                molecule.get_structure()
            except Exception as error:
                disp("There was an error with the PDB import:\n{}".format(error.message))
            else:
                do = False
        return molecule

    def safe_rmtree(self, directory):
        if standards.is_subdirectory(directory, self.get_temp()):
            standards.safe_rmtree(directory)
        else:
            raise OSError("{} cannot remove directory trees outside of {}".format(self.name, self.get_temp()))

    def purge_temp(self):
        standards.safe_rmtree(self.temp)

    def completed(self, args):
        disp("{} has finished running. Please view the results in {}".format(self.name, self.directory))


class InteractionExperiment(Experiment):
    """ This Experiment simply calculates the interaction energy between two groups of chains, performs a relaxation, then calculates the energy again. """
    def __init__(self):
        Experiment.__init__(self)
        self.topology_files = [standards.DefaultTopologyFile]
        self.parameter_files = [standards.DefaultParameterFile]
        self.solvation_files = [standards.DefaultSolvationFile]
        self.charmm_energy_terms = standards.DefaultCharmmEnergyTerms
        self.charmm_iterations = standards.DefaultCharmmIterations
        self.walltime = standards.DefaultWalltime
        self.molecule = self.get_molecule("input", "Please enter the name of the molecule for {}: ".format(self.name))
        self.molecule_file = self.molecule.file
        # Select chains.
        self.chain_ids = [chain.get_id() for chain in self.molecule.get_chains()]
        chains_left = list(self.chain_ids)
        self.chain_groups = list()
        for i in [1, 2]:
            group = user_input.select_from_list("Please select the chain(s) for group {}: ".format(i), chains_left, min_number=1, max_number=(len(chains_left) + i - 2))
            for chain_id in group:
                chains_left.pop(chains_left.index(chain_id))
            self.chain_groups.append(group)
        self.status = 0
        self.save_and_submit()

    def get_task(self, status):
        return [(self.relaxation, "Relaxation")][status]

    def relaxation(self, args):
        garbage = list()
        # Calculate interaction energy before relaxation.
        group1, group2 = self.molecule.disassemble_into_chains(self.chain_groups)
        self.energies = dict()
        with charmm.InteractionEnergy(self, [group1, group2]) as e:
            self.energies["Unrelaxed"] = e.energy
            if self.benchmarking:
                self.add_drive_usage()
        # Relax the complex.
        self.relaxed_file = os.path.join(self.structure_directory, "relaxed.pdb")
        self.molecule.relax(relaxed_file=self.relaxed_file, in_place=True)
        # Calculate interaction energy after relaxation.
        group1, group2 = self.molecule.disassemble_into_chains(self.chain_groups)
        garbage.extend([group1.file, group2.file])
        with charmm.InteractionEnergy(self, [group1, group2]) as e:
            self.energies["Relaxed"] = e.energy
            if self.benchmarking:
                self.add_drive_usage()
        for fn in garbage:
            try:
                os.remove(fn)
            except OSError:
                pass
        self.create_report()
        self.safe_rmtree(self.get_temp())
    
    def create_report(self):
        summary = [
            ("Experiment name", self.name),
            ("Molecule input file", self.molecule_file),
            ("Group 1 chains", ", ".join(self.chain_groups[0])),
            ("Group 2 chains", "; ".join(self.chain_groups[1]))
        ]
        self.summary_file = os.path.join(self.directory, "Summary.txt")
        with open(self.summary_file, "w") as f:
            f.write("\n".join(["{}: {}".format(field, value) for field, value in summary]))
        self.results_file = os.path.join(self.directory, "Results.csv")
        with open(self.results_file, "w") as f:
            writer = csv.DictWriter(f, ["Unrelaxed", "Relaxed"])
            writer.writeheader()
            writer.writerow(self.energies)


class ResidueContactExperiment(Experiment):
    def __init__(self):
        Experiment.__init__(self)
        self.molecule = self.get_molecule("input", "Please enter the name of the molecule for {}: ".format(self.name))
        self.chain_ids = [chain.get_id() for chain in self.molecule.get_chains()]
        chains_left = list(self.chain_ids)
        self.chain_groups = list()
        self.group_numbers = [1, 2]
        for i in self.group_numbers:
            group = user_input.select_from_list("Please select the chain(s) for group {}: ".format(i), chains_left, min_number=1, max_number=(len(chains_left) + i - 2))
            for chain_id in group:
                chains_left.pop(chains_left.index(chain_id))
            self.chain_groups.append(group)
        g1, g2 = self.chain_groups
        group_from_chain = dict()
        for group_number, group in zip(self.group_numbers, self.chain_groups):
            for chain in group:
                group_from_chain[chain] = group_number
        self.cutoff = self.ask("cutoff distance", number=True)
        contacts = self.molecule.interchain_residue_contacts(g1, g2, self.cutoff)
        # Restructure the contact map so that it looks like {group_number: {chain: [residue, residue, ...], chain: ...}, group_number: ...}
        self.contacts_from_group = defaultdict(lambda: defaultdict(list))
        for chain_pair, chain_pair_contacts in contacts.iteritems():
            for chain in chain_pair:
                group = group_from_chain[chain]
                residues = self.contacts_from_group[group][chain]
                for contact in chain_pair_contacts:
                    residue = contact[chain]
                    if residue not in residues:
                        residues.append(residue)
        # Sort the residues by their ids, and then convert them into residue codes.
        for group, chains in self.contacts_from_group.iteritems():
            for chain, residues in chains.iteritems():
                residues.sort(key=Residue.Residue.get_id)
                for i in range(len(residues)):
                    residues[i] = molecules.residue_code(residues[i])
        self.create_report()

    def create_report(self):
        summary = [
            ("Experiment name", self.name),
            ("Molecule input file", self.molecule.file),
            ("Group 1 chains", ", ".join(self.chain_groups[0])),
            ("Group 2 chains", "; ".join(self.chain_groups[1])),
            ("Cutoff distnce", self.cutoff)
        ]
        self.summary_file = os.path.join(self.directory, "Summary.txt")
        with open(self.summary_file, "w") as f:
            f.write("\n".join(["{}: {}".format(field, value) for field, value in summary]))
        self.results_file = os.path.join(self.directory, "Results.csv")
        with open(self.results_file, "w") as f:
            fields = ["Group", "Chain", "Residues"]
            writer = csv.DictWriter(f, fields)
            writer.writeheader()
            for group in self.group_numbers:
                for chain in self.chain_ids:
                    if chain in self.contacts_from_group[group]:
                        info = {
                            "Group": str(group),
                            "Chain": chain,
                            "Residues": ",".join(self.contacts_from_group[group][chain])
                        }
                        writer.writerow(info)


class TransformMoleculeExperiment(Experiment):
    """ This is a simple Experiment that just translates and rotates a molecule. """
    def __init__(self):
        Experiment.__init__(self)
        self.molecule = self.get_molecule("input", "Please enter the name of the molecule for {}: ".format(self.name))
        # Get the translation vector.
        do = True
        while do:
            try:
                translation_vector = np.array(map(float, self.ask("translation vector").split()))
            except ValueError:
                pass
            else:
                do = False
        # Get the rotation angle.
        rotation_degrees = self.ask("rotation angle (in degrees)", number=True)
        rotation_degrees -= np.floor(rotation_degrees / 360.0)
        # Get the rotation axis.
        do = not np.isclose(rotation_degrees, 0.0)
        rotation_axis = None
        while do:
            try:
                rotation_axis = np.array(map(float, self.ask("rotation axis").split()))
            except ValueError:
                pass
            else:
                do = False
        self.output_file = os.path.join(self.structure_directory, self.ask("output molecule", valid_path=True))
        self.molecule.translate(translation_vector, file_name=self.output_file, in_place=True)
        if rotation_axis is not None:
            self.molecule.rotate(standards.rotate_axis_angle(rotation_axis, rotation_degrees, degrees=True), in_place=True)
        self.completed(None)


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
        self.antigen_chain_ids = [chain.get_id() for chain in antigen_input_chains]
        self.select_epitope_residues(antigen_input_chains)
        self.write_antigen_input_residues(entire_input_model)
        self.status = 0
        self.report_directory()
        self.save_and_submit()

    def change_status(self, new_status=None):
        if new_status is None:
            new_status = self.status + 1
        self.status = new_status
        self.save()

    def get_task(self, status):
        return [
            (self.relax_antigen, "Antigen relaxation, positioning, and grid search"),
            (self.maps_energy_batch, "MAPs energy calculations"),
            (self.collect_maps_energies, "Prepare for MILP design"),
            (self.select_parts_batch, "MILP design"),
            (self.select_antibodies, "Select designs"),
            (self.relax_complexes_batch, "Relaxing designs"),
            (self.create_report, "Creating report"),
            (self.completed, "Completed")
        ][status]

    def create_report(self, args):
        summary = [
            ("Experiment name", self.name),
            ("Antigen input file", self.entire_input_file),
            ("Antigen chains", ", ".join(self.antigen_chain_ids)),
            ("Epitope residues", "; ".join([", ".join(["{}-{}".format(chain, "".join(map(str, resid))) for resid in residues]) for chain, residues in self.epitope_residue_ids.items()])),
            ("Total positions", len(self.positions)),
            ("Selected designs", self.select_number)
        ]
        self.summary_report = os.path.join(self.directory, "Summary.txt")
        with open(self.summary_report, "w") as f:
            f.write("\n".join(["{}: {}".format(field, value) for field, value in summary]))
        result_report_info = list()
        result_report_fields = list()
        for index in range(self.select_number):
            with open(self.results_pickle_file(index)) as f:
                result = pkl.load(f)
            result_info = result.output()
            result_report_info.append(result_info)
        for result in result_report_info:
            for field in result:
                if field not in result_report_fields:
                    result_report_fields.append(field)
        self.result_report = os.path.join(self.directory, "Results.csv")
        with open(self.result_report, "w") as f:
            writer = csv.DictWriter(f, fieldnames=result_report_fields)
            writer.writeheader()
            for result in result_report_info:
                writer.writerow(result)
        if self.benchmarking:
            self.benchmarking_file = os.path.join(self.directory, "Benchmarking.csv")
            with open(self.benchmarking_file, "w") as f:
                writer = csv.DictWriter(f, fieldnames=standards.BenchmarkingFields)
                totals = {field: 0.0 for field in standards.UnixTimeCodes.keys() + ["Drive Usage"]}
                totals["Type"] = "Total"
                writer.writeheader()
                benchmark_files = [os.path.join(self.benchmark_directory, fn) for fn in os.listdir(self.benchmark_directory)]
                benchmarks = list()
                for fn in benchmark_files:
                    with open(fn) as f:
                        benchmarks.append(pkl.load(f))
                benchmarks.sort(key=benchmarking.Task.get_time_stamp)
                tried_current_time_file = False
                for benchmark in benchmarks:
                    try:
                        d = benchmark.to_dict()
                    except benchmarking.BlankTimeFileError:
                        # Expect to find exactly one blank time file: the one measuring the current process.
                        if tried_current_time_file:
                            raise ValueError("Found two blank time files.")
                        tried_current_time_file = True
                    else:
                        writer.writerow(d)
                        if isinstance(benchmark, benchmarking.Time):
                            for field in standards.UnixTimeCodes:
                                totals[field] += d[field]
                        elif isinstance(benchmark, benchmarking.DriveUsage):
                            totals["Drive Usage"] = max(d["Drive Usage"], totals["Drive Usage"])
                # Remove the temporary directory, then add the final drive usage.
                #FIXME self.safe_rmtree(self.temp)
                du = self.add_drive_usage().to_dict()
                writer.writerow(du)
                totals["Drive Usage"] = max(du["Drive Usage"], totals["Drive Usage"])
                writer.writerow(totals)
        else:
            pass
            #FIXME: remove temp directory
            #self.safe_rmtree(self.temp) 
        self.change_status()

    def relax_antigen(self, args):
        antigen_molecule = molecules.Molecule(self.antigen_input_name, self.antigen_input_file, self)
        self.antigen_relaxed_name = "relaxed_antigen"
        self.antigen_relaxed_file = os.path.join(self.structure_directory, "antigen_relaxed.pdb")
        antigen_molecule.relax(self.antigen_relaxed_file)
        self.minimize_epitope_z_coordinates(antigen_molecule)

    def minimize_epitope_z_coordinates(self, antigen_molecule, grid_search=True):
        self.epitope_zmin_file = self.antigen_relaxed_file
        antigen_molecule.position_antigen(0, 0, 0, 0, in_place=True)
        # Position the antigen using a grid search.
        if grid_search:
            self.grid_search()
        
    def grid_search(self):
        self.positions_file = os.path.join(self.temp, "positions.dat")
        with klaus.GridSearch(self) as x:
            if self.benchmarking:
                self.add_drive_usage()
        self.change_status()
        self.maps_energies_all(None)

    def get_maps_part_energy_directory_finished(self, part):
        return os.path.join(self.maps_energies_directory, part)

    def get_maps_part_energy_file_finished(self, part):
        return os.path.join(self.get_maps_part_energy_directory_finished(part), "{}_energies.dat".format(part))

    def maps_energies_all(self, args):
        """ Calculate the interacton energy between the antigen and all MAPs parts. """
        self.maps_energies_directory = os.path.join(self.temp, "maps_energies")
        try:
            os.mkdir(self.maps_energies_directory)
        except OSError:
            pass
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
        with klaus.MapsEnergies(self, part, self.get_maps_part_energy_file_finished(part)) as energies:
            if self.benchmarking:
                self.add_drive_usage()
    
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
        self.select_parts_all(maps_energies)

    def get_select_parts_directory(self, index):
        return os.path.join(self.select_parts_directory, "position_{}".format(index))

    def get_select_parts_energy_file(self, index):
        return os.path.join(self.get_select_parts_directory(index), "energies.pickle")

    def get_select_parts_file_finished(self, index):
        return os.path.join(self.get_select_parts_directory(index), "parts.pickle")

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
        if self.benchmarking:
            self.add_drive_usage()
        self.save()
        # Remove maps energies directory: all energies are saved in pickle files.
        self.safe_rmtree(self.maps_energies_directory)
        if self.benchmarking:
            self.add_drive_usage()
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
        if self.benchmarking:
            self.add_drive_usage()

    def select_antibodies(self, args):
        # Cluster the antibodies based on their coordinates.
        antibodies = defaultdict(list)
        for design in map(self.get_select_parts_file_finished, self.positions):
            with open(design) as f:
                antibody = pkl.load(f)
            antibodies[antibody.light_chain].append(antibody)
        clusters = {chain: [sorted(cluster) for cluster in kmeans.optimal_kmeans(chain_abs)] for chain, chain_abs in antibodies.iteritems()}
        # Select the best antibodies from the clusters.
        self.highest_ranked_designs = list()
        cluster_depth = 0
        self.select_number = min(self.number_of_designs, len(self.positions))
        while len(self.highest_ranked_designs) < self.select_number:
            # Collect the best unused design from each cluster that has not been exhausted.
            cluster_heads = list()
            for light_chain, light_chain_clusters in clusters.iteritems():
                for i, cluster in enumerate(light_chain_clusters):
                    if len(cluster) > cluster_depth:
                        antibody = cluster[cluster_depth]
                        antibody.set_cluster(i)
                        cluster_heads.append(antibody)
            # Add these designs to the list of best designs, in order of increasing energy.
            cluster_heads.sort()
            while len(self.highest_ranked_designs) < self.select_number and len(cluster_heads) > 0:
                self.highest_ranked_designs.append(cluster_heads.pop(0))
            cluster_depth += 1
        self.unrelaxed_complex_directory = os.path.join(self.get_temp(), "unrelaxed_complexes")
        try:
            os.mkdir(self.unrelaxed_complex_directory)
        except OSError:
            pass
        self.unrelaxed_complex_files = [os.path.join(self.unrelaxed_complex_directory, "complex_{}.pdb".format(i)) for i in range(self.select_number)]
        [ab.to_molecule("unrelaxed", _file, self) for ab, _file in zip(self.highest_ranked_designs, self.unrelaxed_complex_files)]
        if self.benchmarking:
            self.add_drive_usage()
        # Remove the select parts directory; it is no longer needed.
        self.safe_rmtree(self.select_parts_directory)
        if self.benchmarking:
            self.add_drive_usage()
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
            if self.benchmarking:
                self.add_drive_usage()
        # Relax the complex.
        _complex.relax(relaxed_file=self.relaxed_complex_file(index), in_place=True)
        # Calculate interaction energy after relaxation.
        antigen, antibody = _complex.disassemble_into_chains([_complex.antigen_chains, _complex.antibody_chains])
        garbage.extend([antigen.file, antibody.file])
        with charmm.InteractionEnergy(self, [antigen, antibody]) as e:
            relaxed_energy = e.energy
            if self.benchmarking:
                self.add_drive_usage()
        result = OptmavenResult(self.result_name(index), self.results_directory(index), _complex, self.highest_ranked_designs[index], unrelaxed_energy, relaxed_energy) 
        with open(self.results_pickle_file(index), "w") as f:
            pkl.dump(result, f)
        
    def get_entire_input_model(self):
        # Select the antigen file.
        self.entire_input_name = "entire_input_file"
        do = True
        while do:
            self.entire_input_file = user_input.get_file("Please name the file containing the antigen: ", standards.PDBDirectory, fetch_pdb=True)
            try:
                entire_input_model = molecules.Molecule(self.entire_input_name, self.entire_input_file, self, exclude_hetero_ask=True).get_model()
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
            if any(chain.get_id() in standards.MapsChains for chain in selected_chains):
                disp("The following are not valid names for antigen chains: {}".format(", ".join(standards.MapsChains)))
            else:
                do = False
        return selected_chains

    def select_epitope_residues(self, chains):
        self.epitope_residue_ids = dict()
        for chain in chains:
            self.epitope_residue_ids[chain.get_id()] = user_input.select_from_list("Please select the epitope residues from chain {}: ".format(chain.get_id()), [residue.get_id() for residue in chain], 1, None, map(molecules.residue_code, chain))
        #return antigen_input_residues
    
    def write_antigen_input_residues(self, entire_input_structure):
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
        info = OrderedDict()
        info["Result"] = self.name
        # Write the molecule.
        self.molecule_file = os.path.join(self.directory, "{}.pdb".format(self.name))
        shutil.copyfile(self.molecule.file, self.molecule_file)
        info["PDB file"] = self.molecule_file
        # Get the sequence and output it as a fasta.
        self.fasta_file = os.path.join(self.directory, "{}.fasta".format(self.name))
        records = list(SeqIO.parse(self.molecule_file, "pdb-atom"))
        SeqIO.write(records, self.fasta_file, "fasta")
        for record in records:
            chain = record.id.split(":")[1]
            info[chain] = str(record.seq)
        for cdr, part in self.proto_antibody.get_namesake_parts().items():
            info[cdr] = part
        for dimension, coord in self.proto_antibody.get_labeled_position().items():
            info[dimension] = coord
        info["Cluster"] = self.proto_antibody.cluster
        info["MILP energy (kcal/mol)"] = self.proto_antibody.energy
        info["Unrelaxed energy (kcal/mol)"] = self.unrelaxed_energy
        info["Relaxed energy (kcal/mol)"] = self.relaxed_energy
        return info


class CreateAntigenAntibodyComplexExperiment(OptmavenExperiment):
    def __init__(self):
        Experiment.__init__(self)
        entire_input_model = OptmavenExperiment.get_entire_input_model(self)
        antigen_input_chains = OptmavenExperiment.select_antigen_input_chains(self, entire_input_model)
        self.antigen_chain_ids = [chain.get_id() for chain in antigen_input_chains]
        OptmavenExperiment.select_epitope_residues(self, antigen_input_chains)
        # Position the antigen.
        do = True
        while do:
            try:
                zAngle, x, y, z = map(float, self.ask("antigen position").split())
            except ValueError:
                pass
            else:
                do = False
        position = (zAngle, x, y, z)
        # Get the chain loci.
        heavy = user_input.select_one_from_list("Please specify the heavy chain: ", standards.MapsHeavyChains)
        light = user_input.select_one_from_list("Please specify the light chain: ", standards.MapsLightChains)
        # Get the MAPs parts.
        parts = dict()
        for cdr in standards.MapsCdrs:
            chain = maps.get_chain(cdr)
            if chain not in [heavy, light]:
                continue
            do = True
            while do:
                try:
                    number = int(self.ask(cdr))
                    part = maps.join(cdr, number)
                except ValueError as e:
                    disp(e.message)
                else:
                    do = False
            parts[cdr] = number
        energy = None
        output_file = user_input.get_file("Please specify the output file: ", self.directory, new_file=True)
        # Write the antigen input file.
        OptmavenExperiment.write_antigen_input_residues(self, entire_input_model)
        antigen_molecule = molecules.Molecule(self.antigen_input_name, self.antigen_input_file, self) 
        # Center and point the epitope downward.
        OptmavenExperiment.minimize_epitope_z_coordinates(self, antigen_molecule, grid_search=False)
        # Assemble the complex.
        self.proto_antibody = molecules.ProtoAntibody(parts, position, energy)
        self.proto_antibody.to_molecule("complex", output_file, self)
        clear()
        self.report_directory()


def info_format(field, value):
    return "{}\t{}".format(field, value)
        

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise OSError("Usage: python experiments.py /path/to/experiment.pickle")
    experiment_file = sys.argv[1]
    with open(experiment_file) as f:
        experiment = pkl.load(f)
    experiment.run(sys.argv)
