__doc__ = """ This module is the interface between IPRO and Klaus Schulten
programs VMD and NAMD and contains a built-in function for performing each task
 with VMD and NAMD. """

from collections import OrderedDict 
import itertools
import os
import shutil
import subprocess
import tempfile
import time

import numpy as np

import maps
import molecules
import charmm
import standards
import submitter

# Define how to run VMD.
vmd_command = standards.VmdCommand


class VmdProc(object):
    """ This class is the base class for managing all VMD procedures. """
    def __init__(self, experiment, directory=None):
        self.experiment = experiment
        self.directory = directory
        self.has_run = False

    def make_command(self):
        """ Create a VMD command. """
        command = [standards.VmdCommand, standards.VmdDisp, standards.VmdNoDisp, standards.VmdExec, os.path.join(standards.SourceDirectory, self.script)]
        try:
            command.append("{} {}".format(standards.VmdMolecules, " ".join(map(str, self.molecules))))
        except AttributeError:
            pass
        try:
            command.append("{} {}".format(standards.VmdFrames, " ".join(map(str, self.frames))))
        except AttributeError:
            pass
        try:
            command.append("{} {}".format(standards.VmdArgs, " ".join(map(str, self.args))))
        except AttributeError:
            pass
        return command

    def vmd(self):
        arguments = self.make_command()
        command = " ".join(arguments)
        status = os.system(command)
        if status != 0:
            raise RuntimeError("Running VMD with the following command has failed:\n{}".format(command))

    def collect_garbage(self):
        if not self.experiment.args.keeptemp:
            try:
                os.remove(self.vmd_functions_file)
            except OSError:
                pass
            self.experiment.safe_rmtree(self.directory)

    def write_epitope_file(self):
        self.epitope_file = os.path.join(self.directory, "epitope.txt")
        with open(self.epitope_file, "w") as f:
            f.write(atomselect_residues(self.experiment.epitope_residue_ids))

    def __enter__(self):
        if self.directory is None:
            self.directory = tempfile.mkdtemp(prefix="vmd_", dir=self.experiment.temp)
        else:
            os.mkdir(self.directory)
        # Make a link to the VMD functions file.
        self.vmd_functions_file = os.path.join(self.directory, os.path.basename(standards.VmdFunctions))
        os.symlink(standards.VmdFunctions, self.vmd_functions_file)
        self.previous_directory = os.getcwd()
        os.chdir(self.directory)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.previous_directory)
        self.collect_garbage()


class PositionAntigen(VmdProc):
    """ Move an antigen to a new position. """
    def __init__(self, experiment, input_file, output_file, zAngle, x, y, z):
        VmdProc.__init__(self, experiment)
        self.molecules = [input_file]
        self._args = [experiment.antigen_chain_ids[0], output_file, zAngle, x, y, z]

    def __enter__(self):
        VmdProc.__enter__(self)
        self.write_epitope_file()
        self.args = [self.epitope_file] + self._args
        self.script = "position_antigen.tcl"
        self.vmd()
        return self


class GridSearch(VmdProc):
    """ Generate non-clashing antigen positions using a grid search. """
    def __init__(self, experiment):
        VmdProc.__init__(self, experiment)

    def write_grid_file(self):
        self.grid_file = os.path.join(self.directory, "grid.dat")
        with open(self.grid_file, "w") as f:
            for dimension, levels in [
                (standards.xLabel, self.experiment.grid_x),
                (standards.yLabel, self.experiment.grid_y),
                (standards.zLabel, self.experiment.grid_z),
                (standards.zAngleLabel, self.experiment.grid_zAngle)
            ]:
                f.write("{}: {}\n".format(dimension, " ".join(map(str, levels))))

    def __enter__(self):
        VmdProc.__enter__(self)
        self.write_grid_file()
        self.write_epitope_file()
        self.positions_file = os.path.join(self.directory, "positions.dat")
        self.molecules = [self.experiment.epitope_zmin_file, standards.ScaffoldAntibodies["H"], standards.ScaffoldAntibodies["K"]]
        self.args = [self.grid_file, self.epitope_file, self.experiment.antigen_chain_ids[0], self.positions_file, self.experiment.clash_cutoff]
        self.script = "grid_search.tcl"
        self.vmd()
        shutil.move(self.positions_file, self.experiment.positions_file)
        return self
        

class MergeAntigenMaps(VmdProc):
    """ Calculate the interaction energies between the antigen and a MAPs part in every position. """
    def __init__(self, experiment, maps_part, destination_directory):
        VmdProc.__init__(self, experiment)
        self.part = maps_part
        self.destination_directory = destination_directory

    def __enter__(self):
        VmdProc.__enter__(self)
        self.prefix = os.path.join(self.directory, "merged")
        pdb = "{}.pdb".format(self.prefix)
        psf = "{}.psf".format(self.prefix)
        self.args = [self.experiment.epitope_zmin_file, maps.parts[self.part], self.prefix]
        self.args.extend(self.experiment.topology_files)
        self.script = "merge_antigen_part.tcl"
        self.vmd()
        shutil.move(pdb, self.destination_directory)
        shutil.move(psf, self.destination_directory)
        self.pdb = os.path.join(self.destination_directory, os.path.basename(pdb))
        self.psf = os.path.join(self.destination_directory, os.path.basename(psf))
        return self


class MapsEnergies(VmdProc):
    """ Calculate the interaction energies between the antigen and a MAPs part in every position. """
    def __init__(self, experiment, maps_part, energy_finished):
        VmdProc.__init__(self, experiment)
        self.part = maps_part
        self.energy_finished = energy_finished

    def __enter__(self):
        VmdProc.__enter__(self)
        with MergeAntigenMaps(self.experiment, self.part, self.directory) as merge:
            self.pdb = merge.pdb
            self.psf = merge.psf
        energy_temp = os.path.join(self.directory, "energies_temp.dat")
        self.write_epitope_file()
        self.args = [self.psf, self.pdb, self.experiment.positions_file, self.epitope_file, self.experiment.antigen_chain_ids[0], energy_temp]
        self.args.extend(self.experiment.parameter_files)
        self.script = "interaction_energies.tcl"
        self.vmd()
        try:
            os.makedirs(os.path.dirname(self.energy_finished))
        except OSError:
            pass
        shutil.move(energy_temp, self.energy_finished)
        return self


class CreateAntibody(VmdProc):
    """ Merge six MAPs parts into an antibody. """
    def __init__(self, experiment, maps_parts, output_file):
        VmdProc.__init__(self, experiment)
        self.maps_parts = dict()
        self.args = [output_file]
        for part in maps_parts:
            cdr, number = maps.split(part)
            self.maps_parts[maps.translate_chain_namesake(cdr)] = maps.parts[part]
        if sorted(self.maps_parts) != sorted(standards.MapsNamesakeCdrs):
            raise ValueError("An antibody needs one of each CDR, not: {}".format(", ".join(self.maps_parts)))
    
    def __enter__(self):
        VmdProc.__enter__(self)
        self.molecules = [self.maps_parts[cdr] for cdr in standards.MapsNamesakeCdrs]
        self.script = "create_antibody.tcl"
        self.vmd()
        return self


def MAPs_part_clashes(cut_file, clash_cutoff=1.0):
    """ This function is only used for the pre-processing step of
    determining which pairs of MAPs parts clash sterically and
    developing integer cuts. """
    parts = MAPs.get_parts()
    # Choose all combinations of categories.
    # Do not check for clashes between parts of the same category.
    finished = True
    text = ""
    for cat1, cat2 in itertools.combinations(parts, 2):
        # Also, do not check for clashes between one kappa and one lambda part.
        if "".join(sorted(cat1[0] + cat2[0])) != "KL":
            result = "{}_{}_{}.txt".format(cut_file, cat1, cat2)
            if os.path.isfile(result):
                with open(result) as f:
                    text += f.read().strip() + "\n"
            else:
                script = "{}_{}_{}.sh".format(cut_file, cat1, cat2)
                command = "python {}/databases/MAPs/Make_Integer_Cuts.py {} {} {} {}".format(STANDARDS.InstallFolder, cat1, cat2, result, clash_cutoff)
                SUBMITTER.experiment_script(None, script, command, self_destruct=True)
                finished = False
    if finished:
        with open(cut_file + ".txt", "w") as f:
            f.write(text)
        for cat1, cat2 in itertools.combinations(parts, 2):
            if "".join(sorted(cat1[0] + cat2[0])) != "KL":
                result = "{}_{}_{}.txt".format(cut_file, cat1, cat2)
                try:
                    os.remove(result)
                except OSError:
                    pass
    

def MAPs_part_category_clashes(cat1, cat2, cut_file, clash_cutoff=1.0):
    parts = MAPs.get_parts()
    f = open(cut_file, "w")
    f.close()
    for part1, part2 in itertools.product(parts[cat1], parts[cat2]):
        pdbs = map(MAPs.get_part_file, [part1, part2])
        handle, clash_file = tempfile.mkstemp()
        run_vmd_script("count_clashes.tcl", molecules=pdbs,
                args=[clash_file, clash_cutoff])
        with open(clash_file) as f2:
            clashes = int(f2.read())
        if clashes > 0:
            with open(cut_file, "a") as f:
                f.write("{} {}\n".format(part1, part2).replace("_", " "))
        os.remove(clash_file)


def atomselect_residues(residue_ids):
    return " or ".join(["(chain {} and resid {})".format(c_id, "{}{}{}".format(*r_id)) for c_id, r_ids in residue_ids.items() for r_id in r_ids])
