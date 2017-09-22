from collections import OrderedDict
import os
import shutil
import string
import subprocess
import sys
import tempfile

from Bio.PDB import Selection

import molecules
import standards


class CharmmProc(object):
    """ This class is the base class for managing all CHARMM procedures. """
    def __init__(self, experiment, molecule_list):
        self.experiment = experiment
        self.molecules = molecule_list
        self.resnums = OrderedDict()
        for molecule in self.molecules:
            molecule_resnums = OrderedDict()
            for model in molecule.get_structure():
                model_resnums = OrderedDict()
                for chain in model:
                    model_resnums[chain.get_id()] = [residue.get_id()[1] for residue in chain]
                molecule_resnums[model.get_id()] = model_resnums
            self.resnums.update(molecule_resnums)
        self.has_run = False

    def begin_script(self):
        self.lines = introduction()

    def load_input_files(self):
        """ Load the Topology and Parameter input files in a CHARMM script. """
        # Loop through topology and parameter data.
        input_files = [
            {"title": "topology" , "files": self.experiment.topology_files,  "cmd": "rtf" },
            {"title": "parameter", "files": self.experiment.parameter_files, "cmd": "para"}
        ]
        for file_type in input_files:
            self.lines.append("! Load the {} file(s)".format(file_type["title"]))
            for count, file_ in enumerate(file_type["files"]):
                # Make a symbolic link in the current directory to the file.
                if not os.path.isfile(file_):
                    raise OSError("File does not exist: {}".format(file_))
                file_name = os.path.basename(file_).lower()
                symlink = os.path.join(self.directory, file_name)
                os.symlink(file_, symlink)
                self.lines.append("open read unit 10 form name {} card".format(file_name))
                # If there is more than one file, indicate that the subsequent files are being appended.
                append = " append" * int(count > 0)
                self.lines.append("read {} card unit 10{}".format(file_type["cmd"], append))
                self.lines.append("close unit 10")

    def load_molecules(self):
        """ Create text to load Molecules in a CHARMM script. """
        # Merge all of the molecules.
        merged_structure = molecules.merge(self.molecules)
        # Convert the molecule to CHARMM format.
        charmm_format_name = "charmm_format"
        handle, charmm_format_file = tempfile.mkstemp(suffix=".pdb", dir=self.directory)
        merged_molecule = molecules.Molecule(charmm_format_name, charmm_format_file, self.experiment)
	chains, files = merged_molecule.disassemble_for_CHARMM(structure=merged_structure, directory=self.directory)
        try:
            os.remove(charmm_format_file)
        except OSError:
            pass
        self.chain_ids = list()
        for chain, _file in zip(chains, files):
            _id = chain.get_id()
            self.chain_ids.append(_id)
            segment = make_segment_name(_id)
            base_name = os.path.basename(_file)
            self.lines.extend(["! Load Chain {}".format(_id),
                "open read unit 10 form name {}".format(base_name),
                "read sequ pdb offi unit 10",
                "close unit 10",
                "gene {} setup".format(segment),
                "open read unit 10 form name {}".format(base_name),
                "read coor pdb unit 10",
                "close unit 10"])
        self.lines.extend(ic_fill())

    def calculate_energy(self):
        self.energy_file = os.path.join(self.directory, "energy.dat")
        self.lines.extend(["! Calculate the energy.",
            "ener",
            "set tot ?ener",
            "open write card unit 10 name {}".format(os.path.basename(self.energy_file)),
            "write title unit 10",
            "*@tot",
            "close unit 10"])

    def relax(self):
        self.lines.extend(["! Carry out an energy minimization",
            "nbon nbxm 5",
            energy_line(self.experiment),
            "mini abnr nstep {} nprint 50 -".format(self.experiment.charmm_iterations),
            "tolgrd 0.01 tolenr 0.0001 tolstp 0.00"])

    def output_molecules(self):
        self.output_molecule_files = list()
        for _id in self.chain_ids:
            base_name = make_pdb_name(_id, "out")
            self.output_molecule_files.append(os.path.join(self.directory, base_name))
            segment = make_segment_name(_id)
            self.lines.extend(["open write unit 10 name {} card".format(base_name),
                "write coor sele segi {} end pdb unit 10 card".format(segment),
                "close unit 10"])

    def end_script(self):
        self.lines.append("STOP")

    def charmm(self):
        self.script_file = os.path.join(self.directory, "charmm_input.txt")
        with open(self.script_file, "w") as f:
            f.write("\n".join(self.lines))
        self.output_file = os.path.join(self.directory, "charmm_output.txt")
        with open(self.script_file) as fi, open(self.output_file, "w") as fo:
            proc = subprocess.Popen(standards.CharmmCommand, stdin=fi, stdout=fo)
            proc.wait()
        self.has_run = True

    def collect_garbage(self):
        self.experiment.safe_rmtree(self.directory)

    def __enter__(self):
        self.directory = tempfile.mkdtemp(prefix="charmm_", dir=self.experiment.get_temp())
        self.previous_directory = os.getcwd()
        os.chdir(self.directory)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.previous_directory)
        self.collect_garbage()


class Energy(CharmmProc):
    def __enter__(self):
        CharmmProc.__enter__(self)
        self.begin_script()
        self.load_input_files()
        self.load_molecules()
        self.calculate_energy()
        self.end_script()
        self.charmm()
        with open(self.energy_file) as f:
            try:
                self.energy = float(f.read())
            except ValueError:
                raise IOError("There was a problem with the CHARMM energy file {}".format(self.energy_file))
        return self


class InteractionEnergy(CharmmProc):
    def __init__(self, experiment, molecule_list):
        CharmmProc.__init__(self, experiment, molecule_list)
        if len(self.molecules) == 2:
            self.mol1, self.mol2 = self.molecules
        else:
            raise ValueError("An interaction energy calculation needs exactly two molecules.")

    def __enter__(self):
        CharmmProc.__enter__(self)
        with Energy(self.experiment, [self.mol1]) as e1:
            energy1 = e1.energy
        with Energy(self.experiment, [self.mol2]) as e2:
            energy2 = e2.energy
        with Energy(self.experiment, self.molecules) as e12:
            energy12 = e12.energy
        self.energy = energy12 - (energy1 + energy2)
        return self


class Relaxation(CharmmProc):
    def __init__(self, experiment, molecules, relaxed_file):
        CharmmProc.__init__(self, experiment, molecules)
        self.relaxed_file = relaxed_file

    def __enter__(self):
        CharmmProc.__enter__(self)
        self.begin_script()
        self.load_input_files()
        self.load_molecules()
        self.relax()
        self.output_molecules()
        self.end_script()
        self.charmm()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        #output_molecules = [molecules.Molecule(_id, _file, self.experiment) for _id, _file in zip(self.chain_ids, self.output_molecule_files)]
        relaxed_name = "relaxed"
        #output_molecule = molecules.merge(output_molecules, True, merged_name, self.relaxed_file, write_pdb=True)
        relaxed_molecule = molecules.Molecule(relaxed_name, self.relaxed_file, self.experiment)
        relaxed_molecule.assemble_from_CHARMM(self.chain_ids, self.output_molecule_files, self.resnums)
        CharmmProc.__exit__(self, exc_type, exc_value, traceback)


def introduction(purpose = None):
    """Create a header to start a CHARMM script."""
    return ["wrnl -2",
        "!prnl -2",
        "bomb -2"]
        

def ic_fill():
    """Tell CHARMM to include missing Atoms."""
    return ["! Add missing Atoms and assign them coordinates",
        "ic fill preserve",
        "ic param",
        "ic build",
        "hbuild"]


def energy_line(experiment):
    """Generate a string saying how energy should be calculated in CHARMM."""
    terms = experiment.charmm_energy_terms
    line = "skip all excl {}".format(" ".join(terms))
    return line


def make_segment_name(_id):
    return "ml{}".format(_id).lower()


def make_pdb_name(_id, io):
    if io not in ["in", "out"]:
        raise ValueError("The IO parameter must be 'in' or 'out,' not '{}.'".format(io))
    return "{}chain{}.pdb".format(io, _id).lower()
