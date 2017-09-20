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
        self.garbage = list()

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
                self.garbage.append(symlink)
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
            self.garbage.append(_file)
            self.lines.extend(["! Load Chain {}".format(_id),
                "open read unit 10 form name {}".format(base_name),
                "read sequ pdb offi unit 10",
                "close unit 10",
                "gene {} setup".format(segment),
                "open read unit 10 form name {}".format(base_name),
                "read coor pdb unit 10",
                "close unit 10"])
        self.lines.extend(ic_fill())

    def energy(self):
        handle, self.energy_file = tempfile.mkstemp(dir=self.directory, prefix="energy_", suffix=".dat")
        os.close(handle)
        self.garbage.append(self.energy_file)
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
            #"tolgrd 0.01 tolenr 0.0001 tolstp 0.00"])
            "tolgrd 0.1 tolenr 0.1 tolstp 0.00"])
    #FIXME: change tolgrd and tolenr to 0.01 and 0.0001

    def output_molecules(self):
        self.output_molecule_files = list()
        for _id in self.chain_ids:
            base_name = make_pdb_name(_id, "out")
            self.output_molecule_files.append(os.path.join(self.directory, base_name))
            segment = make_segment_name(_id)
            self.lines.extend(["open write unit 10 name {} card".format(base_name),
                "write coor sele segi {} end pdb unit 10 card".format(segment),
                "close unit 10"])
        self.garbage.extend(self.output_molecule_files)

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
        self.garbage.extend([self.script_file, self.output_file])
        self.has_run = True

    def collect_garbage(self):
        for _file in self.garbage:
            try:
                os.remove(_file)
            except OSError:
                pass
        try:
            os.rmdir(self.directory)
        except OSError:
            pass

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
        self.energy()
        self.end_script()
        self.charmm()
        with open(self.energy_file) as f:
            self.energy = float(f.read())
        return self


class InteractionEnergy(CharmmProc):
    def __init__(self, experiment, molecule_list):
        CharmmProc.__init__(self, experiment, molecule_list)
        if len(self.molecules) == 2:
            self.mol1, self.mol2 = self.molecules
        else:
            raise ValueError("An interaction energy calculation needs exactly two molecules.")

    def __enter__(self):
        with Energy(self.experiment, [self.mol1]) as e1:
            energy1 = e1.energy
        with Energy(self.experiment, [self.mol2]) as e2:
            energy2 = e2.energy
        with Energy(self.experiment, self.molecules) as e12:
            energy12 = e12.energy
        self.energy = e12 - (e1 + e2)
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


# Create the actual functions that are run from outside of this file
def Missing_Atoms(molecules, experiment = None):
    """Add missing Atoms to Molecules."""
    # Validate the Molecules
    molecules, gn = validate_molecules(molecules)
    # Declare the procedure, making sure it is OK (it is, but whatever)
    procedure = validate_procedure("add_missing_atoms")
    # Determine who is running this experiment
    try:
        user = experiment["User"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        user = defaultUser
    # Create the CHARMM script
    script = introduction("Add missing Atoms to Molecules.")
    script += load_input_files(experiment)
    script += load_molecules(molecules, procedure, user, "all")
    script += ic_fill() 
    script += output_molecules(molecules, procedure, "all") + "stop\n"
    # Run CHARMM
    input, output = execute_CHARMM_script(script, procedure, gn)
    # Load the new structures of the Molecules
    load_structures(molecules, procedure, "all")
    # Clean up after the procedure
    clean_up(molecules, input, output, procedure)

'''
def Relaxation(molecules, experiment = None):
    """Use CHARMM to run an energy minimization."""
    # Validate the molecules
    molecules, gn = validate_molecules(molecules)
    # Declare the procedure
    procedure = validate_procedure("relaxation")
    # Determine who is running the experiment
    try:
        user = experiment["User"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        user = defaultUser
    # Create the CHARMM script
    if gn:
        script = introduction("The Relaxation of Design Group " + str(gn))
    else:
        script = introduction("A Structure Relaxation")
    script += load_input_files(experiment)
    script += load_molecules(molecules, procedure, user, "all")
    script += ic_fill() + minimize(molecules, experiment, procedure, gn)
    script += output_molecules(molecules, procedure, "all") + "stop\n"
    # Run CHARMM
    input, output = execute_CHARMM_script(script, procedure, gn)
    # Load the new structures of the Molecules
    load_structures(molecules, procedure, "all")
    # Clean everything up
    clean_up(molecules, input, output, procedure)
'''

def Perturbation(molecules, angles, experiment = None):
    """Use CHARMM to do a backbone perturbation."""
    # Validate the molecules
    molecules, gn = validate_molecules(molecules)
    # Declare the procedure
    procedure = validate_procedure("perturbation")
    # Determine who is running the experiment
    try:
        user = experiment["User"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        user = defaultUser
    # Create the CHARMM script
    if isinstance(gn, int):
        script = introduction("The Perturbation of Design Group " + str(gn))
    else:
        script = introduction("A Backbone Perturbation")
    script += load_input_files(experiment)
    script += load_molecules(molecules, procedure, user, "all")
    script += ic_fill()
    # Specify the names of the Atoms to use in the restraints based on the file
    # format
    if molecules[0].fileFormat == "PDB":
        A1 = ' C '
        A2 = ' N '
        A3 = ' CA '
        A4 = ' C '
        A5 = ' N '
    # If the file format isn't supported
    else:
        text = "The CHARMM Perturbation function does not support the "
        text += str(molecules[0].fileFormat) + " file format."
        raise CHARMM_Error(text)
    # Loop through the Molecules
    for molecule in molecules:
        # Calculate the Molecule's dihedral angles
        molecule.calculate_dihedrals()
        # Loop through the Residues by index
        for i in range(len(molecule)):
            # Get the Residue, the Residue's Molecule's name, the Residue's
            # name, and the Residue's number
            res2 = molecule[i]
            mn = res2.moleculeName
            rn = res2.name
            n = res2.number
            # If this Residue has perturbed angles
            if mn in angles and rn in angles[mn]:
                # Do the phi angle - iff there is a phi value listed for this
                # Residue, it is not the first Residue in the Molecule, the
                # listed value is a floating point number, and the Residue has a
                # phi dihedral angle
                if "phi" in angles[mn][rn] and i != 0 and res2.phi != None and \
                isinstance(angles[mn][rn]["phi"], float):
                    # Try to access the appropriate Atoms to make sure they 
                    # exist and that at least one of them can move
                    error = False
                    try:
                        # Get the previous Residue
                        res1 = molecule[i-1]
                        # Make a list of the Residues / Atoms to validate the
                        # existance of and make sure they're not all fixed in
                        # place
                        sets = [[res1, A1], [res2, A2], [res2, A3], [res2, A4]]
                        count = 0
                        for set in sets:
                            atom = set[0][set[1].strip()]
                            if set[0].freedom == "FIXED" or atom.name in \
                            set[0].fixedAtoms:
                                count += 1
                        # If all 4 Atoms are fixed in place, skip this restraint
                        if count == 4:
                            error = True
                    # If there was an error accessing any of the Atoms, skip the
                    # restraint
                    except MOLECULES.MoleculeError:
                        error = True
                    # Get the new phi angle
                    if not error:
                        phi = res2.phi + angles[mn][rn]["phi"]
                        if phi > 180.0:
                            phi -= 360.0
                        if phi < -180.0:
                            phi += 360.0
                        phi = format(phi, '.3f') + "\n"
                        # Create the restraint
                        script += "cons dihe " + str(n-1) + A1+str(n)+A2+str(n)
                        script += A3 + str(n) + A4 + "force 32800.0 min "+phi
                # Do the same for psi
                if "psi" in angles[mn][rn] and i != len(molecules) - 1 and \
                res2.psi != None and isinstance(angles[mn][rn]["psi"], float):
                    error = False
                    try:
                        res3 = molecule[i+1]
                        sets = [[res2, A2], [res2, A3], [res2, A4], [res3, A5]]
                        count = 0
                        for set in sets:
                            atom = set[0][set[1].strip()]
                            if set[0].freedom == "FIXED" or atom.name in \
                            set[0].fixedAtoms:
                                count += 1
                        if count == 4:
                            error = True
                    except MOLECULES.MoleculeError:
                        error = True
                    if not error:
                        psi = res2.psi + angles[mn][rn]["psi"]
                        if psi > 180.0:
                            psi -= 360.0
                        if psi < -180.0:
                            psi += 360.0
                        script += "cons dihe " + str(n) + A2 + str(n) + A3
                        script += str(n) + A4 + str(n+1) + A5 + " force 32800.0"
                        script += " min " + format(psi, '.3f') + "\n"
    # Now finish the script by minimizing the structure and outputting the
    # molecules
    script += minimize(molecules, experiment, procedure, gn)
    script += output_molecules(molecules, procedure, "all") + "stop\n"
    # Run CHARMM
    input, output = execute_CHARMM_script(script, procedure, gn)
    # Load the Molecules' new structures
    load_structures(molecules, procedure, "all")
    # Clean up the files
    clean_up(molecules, input, output, procedure)

def Energy(molecules, experiment = None, which = "all"):
    """Use CHARMM to calculate the complex energy of a group of Molecules."""
    # Validate the molecules
    molecules, gn = validate_molecules(molecules)
    # Create the procedure
    procedure = validate_procedure("energy")
    # Determine who is doing the calculation
    try:
        user = experiment["User"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        user = defaultUser
    # Determine if solvation is or should be used
    solvation = SOLVATION.get_string(experiment, procedure)
    # Generate the CHARMM script to do the energy calculation
    if isinstance(gn, int):
        script = introduction("Energy Calculation for Design Group " + str(gn))
    else:
        script = introduction("A CHARMM energy calculation")
    script += load_input_files(experiment)
    script += load_molecules(molecules, procedure, user, which)
    script += ic_fill() + solvation
    script += energy_line(experiment, procedure, '', '', '', solvation)
    # Declare the name of the output file
    energyFile = "energy_values.txt"
    script += "\nener\n\nset tot ?ener\n\nopen write card unit 10 name "
    script += energyFile + "\n\nwrite title unit 10\n* @tot\n*\nclose unit 10\n"
    script += "stop\n"
    # Run CHARMM
    input, output = execute_CHARMM_script(script, procedure, gn)
    # Read in the calculated energy
    try:
        file = open(energyFile, "r")
        energy = float(file.readline())
        file.close()
    except IOError:
        text = "A CHARMM energy calculation failed to calculate an energy. "
        text += "Please review " + output + " to determine why."
        raise CHARMM_Error(text)
    # Clean up the folder
    #clean_up(molecules, input, output, procedure, energyFile)
    # Return the calculated energy
    return energy
