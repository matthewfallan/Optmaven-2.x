""" This module stores information that is shared between the modules. """

import os

import numpy as np

OptmavenVersion = "2.1.0"
OptmavenName = "Optmaven-{}".format(OptmavenVersion)

SourceDirectory = os.path.dirname(os.path.realpath(__file__))
OptmavenDirectory = os.path.dirname(SourceDirectory)
ExperimentsDirectory = os.path.join(OptmavenDirectory, "experiments")
DataDirectory = os.path.join(OptmavenDirectory, "data")
PDBDirectory = os.path.join(DataDirectory, "pdbs")
InputsDirectory = os.path.join(DataDirectory, "input_files")
AntibodiesDirectory = os.path.join(DataDirectory, "antibodies")
ScaffoldsDirectory = os.path.join(AntibodiesDirectory, "scaffolds")
MAPsDirectory = os.path.join(AntibodiesDirectory, "maps")

# Check to make sure all of the directories and files are present. If not, the installation is malformed.
if os.path.basename(OptmavenDirectory) != OptmavenName:
    raise OSError("This installation of Optmaven is located in {} but needs to be located in {}.".format(OptmavenDirectory, OptmavenName))

contents = {
    "modules": ["EXPERIMENT.py", "STANDARDS.py"],
    "pdbs": None,
    "programs": ["Optmaven.py"]
}

# FIXME: implement file checking


if not os.path.isdir(ExperimentsDirectory):
    os.mkdir(ExperimentsDirectory)

# Maximum line length in console.
ConsoleWidth = 80
# Maximum number of list items to print.
MaxListDisplay = 20
SelectAll = "all"
SelectNone = "none"
EscapeCharacter = "\\"

# Optmaven grid settings.
DefaultOptmavenGrid_x = np.linspace(-5, 10, 7)
DefaultOptmavenGrid_y = np.linspace(-10, 5, 7)
DefaultOptmavenGrid_z = np.linspace(3.75, 16.25, 11)
DefaultOptmavenGrid_zAngle = np.linspace(0, 300, 6)

# Input files.
DefaultTopologyFile = os.path.join(InputsDirectory, "top_all27_prot_na.rtf")
DefaultParameterFile = os.path.join(InputsDirectory, "par_all27_prot_na.prm")
DefaultSolvationFile = os.path.join(InputsDirectory, "solvation.dat")

# CHARMM configuration.
DefaultCharmmEnergyTerms = ['angl', 'bond', 'dihe', 'elec', 'impr', 'urey', 'vdw']
DefaultCharmmIterations = 5000
CharmmCommand = "/gpfs/group/cdm/c34b1.xj.gnu/exec/gnu/charmm.serial.xlarge" 

# VMD configuration.
VmdCommand = "vmd"
VmdDisp = "-dispdev"
VmdNoDisp = "none"
VmdExec = "-e"
VmdMolecules = "-m"
VmdFrames = "-f"
VmdArgs = "-args"

# Optmaven files.
ScaffoldChains = ["H", "K"]
ScaffoldAntibodies = {chain: os.path.join(ScaffoldsDirectory, "Molecule{}.pdb".format(chain)) for chain in ScaffoldChains}

# Antigen positioning parameters.
DefaultClashCutoff = 1.25  # Angstroms

# Rotation matrix routines.
def rotate_vi_to_vf(vi, vf):
    """ Compute a rotation matrix to rotate the 3D coordinate vi to vf. """
    # Normalize both vectors.
    vi, vf = np.array(vi), np.array(vf)
    vi /= np.linalg.norm(vi)
    vf /= np.linalg.norm(vf)
    # Find the axis of rotation, which is perpendicular to both vectors.
    ax = np.cross(vf, vi)
    # Calculate the sine and cosine of the angle of rotation.	
    sin = np.linalg.norm(ax)
    cos = np.dot(vf, vi)
    # Normalize the axis.
    ax = ax / sin
    # Compute the cross product matrix of the axis.
    cpm = np.array([
    	[     0, -ax[2],  ax[1]],
    	[ ax[2],      0, -ax[0]],
    	[-ax[1],  ax[0],      0]
    ])
    # Compute the tensor product matrix of the axis.
    tpm = ax.reshape(1, 3) * ax.reshape(3, 1)
    # Return the rotation matrix.
    return cos * np.eye(3) + sin * cpm + (1 - cos) * tpm
