""" This module stores information that is shared between the modules. """

import itertools
import os
import random
import string

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
MapsDirectory = os.path.join(AntibodiesDirectory, "maps")

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

# Experimental configuration.
DefaultNumberOfDesigns = 50

# Default PBS settings.
DefaultWalltime = 86399  # 23 hours, 59 minutes, 59 seconds
DefaultBatchSize = 10
PbsQueue = "lionxv"
PbsQsub = "qsub"
PbsArrayId = "$PBS_ARRAYID"
PbsJobFilePrefix = "job-"

# Python settings.
PythonCommand = "python"
PythonModule = "python/2.7.5.ucs4"

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
VmdFunctions = os.path.join(SourceDirectory, "vmd_functions.tcl")

# CPLEX configuration.
CplexDirectory = "/usr/global/ilog/CPLEX_Studio124/cplex/python/x86-64_sles10_4.1"

# Optmaven files.
ScaffoldChains = ["H", "K"]
ScaffoldAntibodies = {chain: os.path.join(ScaffoldsDirectory, "Molecule{}.pdb".format(chain)) for chain in ScaffoldChains}

# Antigen positioning parameters.
DefaultClashCutoff = 1.25  # Angstroms

# MAPs database standards.
MapsHeavyChains = ["H"]
MapsLightChains = ["K", "L"]
MapsNamesakeHeavy = "H"
MapsNamesakeLight = "L"
MapsChains = MapsHeavyChains + MapsLightChains
MapsNamesakeChains = [MapsNamesakeHeavy, MapsNamesakeLight]
MapsRegions = ["V", "CDR3", "J"]
MapsCdrs, MapsNamesakeCdrs = [["{}{}".format(chain, cdr) for chain, cdr in itertools.product(chains, MapsRegions)] for chains in [MapsChains, MapsNamesakeChains]]
MapsIntegerCutsFile = os.path.join(MapsDirectory, "MAPs_Integer_Cuts.txt")

# Antibody coordinate standards.
MapsCoordSep = "."
MapsCoordDimension = 3
def make_maps_coord(cdr, dimension):
    return "{}{}{}".format(cdr, MapsCoordSep, dimension)
def to_maps_coord(item):
    if MapsCoordSep in item:
        split = item.split(MapsCoordSep)
        if len(split) == 2:
            cdr, dimension = split
            if cdr in MapsCdrs and str(dimension) in map(str, range(MapsCoordDimension)):
                return cdr, int(dimension)
AngleLabel = "Angle"
xLabel = "x"
yLabel = "y"
zLabel = "z"
zAngleLabel = zLabel + AngleLabel
SinLabel = "Sin"
CosLabel = "Cos"
PositionOrder = [zAngleLabel, xLabel, yLabel, zLabel]
AtomCoordOrder = [xLabel, yLabel, zLabel]
PositionCoordOrder = list(itertools.chain(*[["{}{}".format(item, label) for label in [SinLabel, CosLabel]] if AngleLabel in item else [item] for item in PositionOrder]))
MapsCoordOrder = [make_maps_coord(cdr, dim) for cdr, dim in itertools.product(MapsNamesakeCdrs, range(MapsCoordDimension))]
EnergyLabel = "energy"
CoordOrder = PositionCoordOrder + MapsCoordOrder
DefaultGapPenalty = 8

# Kmeans standards.
DefaultKmeansMaxIterations = 1000
DefaultKmeansTolerance = 0.001
DefaultKmeansOptKThreshold = 0.5

# Rotation matrix routines.

_ROTATION_MATRIX_DIMENSION = len(AtomCoordOrder)
dim = _ROTATION_MATRIX_DIMENSION
xAxis, yAxis, zAxis = np.eye(dim)

def rotate_vi_to_vf(vi, vf):
    """ Compute a rotation matrix to rotate the 3D coordinate vi to vf. """
    if len(vi) != dim or len(vf) != dim:
        raise ValueError("The rotation matrix function requires 3D coordinates.")
    # Normalize both vectors.
    vi, vf = np.array(vi), np.array(vf)
    vin = np.linalg.norm(vi)
    vfn = np.linalg.norm(vf)
    if np.isclose(vin, 0) or np.isclose(vfn, 0):
        # Abort rotation if either vector is almost zero.
        return np.eye(dim)
    vi /= vin
    vf /= vfn
    # Find the axis of rotation, which is perpendicular to both vectors.
    ax = np.cross(vf, vi)
    # Calculate the sine and cosine of the angle of rotation.	
    sin = np.linalg.norm(ax)
    cos = np.dot(vf, vi)
    return rotate_axis_angle(ax, sin=sin, cos=cos)


def rotate_axis_angle(axis, angle=None, sin=None, cos=None):
    """ Create a rotation matrix to rotate by an angle around an axis. """
    if len(axis) != dim:
        raise ValueError("The rotation matrix function requires 3D coordinates.")
    if sin is None:
        sin = np.sin(angle)
    if cos is None:
        cos = np.cos(angle)
    axn = np.linalg.norm(axis)
    if np.isclose(axn, 0):
        # If the rotation axis is almost zero, then the rotation angle is either almost 0 or 180.
        # If the angle is 0, then cos = 1. Return identity matrix.
        # If the angle is 180, then cos = -1. Return the negative identity matrix.
        return cos * np.eye(dim)
    # Normalize the axis.
    ax = axis / axn
    # Compute the cross product matrix of the axis.
    cpm = np.array([
    	[     0, -ax[2],  ax[1]],
    	[ ax[2],      0, -ax[0]],
    	[-ax[1],  ax[0],      0]
    ])
    # Compute the tensor product matrix of the axis.
    tpm = ax.reshape(1, dim) * ax.reshape(dim, 1)
    # Return the rotation matrix.
    return cos * np.eye(dim) + sin * cpm + (1 - cos) * tpm


def random_string(length, alphabet=None):
    if alphabet is None:
        alphabet = string.letters
    return "".join([random.choice(alphabet) for i in range(length)])


def is_path_component(file_name):
    """ Determine if a file name is a single file without whitespace. """
    fn = os.path.basename(file_name).split()
    return len(fn) > 0 and fn[0] == file_name
