# This VMD script tests all of the combinations of levels in the OptMAVEn grid and outputs the combinations that do not cause clashes between the antigen and both prototype antibodies.
# The following arguments are required:
# -e: the name of this script.
# -m: 0: antigen/coords.pdb 
#     1: MoleculeH/coords.pdb
#     2: MoleculeK/coords.pdb
# -args: 0: grid/file.dat
#        1: epitope/file.txt
#        2: antigenFirstChain
#        3: output/file.dat
#        4: clashCutoff

# Load the VMD functions.
source vmd_functions.tcl

# Select every atom in each molecule.
set Ag [atomselect 0 "all"]
set IgH [atomselect 1 "all"]
set IgK [atomselect 2 "all"]

# Make sure each molecule exists.
if {[llength [$Ag get name]] == 0 || [llength [$IgH get name]] == 0 || [llength [$IgK get name]] == 0} {
    puts "One or more molecules failed to load."
    exit 1
}

# Load the arguments.
if {[llength $argv] != 5} {
    puts "Usage: -args grid/file.dat epitope/file.txt antigenFirstChain output/file.dat clashCutoff"
    exit 1
}

set gridFile [lindex $argv 0]
set epitopeFile [lindex $argv 1]
set antigenFirstChain [lindex $argv 2]
set positionsFile [lindex $argv 3]
set clashCutoff [lindex $argv 4]
set clashesPermitted 0

# Read the grid levels.
set xLevels [split [readAttribute $gridFile "x" 1]]
set yLevels [split [readAttribute $gridFile "y" 1]]
set zLevels [split [readAttribute $gridFile "z" 1]]
set zAngleLevels [split [readAttribute $gridFile "zAngle" 1]]

# Select the epitope.
set epitope [atomselect 0 [readFile $epitopeFile]]

# Determine the maximum z coordinate of any atom in either antibody.
set maxZH [lindex [measure minmax $IgH] {1 2}]
set maxZK [lindex [measure minmax $IgK] {1 2}]
set maxZIg [expr max($maxZH, $maxZK)]


# Test all positions and record those that do not clash.
set positions [open $positionsFile "w"]
foreach zAngle $zAngleLevels {
    foreach z $zLevels {
        # Determine the minimum z coordinate of any atom in the antigen.
        set minZAg [lindex [measure minmax $Ag] {0 2}]
        # Select all of the atoms in the antigen with z coordinates within the cutoff z coordinate of either antibody.
        set AgNear [atomselect 0 "z < $maxZIg + $clashCutoff"]
        # Select the atoms in the antibody with z coordinates within the cutoff z coordinates of the antigen.
        set IgHNear [atomselect 1 "z > $minZAg - $clashCutoff"]
        set IgKNear [atomselect 2 "z > $minZAg - $clashCutoff"]
        foreach y $yLevels {
            foreach x $xLevels {
                # Move the antigen to the correct position.
                repositionAntigen $Ag $epitope $antigenFirstChain $zAngle $x $y $z
                # Count the number of clashes between the antigen and the antibody.
                set clashes [llength [lindex [measure contacts $clashCutoff $AgNear $IgHNear] 0]]
                if {$clashes <= $clashesPermitted} {
                    set clashes [llength [lindex [measure contacts $clashCutoff $AgNear $IgKNear] 0]]
                    if {$clashes <= $clashesPermitted} {
                        # Write non-clashing positions to the file.
                        puts $positions "$zAngle $x $y $z"
                    }
                }
            }
        }
    }
}
close $positions

exit 0
