# This VMD script is the NEWEST one (2.1.0) that  calculates the interaction energy between the antigen and the MAPs part for each antigen position.
# The following arguments are required:
# -e: the name of this script
# -args: 0: the antigen+part PSF
#        1: the antigen+part PDB
#        2: the file of positions to use
#        3: the name of the output file
#        4: the first parameter file
#        5: (optional) the second parameter file
#        etc.

# Load the VMD functions.
source vmd_functions.tcl
package require namdenergy

# Load the arguments.
if {[llength $argv] < 5} {
	puts "Usage: structure.psf coordinates.pdb positions/file.dat output/energies.dat experiment/details.txt parameter/file1.prm parameter/file2.prm (optional) ..."
	exit 1
}
set inStruct [lindex $argv 0]
set inCoords [lindex $argv 1]
set positionsFile [lindex $argv 2]
set outEnergy [lindex $argv 3]
set parameterFiles []
for {set i 4} {$i < [llength $argv]} {incr i} {
	lappend parameterFiles [lindex $argv $i]
}

# Load the antigen and the MAPs part.
mol new $inStruct
mol addfile $inCoords

set agSeg [antigenSegment]
set mapsSeg [mapsSegment]

# Select every atom in each segment. These segments were named by merge_antigen_part.tcl.
set Ag [atomselect 0 "segname $agSeg"]
set MAPsPart [atomselect 0 "segname $mapsSeg"]

# Define the initial position of the antigen, which is zero degrees and point (0, 0, 0).
set faceAngle 0
set pos "0 0 0"

# Make a file of interaction energies.
set eFile [open $outEnergy "w"]
close $eFile

# Calculate the energy at all positiions.
set positions [open $positionsFile]
while {[gets $positions P] >= 0} {
	# Get the difference between the needed and current coordinates.
	set dP [vecsub $P "$faceAngle $pos"]
	set dzRot [lindex $dP 0]
	set dx [lindex $dP 1]
	set dy [lindex $dP 2]
	set dz [lindex $dP 3]
	set dxyz "$dx $dy $dz"
	# Rotate the antigen to the specified direction.
	if {$dzRot != 0} {
		rotateZAxis $Ag $pos $dzRot
		set faceAngle [expr $faceAngle + $dzRot]
	}
	# Translate the antigen by the amount needed.
	$Ag moveby $dxyz
	set pos [vecadd $pos $dxyz]
	# FIXME: write files for debugging purposes only
	#$Ag writepdb "zr$faceAngle.pos$pos.pdb"
	puts " "
	puts $pos
	#puts [getAntigenPosition $Ag $epitope 0 $antigenSegment]
	# Calculate the interaction energy (electrostatic and van der Waals) between the antigen and the MAPs part.
	set energyList [lindex [namdenergy -sel $Ag $MAPsPart -elec -vdw -par $parameterFiles] 0]
        puts "Result: $energyList"
	# Read the energy from the fourth position in the list.
	set energy [lindex [split $energyList] 4]
        puts "Energy: $energy"
	# Write the energy to the file of all energies.
	set eFile [open $outEnergy "a"]
	puts $eFile "$P $energy"
	close $eFile
}
close $positions

exit 0
