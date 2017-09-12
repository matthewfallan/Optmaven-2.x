# This VMD script renumbers residues sequentially.
# The following arguments are required:
# -e: the name of this script
# -m: the files of the molecules
# -args: 0: the first atom number
#        1: the first residue number
#        2: the name of the output pdb

# Load the arguments.
if {[llength $argv] < 3} {
	puts "Usage: first_atom_number first_residue_number output/file.pdb"
	exit 1
}
set atom1 [lindex $argv 0]
set res1 [lindex $argv 1]
set outFile [lindex $argv 2]

set sel [atomselect $i all]
# Renumber atoms.
$sel set index [vecadd atom1 [$sel get index]]
# Renumber residues.
$sel set resid [vecadd res1 [$sel get residue]]

$sel writepdb $outFile
exit 0
