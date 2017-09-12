# This VMD script renumbers residues sequentially.
# The following arguments are required:
# -e: the name of this script
# -m: the files of the molecules
# -args: 0: the number of molecules
#        1: the first atom number
#        2: the first residue number


# Load the arguments.
if {[llength $argv] < 3} {
	puts "Usage: number_of_molecules first_atom_number first_residue_number"
	exit 1
}
set nMols [lindex $argv 0]
set atom1 [lindex $argv 1]
set res1 [lindex $argv 2]


for {set i 0} {$i < $nMols} {incr $i} {
    set sel [atomselect $i all]
    # Renumber atoms.
    $sel set index [vecadd atom1 [$sel get index]]
    set atomNums [$sel get index]
    set atom1 [expr {max($atomNums)}]
    # Renumber residues.
    $sel set resid [vecadd res1 [$sel get residue]]
    set resNums [$sel get resid]
    set res1 [expr {max($resNums)}]
}


