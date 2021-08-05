#!/usr/bin/env bash
# the next line restart using vmd \
exec vmd "-dispdev text -e  $0"


set NRES 6 

################################################################################
#                                                                              #
# Generate the index file for calculating the Phi/Psi angles of the residue    #
# by using Gromacs g_dihedral programs.                                        #
#                                                                              #
# This script is designed for cyclic peptide, since the g_rama fails to work   #
# for the calculation of N-terminal Phi angle and C-terminal Psi angle.        #
#                                                                              #
################################################################################


proc getNextResID { nres resi } {
  if { $resi == $nres } {
    return 1
  } else {
    return [expr $resi + 1 ] 
  }
}

proc getPreviousResID { nres resi } {
  if { $resi == 1 } {
    return $nres 
  } else {
    return [expr $resi - 1 ] 
  }
}


proc getAtomID {molid resid name} {
  set sel [atomselect $molid "resid $resid and name $name"]
  set atomid [$sel get serial]
#  set atomid [$sel get index]
  return $atomid
}


proc getPhiAtoms {molid resi} {
  global NRES 
  set resi0 [getPreviousResID $NRES $resi ]
  set atomidList [list]
  lappend atomidList [getAtomID $molid $resi0 "C"]
  lappend atomidList [getAtomID $molid $resi "N"]
  lappend atomidList [getAtomID $molid $resi "CA"]
  lappend atomidList [getAtomID $molid $resi "C"]
  return $atomidList
}


proc getPsiAtoms {molid resi} {
  global NRES 
  set resi1 [getNextResID $NRES $resi ]
  set atomidList [list]
  lappend atomidList [getAtomID $molid $resi "N"]
  lappend atomidList [getAtomID $molid $resi "CA"]
  lappend atomidList [getAtomID $molid $resi "C"]
  lappend atomidList [getAtomID $molid $resi1 "N"]
  return $atomidList
}



proc writeAtomList {ofp atomidList} {
  lassign $atomidList a0 a1 a2 a3 a4 
  puts $ofp "$a0 $a1 $a2 $a3"
  puts $ofp "$a1 $a2 $a3 $a4"
}


#############################################
set GRO "emin.gro"
set molid [ mol new $GRO type gro ]
#set PDB "prot.pdb"
#set molid [mol new $PDB type pdb]

# Set up residue list based on CA atoms
set calpha [atomselect $molid "name CA"]
set resids [$calpha get resid]


set ofp [open "BE_2D_MetaD_Input.dat" w]
puts $ofp {RANDOM_EXCHANGES}
puts $ofp ""

set icv 0
set irep 0


foreach resi $resids {
  set phi  [getPhiAtoms $molid $resi ]
  set psi  [getPsiAtoms $molid $resi ]

  puts $ofp "#Rep$irep; Res$resi:phi/psi"
  puts $ofp "cv$icv: TORSION ATOMS=[join $phi ","]" 
  set icv [expr $icv+1]
  puts $ofp "cv$icv: TORSION ATOMS=[join $psi ","]" 
  set icv [expr $icv+1]
  set irep [expr $irep+1]
  puts $ofp " "
}


foreach resi $resids {
  set resi1 [getNextResID $NRES $resi]
  set psi  [getPsiAtoms $molid $resi ]
  set phi1 [getPhiAtoms $molid $resi1 ]

  puts $ofp "#Rep$irep; Res$resi:psi/Res$resi1:phi"
  puts $ofp "cv$icv: TORSION ATOMS=[join $psi ","]" 
  set icv [expr $icv+1]
  puts $ofp "cv$icv: TORSION ATOMS=[join $phi1 ","]" 
  set icv [expr $icv+1]
  set irep [expr $irep+1]
  puts $ofp " "
}

puts $ofp "METAD ARG=cvX,cvY SIGMA=0.31416,0.31416 HEIGHT=0.1 PACE=2000 LABEL=metad.cvX FILE=HILLS
PRINT ARG=cvX,cvY STRIDE=500 FILE=COLVAR

ENDPLUMED"


close $ofp
exit 0
