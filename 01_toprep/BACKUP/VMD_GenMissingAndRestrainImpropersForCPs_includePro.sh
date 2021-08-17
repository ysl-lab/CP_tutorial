#!/usr/bin/env bash
# the next line restart using vmd \
exec vmd "-dispdev text -e  $0"


set NRES 5
set IRES 1

################################################################################
#                                                                              #
# Generate the index file for calculating the Phi/Psi angles of the residue    #
# by using Gromacs g_dihedral programs.                                        #
#                                                                              #
# This script is designed for cyclic peptide, since the g_rama fails to work   #
# for the calculation of N-terminal Phi angle and C-terminal Psi angle.        #
#                                                                              #
################################################################################


proc getNextResID { resi } {
  global NRES
  global IRES
  if { $resi == $NRES } {
    return $IRES 
  } else {
    return [expr $resi + 1 ] 
  }
}

proc getPreviousResID {resi } {
  global NRES
  global IRES
  if { $resi == $IRES } {
    return $NRES
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

proc getResName {molid resid name} {
  set sel [atomselect $molid "resid $resid and name $name"]
  set residuename [$sel get resname]
#  set atomid [$sel get index]
  puts $residuename
  return $residuename
}

proc getOmega1Atoms {molid resi} {
  set resi0 [getPreviousResID $resi ]
  set resiname [getResName $molid $resi "CA"]
  set atomidList [list]
  if { $resiname == "PRO" } {
    lappend atomidList [getAtomID $molid $resi0 "C"]
    lappend atomidList [getAtomID $molid $resi  "CA"]
    lappend atomidList [getAtomID $molid $resi "N"]
    lappend atomidList [getAtomID $molid $resi "CD"]
  } else {
    lappend atomidList [getAtomID $molid $resi0 "C"]
    lappend atomidList [getAtomID $molid $resi  "CA"]
    lappend atomidList [getAtomID $molid $resi "N"]
    lappend atomidList [getAtomID $molid $resi "H"]
  }
  return $atomidList
}


proc getOmega2Atoms {molid resi} {
  global NRES 
  set resi1 [getNextResID $resi ]
  set atomidList [list]
  lappend atomidList [getAtomID $molid $resi "CA"]
  lappend atomidList [getAtomID $molid $resi1 "N"]
  lappend atomidList [getAtomID $molid $resi "C"]
  lappend atomidList [getAtomID $molid $resi "O"]
  return $atomidList
}

proc getRestrainImpAtoms {molid resi} {
  global NRES
  set resi1 [getNextResID $resi ]
  set resi1name [getResName $molid $resi1 "CA"]
  set atomidList [list]
  if { $resi1name == "PRO" } {
    lappend atomidList [getAtomID $molid $resi1 "CD"]
    lappend atomidList [getAtomID $molid $resi1 "N"]
    lappend atomidList [getAtomID $molid $resi "C"]
    lappend atomidList [getAtomID $molid $resi "O"]
  } else {
    lappend atomidList [getAtomID $molid $resi1 "H"]
    lappend atomidList [getAtomID $molid $resi1 "N"]
    lappend atomidList [getAtomID $molid $resi "C"]
    lappend atomidList [getAtomID $molid $resi "O"]
  }
  return $atomidList
}

#proc writeAtomList {ofp atomidList} {
#  lassign $atomidList a0 a1 a2 a3 a4
#  puts $ofp "$a0 $a1 $a2 $a3"
#  puts $ofp "$a1 $a2 $a3 $a4"
#}

#############################################
set GRO "02_prot.gro"
set molid [ mol new $GRO type gro ]

# Set up residue list based on CA atoms
set calpha [atomselect $molid "name CA"]
set resids [$calpha get resid]

set ofp [open "Impropers.dat" w]
set Omega2  [getOmega2Atoms $molid $NRES ]
lassign $Omega2 a1 a2 a3 a4
puts $ofp "$a1 $a2 $a3 $a4"
set Omega1  [getOmega1Atoms $molid $IRES ]
lassign $Omega1 a1 a2 a3 a4
puts $ofp "$a1 $a2 $a3 $a4"

foreach resi $resids {
  set RestrainImp  [getRestrainImpAtoms $molid $resi ]
  lassign $RestrainImp a1 a2 a3 a4
  puts $ofp "$a1 $a2 $a3 $a4"
}

close $ofp
exit 0
