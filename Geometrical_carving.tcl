set sel [atomselect top "all and not (x>-21.55 and x<21.55 and y<18.9 and y>-18.9)"]
set atomList [lsort -unique [$sel get {segname resid}]]
$sel delete
package require psfgen
resetpsf
readpsf CNC_11_11_19.psf
coordpdb CNC_11_11_19.pdb
foreach atom $atomList { delatom [lindex $atom 0] [lindex $atom 1] }
writepsf CNC_3_3_20.psf
writepdb CNC_3_3_20.pdb
mol new CNC_3_3_20.psf
mol addfile CNC_3_3_20.pdb

#source:https://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/23399.html
