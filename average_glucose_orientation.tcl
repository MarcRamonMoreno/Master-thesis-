#Orientation of CTA molecules
#By Jordi Faraudo 
#November 2017
# 
# Use with: vmd -dispdev text -e <nameofscript.tcl>
#
#
# PROCEDURE 1: BIG DCD
#First of all I include here the BigDCD command
#

proc bigdcd { script type args } {
    global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame bigdcd_running
  
    set bigdcd_running 1
    set bigdcd_frame 0
    set bigdcd_firstframe [molinfo top get numframes]
    set bigdcd_proc $script

    # backwards "compatibility". type flag is omitted.
    if {[file exists $type]} { 
        set args [linsert $args 0 $type] 
        set type auto
    }
  
    uplevel #0 trace variable vmd_frame w bigdcd_callback
    foreach dcd $args {
        if { $type == "auto" } {
            mol addfile $dcd waitfor 0
        } else {
            mol addfile $dcd type $type waitfor 0
        }
    }
    after idle bigdcd_wait
}

proc bigdcd_callback { tracedvar mol op } {
    global bigdcd_frame bigdcd_proc bigdcd_firstframe vmd_frame
    set msg {}
 
    # If we're out of frames, we're also done 
    # AK: (can this happen at all these days???). XXX
    set thisframe $vmd_frame($mol)
    if { $thisframe < $bigdcd_firstframe } {
        puts "end of frames"
        bigdcd_done
        return
    }
 
    incr bigdcd_frame
    if { [catch {uplevel #0 $bigdcd_proc $bigdcd_frame} msg] } { 
        puts stderr "bigdcd aborting at frame $bigdcd_frame\n$msg"
        bigdcd_done
        return
    }
    animate delete beg $thisframe end $thisframe $mol
    return $msg
}

proc bigdcd_done { } {
    global bigdcd_running
    
    if {$bigdcd_running > 0} then {
        uplevel #0 trace vdelete vmd_frame w bigdcd_callback
        puts "bigdcd_done"
        set bigdcd_running 0
    }
}

proc bigdcd_wait { } {
    global bigdcd_running bigdcd_frame
    while {$bigdcd_running > 0} {
        global bigdcd_oldframe
        set bigdcd_oldframe $bigdcd_frame
        # run global processing hooks (including loading of scheduled frames)
        display update ui
        # if we have read a new frame during then the two should be different.
        if { $bigdcd_oldframe == $bigdcd_frame } {bigdcd_done}
    }
}

#
# Define procedure for selecting molecules and compute angle
#

proc dothejob { frame } {

# definition of global variables
  global sel_N fptop fpmiddle fpbottom mol unitz M_PI
  package require pbctools
  
# show frame number  
#  puts "$frame"
    
# update selections
#  $seleccio frame $frame
  $sel_N frame $frame
#  $seleccio update
  $sel_N update

# count number of atoms and number of molecules
#  set na [$seleccio num] 
#  set n [$sel_N num]
#  puts "frame: $frame ; number molecules: $n" 
#  puts $fp "$frame $n" 

# compute center of mass of Nitrogen atoms and z component
  #set com1 [measure center $sel_N weight mass]
  #set comz [vecdot $unitz $com1]

# select N of top  bottom layers
  set topN [atomselect $mol "type OC3C61 and z>60" frame $frame]
  set middleN [atomselect $mol "type OC3C61 and (z>-20 and z<20)" frame $frame]
  set bottomN [atomselect $mol "type OC3C61 and z<-60" frame $frame]
# check number of top and bottom molecules
  set counttop [$topN num]
  set countbottom [$bottomN num]
  set countmiddle [$middleN num]
  puts "frame: $frame ; num_top: $counttop ; num_bottom $countbottom ; num_middle $countmiddle"
#puts "Top region"
#show progress
#  puts "frame: $frame ; center: $comz ; number molecules: $n"
  set avcosang 0
  set avang 0
  set timesec 0
  foreach restop [$topN get residue] {
		set modulus [expr {abs(fmod($restop, 2))}]
		
# Check the modulus and perform atom selection
		if {$modulus == 0} {
			set cap [atomselect $mol "residue $restop and name C5" frame $frame]
			set cua [atomselect $mol "residue $restop and name C2" frame $frame]
			#puts "Even_frame: $restop"
		} else {
			set cap [atomselect $mol "residue $restop and name C2" frame $frame]
			set cua [atomselect $mol "residue $restop and name C5" frame $frame]
			#puts "Odd_frame: $restop"
		}
   	#take coordinates
    	set capcoords [lindex [$cap get {x y}] 0]
    	set cuacoords [lindex [$cua get {x y}] 0]

        #puts "cap_coords: $capcoords"
        #puts "cua_coords: $cuacoords"

#      	puts "$nc(1) $nc(2)"
#        set nc(1) [lindex $capcoords 0]
#        set nc(2) [lindex $cuacoords 0]
#        puts "$nc(1) $nc(2)"
#        set vdist [vecsub $nc(1) $nc(2)] 
      	set vdist [vecsub $capcoords $cuacoords] 
      	#puts "vdist: $vdist"
      	set dist [veclength $vdist] 
      	#puts "dist: $dist"
      	set normz [veclength $unitz]
      	#puts "normz: $normz"
		set dz [vecdot $unitz $vdist]
      	#puts "dz: $dz"
      	set normprod [expr $dist*$normz]
      	#puts "normprod: $normprod"
      	set cosang [expr $dz/$normprod ]
      	#puts "cosinus: $cosang"
      	set ang [expr acos($cosang) ]
      	#puts "angle: $ang"
      	set angdeg [expr 180.0*$ang/$M_PI]
      	#puts "angle_deg: $angdeg"
#    if {$angdeg > 180} {
#       puts "reducing $angdeg ..." 
#       set angdeg [expr $angdeg - 180.0]
#       puts "now $angdeg" 
#    }
		set avcosang [expr $avcosang+$cosang]
		set avang [expr $avang+$angdeg]

		$cua delete 
		$cap delete     
  }
  	set timesec [expr $frame*2*10**-6]
  	set avcosang [expr $avcosang/$counttop]
	set avang [expr $avang/$counttop]
	puts $fptop "$frame $timesec $avcosang $avang"
  #puts "Middle region"
    set avcosang 0
	set avang 0
	set timesec 0
   foreach resmiddle [$middleN get residue] {
	   
		set modulus [expr {abs(fmod($resmiddle, 2))}]
		
# Check the modulus and perform atom selection
		if {$modulus == 0} {
			set cap [atomselect $mol "residue $resmiddle and name C5" frame $frame]
			set cua [atomselect $mol "residue $resmiddle and name C2" frame $frame]
			#puts "Even_frame: $resmiddle"
		} else {
			set cap [atomselect $mol "residue $resmiddle and name C2" frame $frame]
			set cua [atomselect $mol "residue $resmiddle and name C5" frame $frame]
			#puts "Odd_frame: $resmiddle"
		}
   	#take coordinates
    	set capcoords [lindex [$cap get {x y}] 0]
    	set cuacoords [lindex [$cua get {x y}] 0]

        #puts "cap_coords: $capcoords"
        #puts "cua_coords: $cuacoords"

#      	puts "$nc(1) $nc(2)"
#        set nc(1) [lindex $capcoords 0]
#        set nc(2) [lindex $cuacoords 0]
#        puts "$nc(1) $nc(2)"
#        set vdist [vecsub $nc(1) $nc(2)] 
      	set vdist [vecsub $capcoords $cuacoords] 
      	 #puts "vdist: $vdist"
      	set dist [veclength $vdist] 
      	 #puts "dist: $dist"
      	set normz [veclength $unitz]
      	 #puts "normz: $normz"
		set dz [vecdot $unitz $vdist]
      	 #puts "dz: $dz"
      	set normprod [expr $dist*$normz]
      	 #puts "normprod: $normprod"
      	set cosang [expr $dz/$normprod ]
      	 #puts "cosinus: $cosang"
      	set ang [expr acos($cosang) ]
      	set angdeg [expr 180.0*$ang/$M_PI]
#    if {$angdeg > 180} {
#       puts "reducing $angdeg ..." 
#       set angdeg [expr $angdeg - 180.0]
#       puts "now $angdeg" 
#    }
		set avcosang [expr $avcosang+$cosang]
		set avang [expr $avang+$angdeg]

	$cua delete 
	$cap delete     
  }
	set timesec [expr $frame*2*10**-6]
    set avcosang [expr $avcosang/$countmiddle]
	set avang [expr $avang/$countmiddle]
	puts $fpmiddle "$frame $timesec $avcosang $avang"
	
 #puts "Bottom region"
   set timesec 0
   set avcosang 0
   set avang 0
  foreach resbottom [$bottomN get residue] {
		set modulus [expr {abs(fmod($resbottom, 2))}]
		
# Check the modulus and perform atom selection
		if {$modulus == 0} {
			set cap [atomselect $mol "residue $resbottom and name C5" frame $frame]
			set cua [atomselect $mol "residue $resbottom and name C2" frame $frame]
			#puts "$resbottom"
		} else {
			set cap [atomselect $mol "residue $resbottom and name C2" frame $frame]
			set cua [atomselect $mol "residue $resbottom and name C5" frame $frame]
			#puts "$resbottom"
		}
   	#take coordinates
    	set capcoords [lindex [$cap get {x y}] 0]
    	set cuacoords [lindex [$cua get {x y}] 0]

        #puts "cap_coords: $capcoords"
        #puts "cua_coords: $cuacoords"

#      	puts "$nc(1) $nc(2)"
#        set nc(1) [lindex $capcoords 0]
#        set nc(2) [lindex $cuacoords 0]
#        puts "$nc(1) $nc(2)"
#        set vdist [vecsub $nc(1) $nc(2)] 
      	set vdist [vecsub $capcoords $cuacoords] 
      	#puts "vdist: $vdist"
      	set dist [veclength $vdist] 
      	#puts "dist: $dist"
      	set normz [veclength $unitz]
      	#puts "normz: $normz"
		set dz [vecdot $unitz $vdist]
      	#puts "dz: $dz"
      	set normprod [expr $dist*$normz]
      	#puts "normprod: $normprod"
      	set cosang [expr $dz/$normprod ]
      	#puts "cosinus: $cosang"
      	set ang [expr acos($cosang) ]
      	set angdeg [expr 180.0*$ang/$M_PI]
#    if {$angdeg > 180} {
#       puts "reducing $angdeg ..." 
#       set angdeg [expr $angdeg - 180.0]
#       puts "now $angdeg" 
#    }
		set avcosang [expr $avcosang+$cosang]
		set avang [expr $avang+$angdeg]     
   }
	set timesec [expr $frame*2*10**-6]
    set avcosang [expr $avcosang/$countbottom]
	set avang [expr $avang/$countbottom]
	puts $fpbottom "$frame $timesec $avcosang $avang"
	
# Loop over all selected atoms and save data
#  set count 0
#  foreach r [$sel_N get residue] {
#  incr count
#  set com1 [measure center [atomselect $mol "resname CTA and residue $r" frame $frame] weight mass]
#  puts $fp "$count $com1"
#}

#  foreach coord [$seleccio get {x y z}] {
#    incr count
#    puts $fp "$count $coord"
#    } 
  
}
# ------------------------------------
# Main Program
# ------------------------------------

#open file for output --> MODIFY NAME IF REQUIRED

#Output of top layer data
set fptop [ open "Average_tip4_eq_orientations_TOP_timesec.dat" w ] 
puts $fptop "# Average_Orientation cellulose in top layer"
puts $fptop "# frame cos(angle) angle"

#Output of top layer data
set fpmiddle [ open "Average_tip4_eq_orientations_Middle_timesec.dat" w ] 
puts $fpmiddle "# Average_Orientation cellulose in middle layer"
puts $fpmiddle "# frame cos(angle) angle" 

#Output of bottom layer data
set fpbottom [ open "Average_tip4_eq_orientations_BOTTOM_timesec.dat" w ] 
puts $fpbottom "# Average_Orientation cellulose in bottom layer"
puts $fpbottom "# frame cos(angle) angle"


set unitz {-1.0 1.0}    
#open structure file of simulated system ---> MODIFY NAME OF .psf FILE IF NEEDED
set mol [mol new Tip4_eq_CNC_100ns_aligned.psf type psf waitfor all]

#Define the objects to be selected  --> MODIFY HERE IF NEEDED
#set seleccio [atomselect $mol "resname CTA"]
set sel_N [atomselect $mol "all"]

#Perform calculations over each frame using BigDCD
puts "Calculations ongoing. Please wait..."
#count number of frames

#ATENTION: MODIFY NAMES OF dcd FILES IF NEEDED
bigdcd dothejob auto Tip4_eq_CNC_100ns_aligned.dcd
bigdcd_wait

#close output file          
close $fptop
close $fpmiddle
close $fpbottom 

exit
