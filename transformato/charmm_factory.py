import datetime
from transformato.constants import temperature
from simtk import unit


def charmm_factory(configuration: dict, structure: str, env: str) -> str:
    """Function to build the string needed to create a CHARMM input and streaming file"""

    # get env_dir
    intermediate_filename = configuration["system"][structure][env][
        "intermediate-filename"
    ]

    # tlc = configuration["system"][structure]["tlc"]
    nstep = configuration["simulation"]["parameters"]["nstep"]
    print_frq = configuration["simulation"]["parameters"]["nstout"]
    nstdcd = configuration["simulation"]["parameters"]["nstdcd"]

    try:
        GPU = configuration["simulation"]["GPU"]
    except KeyError:
        GPU = False
        pass

    switch = "VFSWItch"  # hard coded switch

    charmm_str = charmm_string(
        env, intermediate_filename, nstep, print_frq, nstdcd, switch, GPU
    )
    return charmm_str


# toppar file
def build_reduced_toppar(tlc: str) -> str:
    date = datetime.date.today()
    toppar = f"""* Simplified toppar script 
* Version from {date} 
*

! Read Protein Topology and Parameter 
open read card unit 10 name ./toppar/top_all36_prot.rtf 
read  rtf card unit 10 
    
open read card unit 20 name ./toppar/par_all36m_prot.prm 
read para card unit 20 flex 

! Read Nucleic Acids 
open read card unit 10 name ./toppar/top_all36_na.rtf 
read  rtf card unit 10 append 
    
open read card unit 20 name ./toppar/par_all36_na.prm 
read para card unit 20 append flex
    
! Read Carbohydrates 
open read card unit 10 name ./toppar/top_all36_carb.rtf 
read  rtf card unit 10 append 
    
open read card unit 20 name ./toppar/par_all36_carb.prm 
read para card unit 20 append flex 

! Read Lipids 
open read card unit 10 name ./toppar/top_all36_lipid.rtf 
read  rtf card unit 10 append 
    
open read card unit 20 name ./toppar/par_all36_lipid.prm 
read para card unit 20 append flex
    
!Read CGENFF 
open read card unit 10 name ./toppar/top_all36_cgenff.rtf 
read  rtf card unit 10 append 
    
open read card unit 20 name ./toppar/par_all36_cgenff.prm 
read para card unit 20 append flex
    
! Additional topologies and parameters for water and ions 
stream ./toppar/toppar_water_ions.str

! Read {tlc} RTF 
open read unit 10 card name {tlc}_g.rtf 
read rtf card unit 10 append

! Read {tlc} prm 
open read unit 10 card name {tlc}.prm 
read para card unit 10 append flex

! Read dummy_atom RTF 
open read unit 10 card name dummy_atom_definitions.rtf 
read rtf card unit 10 append

! Read dummy_atom prm 
open read unit 10 card name dummy_parameters.prm 
read para card unit 10 append flex

"""
    return toppar


def charmm_string(
    env: str,
    intermediate_filename: str,
    nstep: int,
    print_frq: int,
    nstdcd: int,
    switch: str,
    GPU: bool,
):
    """Body of the CHARMM file with option for gas pahse, waterbox with vswitch and vfswitch"""

    if GPU == True:
        GPU = f"""domdec gpu only"""
    else:
        GPU = ""

    header = f"""*Version September 2020 
*Run script for CHARMM jobs from transformato 
*

! Read topology and parameter files 
stream charmm_toppar.str 

! Read PSF 
open read unit 10 card name {intermediate_filename}.psf 
read psf  unit 10 card

! Read Coordinate 
open read unit 10 card name {intermediate_filename}.crd 
read coor unit 10 card
"""

    ##### gas phase ######

    gas_phase = f"""
!coor orie sele all end ! put the molecule at the origin

MMFP
GEO rcm sphere -
    Xref 0.0 Yref 0.0 Zref 0.0 XDIR 1.0 YDIR 1.0 ZDIR 1.0 -
    harmonic FORCE 1.0 select .not. ( hydrogen .or. resname TIP3 ) end
END

set ctofnb 990.
set ctonnb 980.
set cutnb  1000.

nbonds ctonnb @ctonnb ctofnb @ctofnb cutnb @cutnb -
  atom swit vatom {switch} -
  inbfrq 1 

energy   inbfrq 1
energy   inbfrq 0

mini sd nstep 200

set nstep = {nstep} 
set temp = {temperature.value_in_unit(unit.kelvin)}

scalar fbeta set 5. sele all end  
open write unit 21 file name lig_in_vacuum.dcd
 
DYNA lang leap start time 0.001 nstep @nstep -
    nprint {print_frq} iprfrq {print_frq} -
    iunread -1 iunwri -1 iuncrd 21 iunvel -1 kunit -1 -
    nsavc {nstdcd} nsavv 0 -
    rbuf 0. tbath @temp ilbfrq 0  firstt @temp -
    ECHECK 0
    
stop"""

    ##### waterbox ######
    liquid_phase = f"""
!
! Setup PBC (Periodic Boundary Condition)
!

stream charmm_step3_pbcsetup.str

!
! Image Setup
!

open read unit 10 card name charmm_crystal_image.str
CRYSTAL DEFINE @XTLtype @A @B @C @alpha @beta @gamma
CRYSTAL READ UNIT 10 CARD

!Image centering by residue
!IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen sele resname TIP3 end

!
! Nonbonded Options
!

nbonds atom vatom {switch} bycb -
       ctonnb 10.0 ctofnb 12.0 cutnb 16.0 cutim 16.0 -
       inbfrq -1 imgfrq -1 wmin 1.0 cdie eps 1.0 -
       ewald pmew fftx @fftx ffty @ffty fftz @fftz  kappa .34 spline order 6

energy

!
!use a restraint to place center of mass of the molecules near the origin
!

MMFP
GEO rcm sphere -
    Xref @xcen Yref @ycen Zref @zcen XDIR 1.0 YDIR 1.0 ZDIR 1.0 -
    harmonic FORCE 1.0 select .not. ( hydrogen .or. resname TIP3 ) end
END

shak bonh para fast sele segi SOLV end

mini SD nstep 500
mini ABNR nstep 500

!
! NPT dynamics:
! you can change
! nstep  : number of MD steps
! nprint : print-out frequency
! nsavc  : the trajectory saving frequency
!

! estimate Pmass from SYSmass (total system mass)
! [there could be problems with exreme values, such as  Pmass << SYSmass or Pmass >> SYSmass
scalar mass stat
calc Pmass = int ( ?stot  /  50.0 )

energy
{GPU}
energy

set nstep = {nstep}
set temp = {temperature.value_in_unit(unit.kelvin)}

scalar fbeta set 5. sele all end 
open write unit 13 file name lig_in_waterbox.dcd 

DYNA CPT leap start time 0.001 nstep @nstep -
     nprint {print_frq} iprfrq {print_frq} -
     iunread -1 iunwri -1 iuncrd 13 iunvel -1 kunit -1 -
     nsavc {nstdcd} nsavv 0 -
     PCONSTANT pref   1.0  pmass @Pmass  pgamma   20.0 -
     lang rbuf 0. tbath @temp ilbfrq 0  firstt @temp -
     ECHECK 0

stop"""

    if env == "vacuum":
        charmm_str = f"{header}{gas_phase}"
    elif env == "waterbox":
        charmm_str = f"{header}{liquid_phase}"
    else:
        raise RuntimeError(f"Something went wrong. {env} not availalbe.")

    return charmm_str


def charmm_evaluation(
    env: str,
    intermediate_filename: str,
    switch: str,
):

    date = datetime.date.today()
    header = f"""*Version from {date} 
*Run script for CHARMM jobs from transformato 
*

! Read topology and parameter files 
stream charmm_toppar.str 

! Read PSF 
open read unit 10 card name {intermediate_filename}.psf 
read psf  unit 10 card

! Read Coordinate 
open read unit 10 card name {intermediate_filename}.crd 
read coor unit 10 card"""

    ##### gas phase ######

    gas_phase = f"""
coor orie sele all end ! put the molecule at the origin

set ctofnb 990.
set ctonnb 980.
set cutnb  1000.

nbonds ctonnb @ctonnb ctofnb @ctofnb cutnb @cutnb -
  atom swit vatom vswitch -
  inbfrq 1 

energy

energy inbfrq 0

scalar fbeta set 5. sele all end

open read file unit 41 name traj.dcd
traj query unit 41

set start ?start
set nframes ?nfile
set skip ?skip

set nframes @nframes !?nfile
traj firstu 41 nunit 1 begi @start skip @skip stop @nframes

open form writ unit 51 name ener_vac.log
echu 51
set idx 0
label loop
traj read
energy
echo ?ener
incr idx by 1
if @idx .lt. @nframes goto loop
   
  
stop"""

    ##### solv phase ######

    solv_phase = f"""
stream charmm_step3_pbcsetup.str

!
! Image Setup
!

open read unit 10 card name charmm_crystal_image.str
CRYSTAL DEFINE @XTLtype @A @B @C @alpha @beta @gamma
CRYSTAL READ UNIT 10 CARD

!Image centering by residue
IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen sele resname TIP3 end

!
! Nonbonded Options
!
! cons fix sele segi solv end

nbonds atom vatom {switch} bycb -
       ctonnb 10.0 ctofnb 12.0 cutnb 12.0 cutim 12.0 -
       inbfrq 1 imgfrq 1 wmin 1.0 cdie eps 1.0 -
       ewald pmew fftx @fftx ffty @ffty fftz @fftz  kappa .34 spline order 6

energy

!
!use a restraint to place center of mass of the molecules near the origin
!

open read file unit 41 name traj.dcd
traj query unit 41

set start ?start
set nframes ?nfile
set skip ?skip

set nframes @nframes !?nfile
traj firstu 41 nunit 1 begi @start skip @skip stop @nframes

open form writ unit 51 name ener_solv.log
echu 51
set idx 0
label loop
traj read
energy
echo ?ener
incr idx by 1
if @idx .lt. @nframes goto loop
        
  
stop"""

    if env == "vacuum":
        charmm_evaluation_str = f"{header}{gas_phase}"
    elif env == "waterbox":
        charmm_evaluation_str = f"{header}{solv_phase}"
    else:
        raise RuntimeError(f"Something went wrong. {env} not availalbe.")

    return charmm_evaluation_str