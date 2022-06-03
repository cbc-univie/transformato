import datetime
from transformato.constants import temperature, charmm_gpu
from simtk import unit


class CharmmFactory:
    """Class to build the string needed to create a CHARMM input and streaming file"""

    def __init__(self, configuration: dict, structure: str) -> None:

        self.configuration = configuration
        self.structure = structure
        self.vdw_switching_keyword = self._check_switching_function()
        self.charmm_gpu = charmm_gpu

    def _get_simulations_parameters(self):
        prms = {}
        for key in self.configuration["simulation"]["parameters"]:
            prms[key] = self.configuration["simulation"]["parameters"][key]
        return prms

    def _check_switching_function(self) -> str:
        prms = self._get_simulations_parameters()
        vdw = prms.get("vdw", "Force-switch")  # default is vfswitch
        if vdw.lower() == "force-switch":
            return "vfswitch"
        elif vdw.lower() == "switch":
            return "vswitch"
        elif vdw.lower() == "no-switch":  # not implemented
            raise NotImplementedError()
            return ""
        else:
            raise RuntimeError()

    def generate_CHARMM_postprocessing_files(self, env: str) -> str:

        charmm_postprocessing_script = self._get_CHARMM_postprocessing_header(env)
        if env == "vacuum":
            charmm_postprocessing_script += (
                self._get_CHARMM_vacuum_postprocessing_body()
            )
        elif env == "waterbox" or env == "complex":
            charmm_postprocessing_script += (
                self._get_CHARMM_waterbox_postprocessing_body(env)
            )
        else:
            raise NotImplementedError(f"Something went wrong with {env}.")

        return charmm_postprocessing_script

    def generate_CHARMM_production_files(self, env: str) -> str:
        """Body of the CHARMM file with option for gas pahse, waterbox with vswitch and vfswitch"""

        charmm_production_script = self._get_CHARMM_production_header(env)
        if env == "vacuum":
            charmm_production_script += self._get_CHARMM_vacuum_production_body()
        elif env == "waterbox":
            charmm_production_script += self._get_CHARMM_waterbox_production_body(env)

        elif env == "complex":  ###needs to be adaptet from waterbox to complex
            charmm_production_script += self._get_CHARMM_waterbox_production_body(env)
        else:
            raise NotImplementedError(f"Something went wrong with {env}.")

        return charmm_production_script

    @staticmethod
    def build_reduced_toppar(tlc: str) -> str:
        date = datetime.date.today()
        toppar = f"""* Simplified toppar script
* Version from {date}
*

! protein topology and parameter
open read card unit 10 name ../../toppar/top_all36_prot.rtf
read  rtf card unit 10

open read card unit 20 name ../../toppar/par_all36m_prot.prm
read para card unit 20 flex

! nucleic acids
open read card unit 10 name ../../toppar/top_all36_na.rtf
read  rtf card unit 10 append

open read card unit 20 name ../../toppar/par_all36_na.prm
read para card unit 20 append flex

! carbohydrates
open read card unit 10 name ../../toppar/top_all36_carb.rtf
read  rtf card unit 10 append

open read card unit 20 name ../../toppar/par_all36_carb.prm
read para card unit 20 append flex

! lipids
open read card unit 10 name ../../toppar/top_all36_lipid.rtf
read  rtf card unit 10 append

open read card unit 20 name ../../toppar/par_all36_lipid.prm
read para card unit 20 append flex

! CGENFF
open read card unit 10 name ../../toppar/top_all36_cgenff.rtf
read  rtf card unit 10 append

open read card unit 20 name ../../toppar/par_all36_cgenff.prm
read para card unit 20 append flex

! Interface FF
!open read card unit 10 name ../../toppar/top_interface.rtf
!read  rtf card unit 10 append

!open read card unit 10 name ../../toppar/par_interface.prm
!read para card unit 10 append flex

!stream ../../toppar/toppar_all36_nano_lig.str
!stream ../../toppar/toppar_all36_nanolig_patch.str

! Additional topologies and parameters for synthetic polymer
!stream ../../toppar/toppar_all36_synthetic_polymer.str
!stream ../../toppar/toppar_all36_synthetic_polymer_patch.str
!stream ../../toppar/toppar_all36_polymer_solvent.str

! Additional topologies and parameters for water and ions
stream ../../toppar/toppar_water_ions.str
stream ../../toppar/toppar_dum_noble_gases.str
!stream ../../toppar/toppar_ions_won.str

! Additional topologies and parameters for protein
!stream ../../toppar/toppar_all36_prot_arg0.str
stream ../../toppar/toppar_all36_prot_c36m_d_aminoacids.str
stream ../../toppar/toppar_all36_prot_fluoro_alkanes.str
stream ../../toppar/toppar_all36_prot_heme.str
stream ../../toppar/toppar_all36_prot_retinol.str
!stream ../../toppar/toppar_all36_prot_modify_res.str

! Additional topologies and parameters for nucleic acids
stream ../../toppar/toppar_all36_na_nad_ppi.str
!stream ../../toppar/toppar_all36_na_rna_modified.str

! Additional topologies and parameters for lipids
!stream ../../toppar/toppar_all36_lipid_archaeal.str
stream ../../toppar/toppar_all36_lipid_bacterial.str
stream ../../toppar/toppar_all36_lipid_cardiolipin.str
stream ../../toppar/toppar_all36_lipid_cholesterol.str
!stream ../../toppar/toppar_all36_lipid_dag.str
stream ../../toppar/toppar_all36_lipid_inositol.str
!stream ../../toppar/toppar_all36_lipid_lnp.str
stream ../../toppar/toppar_all36_lipid_lps.str
!!stream ../../toppar/toppar_all36_lipid_mycobacterial.str
stream ../../toppar/toppar_all36_lipid_miscellaneous.str
stream ../../toppar/toppar_all36_lipid_model.str
stream ../../toppar/toppar_all36_lipid_prot.str
stream ../../toppar/toppar_all36_lipid_sphingo.str
!stream ../../toppar/toppar_all36_lipid_tag.str
!stream ../../toppar/toppar_all36_lipid_yeast.str
!stream ../../toppar/toppar_all36_lipid_hmmm.str
!stream ../../toppar/toppar_all36_lipid_detergent.str
!stream ../../toppar/toppar_all36_lipid_ether.str
stream ../../toppar/toppar_all36_prot_na_combined.str


! Additional topologies and parameters for carbohydrates
!stream ../../toppar/toppar_all36_carb_glycolipid.str
!stream ../../toppar/toppar_all36_carb_glycopeptide.str
!stream ../../toppar/toppar_all36_carb_imlab.str

! Additional topologies and parameters for spin/fluorophore labels
stream ../../toppar/toppar_all36_label_spin.str
stream ../../toppar/toppar_all36_label_fluorophore.str


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

    def _get_CHARMM_production_header(self, env: str) -> str:
        intermediate_filename = self.configuration["system"][self.structure][env][
            "intermediate-filename"
        ]

        header = f"""*Version September 2020
*Run script for CHARMM jobs from transformato
*

! Read topology and parameter files
stream charmm_toppar.str

! Read PSF
open read unit 10 card name {intermediate_filename}_corr.psf
read psf  unit 10 card

! Read Coordinate
open read unit 10 card name {intermediate_filename}.crd
read coor unit 10 card
        """
        return header

    def _get_CHARMM_vacuum_postprocessing_body(self) -> str:
        body = f"""
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

open read file unit 41 name ../traj.dcd
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

        return body

    def _get_CHARMM_vacuum_production_body(self) -> str:
        ##### gas phase ######
        nstep = self.configuration["simulation"]["parameters"]["nstep"]
        nstout = self.configuration["simulation"]["parameters"]["nstout"]
        nstdcd = self.configuration["simulation"]["parameters"]["nstdcd"]
        switch = self.vdw_switching_keyword

        body = f"""
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

mini sd nstep 100

set nstep = {nstep}
set temp = {temperature.value_in_unit(unit.kelvin)}

scalar fbeta set 5. sele all end
open write unit 21 file name lig_in_vacuum.dcd

DYNA lang leap start time 0.001 nstep @nstep -
    nprint {nstout} iprfrq {nstout} -
    iunread -1 iunwri -1 iuncrd 21 iunvel -1 kunit -1 -
    nsavc {nstdcd} nsavv 0 -
    rbuf 0. tbath @temp ilbfrq 0  firstt @temp -
    ECHECK 0

stop"""
        return body

    def _get_CHARMM_waterbox_postprocessing_body(self, env: str):
        ##### solv phase ######
        switch = self.vdw_switching_keyword
        if (
            self.configuration["simulation"].get("GPU", False) == True
            and self.charmm_gpu == "domdec-gpu"
        ):
            GPU_domdec = "domdec gpu on"
            GPU_openMM = ""
        elif (
            self.configuration["simulation"].get("GPU", False) == True
            and self.charmm_gpu != "domdec-gpu"
        ):
            GPU_openMM = "omm on"
            GPU_domdec = ""
        else:
            GPU_openMM = ""
            GPU_domdec = ""

        body = f"""
stream charmm_{env}_step3_pbcsetup.str

!
! Image Setup
!

open read unit 10 card name charmm_{env}_crystal_image.str
CRYSTAL DEFINE @XTLtype @A @B @C @alpha @beta @gamma
CRYSTAL READ UNIT 10 CARD

!Image centering by residue
IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen sele resname TIP3 end
IMAGE BYSEGI XCEN @xcen YCEN @ycen ZCEN @zcen sele segid pro* end
IMAGE BYSEGI XCEN @xcen YCEN @ycen ZCEN @zcen sele segid het* end

!
! Nonbonded Options
!
! cons fix sele segi solv end

{GPU_openMM} !omm on

nbonds atom vatom {switch} bycb -
       ctonnb 10.0 ctofnb 12.0 cutnb 12.0 cutim 12.0 -
       inbfrq 1 imgfrq 1 wmin 1.0 cdie eps 1.0 -
       ewald pmew fftx @fftx ffty @ffty fftz @fftz  kappa .34 spline order 6

energy
{GPU_domdec}
energy

open read file unit 41 name ../traj.dcd
traj query unit 41

set start ?start
set nframes ?nfile
set skip ?skip

set nframes @nframes !?nfile
traj firstu 41 nunit 1 begi @start skip @skip stop @nframes

open form writ unit 51 name ener_{env}.log
echu 51
set idx 0
label loop
traj read
energy
echo ?ener
incr idx by 1
if @idx .lt. @nframes goto loop


stop"""
        return body

    def _get_CHARMM_waterbox_production_body(self, env):
        ##### waterbox ######
        switch = self.vdw_switching_keyword

        nstep = self.configuration["simulation"]["parameters"]["nstep"]
        nstout = self.configuration["simulation"]["parameters"]["nstout"]
        nstdcd = self.configuration["simulation"]["parameters"]["nstdcd"]
        if (
            self.configuration["simulation"].get("GPU", False) == True
            and self.charmm_gpu == "domdec-gpu"
        ):
            GPU_domdec = "domdec gpu on"
            GPU_openMM = ""
            dyn = """    lang rbuf 0. tbath @temp ilbfrq 0  firstt @temp -
    ECHECK 0"""
            centering = """
calc xcen = 0. 
calc ycen = 0. 
calc zcen = 0.        
"""
        elif (
            self.configuration["simulation"].get("GPU", False) == True
            and self.charmm_gpu != "domdec-gpu"
        ):
            GPU_openMM = "omm on"
            GPU_domdec = ""
            dyn = """    omm langevin gamma 10 firstt @temp finalt @temp -
     prmc pref 1.0 iprsfrq 15"""

        else:
            GPU_openMM = ""
            GPU_domdec = ""
            dyn = """    lang rbuf 0. tbath @temp ilbfrq 0  firstt @temp -
    ECHECK 0"""
        # always center
        centering = """
calc xcen = @A / 2 
calc ycen = @B / 2 
calc zcen = @C / 2        
"""

        body = f"""
!
! Setup PBC (Periodic Boundary Condition)
!

stream charmm_{env}_step3_pbcsetup.str

!
! Image Setup
!

open read unit 10 card name charmm_{env}_crystal_image.str
CRYSTAL DEFINE @XTLtype @A @B @C @alpha @beta @gamma
CRYSTAL READ UNIT 10 CARD

!Image centering by residue
{centering}
IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen sele resname TIP3 end
IMAGE BYRESID XCEN @xcen YCEN @ycen ZCEN @zcen sele segid IONS end
IMAGE BYSEGI XCEN @xcen YCEN @ycen ZCEN @zcen sele segid pro* end
IMAGE BYSEGI XCEN @xcen YCEN @ycen ZCEN @zcen sele segid het* end

!
! Nonbonded Options
!

{GPU_openMM} !omm on

nbonds atom vatom {switch} bycb -
       ctonnb 10.0 ctofnb 12.0 cutnb 16.0 cutim 16.0 -
       inbfrq -1 imgfrq -1 wmin 1.0 cdie eps 1.0 -
       ewald pmew fftx @fftx ffty @ffty fftz @fftz  kappa .34 spline order 6

energy

shak bonh para fast sele segi SOLV end

mini SD nstep 100
mini ABNR nstep 100

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
{GPU_domdec}
energy

set nstep = {nstep}
set temp = {temperature.value_in_unit(unit.kelvin)}

scalar fbeta set 5. sele all end
open write unit 13 file name lig_in_{env}.dcd

DYNA CPT leap start time 0.001 nstep @nstep -
     nprint {nstout} iprfrq {nstout} -
     iunread -1 iunwri -1 iuncrd 13 iunvel -1 kunit -1 -
     nsavc {nstdcd} nsavv 0 -
     PCONSTANT pref   1.0  pmass @Pmass  pgamma   20.0 -
     {dyn}

stop"""
        return body

    def _get_CHARMM_postprocessing_header(self, env: str) -> str:

        intermediate_filename = self.configuration["system"][self.structure][env][
            "intermediate-filename"
        ]

        date = datetime.date.today()
        header = f"""*Version from {date}
*Run script for CHARMM jobs from transformato
*

! Read topology and parameter files
stream charmm_toppar.str

! Read PSF
open read unit 10 card name {intermediate_filename}_corr.psf
read psf  unit 10 card

! Read Coordinate
open read unit 10 card name {intermediate_filename}.crd
read coor unit 10 card"""
        return header
