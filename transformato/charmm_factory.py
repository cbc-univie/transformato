import os
import datetime
from transformato.utils import get_bin_dir,get_toppar_dir,load_config_yaml

def parser(string, name):
        file_path = os.getcwd()
        file_name = name
        tmp_path = file_path + file_name
        
        try:
            with open(tmp_path, 'w') as f:
                f.write(string)
                f.close()
        except IOError:
            print(f"CHARMM input for {file_name} could not be created.")
            #logger.info(f"Data class could not be created: {file_name}")
            pass

def charmm_factory(configuration,structure):
    """Function to build the string needed to create a CHARMM input and streaming file"""

    tlc = configuration['system'][structure]['tlc']
    vacuum = configuration['system'][structure]['vacuum']['intermediate-filename']
    waterbox = configuration['system'][structure]['waterbox']['intermediate-filename']
    nstep = configuration['simulation']['parameters']['nstep']
    nstout = configuration['simulation']['parameters']['nstout']
    nstdcd = configuration['simulation']['parameters']['nstdcd']
    steps_for_equilibration = configuration['solvation']['steps_for_equilibration']
    
    #building a reduced toppar file and including dummy rtf and prm
    toppar = build_reduced_toppar(tlc)
    parser(toppar,'/toppar_with_dummy.str')

    #gas phase
    env = 'vacuum'
    header = header_string(vacuum)
    body = body_string(env,nstep,nstout,nstdcd,steps_for_equilibration)
    gas_phase_file = f'{header}{body}'
    parser(gas_phase_file,'/run_gasp_md.inp')

    #waterbox
    env = 'waterbox'
    header = header_string(waterbox)
    body_vswi, body_vfswi = body_string(env,nstep,nstout,nstdcd,steps_for_equilibration)
    waterbox_vswitch = f'{header}{body_vswi}'
    waterbox_vfswitch = f'{header}{body_vfswi}'
    parser(waterbox_vswitch,'/run_liqp_md_vswi.inp')
    parser(waterbox_vfswitch,'/run_liqp_md_vfsw.inp')


#toppar file
def build_reduced_toppar(tlc):    
    toppar = """
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

"""
    
    rtf = f'! Read {tlc} RTF \nopen read unit 10 card name {tlc}_g.rtf \nread rtf unit 10 append \n\n'
    prm = f'! Read {tlc} prm \nopen read unit 10 card name {tlc}.prm \nread para unit 10 append flex \n\n'
    dummy_rtf = f'! Read dummy_atom RTF \nopen read unit 10 card name dummy_atom_definitions.rtf \nread rtf unit 10 append \n\n'
    dummy_prm = f'! Read dummy_atom prm \nopen read unit 10 card name dummy_parameters.prm \nread para unit 10 append flex \n\n'
    date = datetime.date.today()
    toppar_file = f'* Simplified toppar script \n* Version from {date} \n* \n\n{toppar}{rtf}{prm}{dummy_rtf}{dummy_prm}'
    
    return toppar_file

def header_string(type): 
    """Header of the CHARMM file with flexible psf and crd input"""   
    version = '*Version September 2020 \n*Run script for CHARMM jobs from transformato \n*\n\n'
    streaming_file = '! Read topology and parameter files \nstream toppar_with_dummy.str \n\n'
    psf = f'! Read PSF \nopen read unit 10 card name {type}.psf \nread psf  unit 10 card\n\n'
    crd = f'! Read Coordinate \nopen read unit 10 card name {type}.crd \nread coor unit 10 card'
    header = f'{version}{streaming_file}{psf}{crd}'
    return header

def  body_string(env,nstep,nstout,nstdcd,steps_for_equilibration): 
    """Body of the CHARMM file with option for gas pahse, waterbox with vswitch and vfswitch"""

    ##### gas phase ######
    gas_phase_1 = """
coor orie sele all end ! put the molecule at the origin

MMFP
GEO rcm sphere -
    Xref 0.0 Yref 0.0 Zref 0.0 XDIR 1.0 YDIR 1.0 ZDIR 1.0 -
    harmonic FORCE 1.0 select .not. ( hydrogen .or. resname TIP3 ) end
END

set ctofnb 990.
set ctonnb 980.
set cutnb  1000.

nbonds ctonnb @ctonnb ctofnb @ctofnb cutnb @cutnb -
  atom swit vatom vswitch -
  inbfrq 1 

energy   inbfrq 1
energy   inbfrq 0

mini sd nstep 10

if @?rand .eq. 0 set rand 1

"""
    gas_phase_2 = f'set nstep = {round(nstep/10)} \nset temp = 300.0'
    gas_phase_3 = """
open write unit 12 card name gas_equil1.@rand.rst
scalar fbeta set 5. sele all end

"""
    gas_phase_4 = f'DYNA leap lang start time 0.001 nstep @nstep -\n    nprint {nstout} iprfrq @nstep -\n    iunread -1 iunwri 12 iuncrd -1 iunvel -1 kunit -1 -\n    nsavc {round(nstep/40)} nsavv 0 iseed @rand @rand @rand @rand -\n    rbuf 0. tbath @temp ilbfrq 0  echeck 0 firstt 250. ! dont overshoot\n'
                    
    gas_phase_5 = """
open read  unit 11 card name gas_equil1.@rand.rst
open write unit 12 card name gasp.@rand.rst
open write unit 21 file name gasp.@rand.dcd

"""
    gas_phase_6 = f'set nstep = {nstep} \nDYNA lang leap restart time 0.001 nstep @nstep -\n    nprint {nstout} iprfrq {round(nstep/20)} -\n    iunread 11 iunwri 12 iuncrd 21 iunvel -1 kunit -1 -\n    nsavc {nstdcd} nsavv 0 -\n    rbuf 0. tbath @temp ilbfrq 0  firstt @temp -\n    echeck 0\n\nstop'  
  
    ##### waterbox ###### 
    liquid_phase_1 ="""
set bxl 30.0

MMFP
GEO rcm sphere -
    Xref 0.0 Yref 0.0 Zref 0.0 XDIR 1.0 YDIR 1.0 ZDIR 1.0 -
    harmonic FORCE 1.0 select .not. ( hydrogen .or. resname TIP3 ) end
END

!shak bonh para fast sele segi WAT end
shak bonh para fast sele segi SOLV end

set ctofnb 12.
set ctonnb 10.
set cutnb  14.
set cutim  @cutnb

crystal define cubic @bxl @bxl @bxl 90.00 90.00 90.00
crystal build Noper 0 cutoff @cutim
image byres xcen 0.0 ycen 0.0 zcen 0.0 sele all end

nbonds ctonnb @ctonnb ctofnb @ctofnb cutnb @cutnb cutim @cutim -
"""
    VSWI = '  atom shif vatom VSWI -\n'
    VFSW = '  atom shif vatom VFSW -\n'
    liquid_phase_2 ="""  inbfrq -1 imgfrq -1 -
  ewald pmew kappa 0.34 spline order 6 fftx 32 ffty 32 fftz 32

energy
domdec gpu only
energy

mini sd nstep 10

! from charmm-gui scripts -- CPT dynamics
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

if @?rand .eq. 0 set rand 1

"""
    liquid_phase_3 = f'set nstep = {round(nstep/10)} \nset temp = 300.0'
    liquid_phase_4 = """
    
open write unit 12 card name equilbox_@{method}.@{rand}.rst
scalar fbeta set 5. sele all end

"""
    liquid_phase_5 = f'DYNA CPT leap start time 0.001 nstep @nstep -\n    nprint {steps_for_equilibration} iprfrq {steps_for_equilibration} ntrfrq {steps_for_equilibration} -\n    iunread -1 iunwri 12 iuncrd -1 iunvel -1 kunit -1 -\n    nsavc {round(nstep/40)} nsavv 0 iseed @rand @rand @rand @rand -\n    PCONSTANT pref   1.0  pmass @Pmass  pgamma   20.0 -\n    HOOVER reft @temp  tmass 2000.0  tbath   @temp  firstt @temp\n'
    liquid_phase_6 = """   
open read  unit 11 card name equilbox_@{method}.@{rand}.rst
open write unit 12 card name prod_@{method}.@{rand}.rst
open write unit 21 file name prod_@{method}.@{rand}.dcd 

"""
    liquid_phase_7 = f'DYNA lang leap restart time 0.001 nstep @nstep -\n    nprint {steps_for_equilibration} iprfrq {round(nstep/200)} -\n    iunread 11 iunwri 12 iuncrd 21 iunvel -1 kunit -1 -\n    nsavc {nstdcd} nsavv 0 -\n    PCONSTANT pref 1.0  pmass @Pmass  pgamma   20.0 -\n    lang rbuf 0. tbath @temp ilbfrq 0  firstt @temp -\n    echeck 0\n\nstop'
    


    if env == 'vacuum':
        body = f'{gas_phase_1}{gas_phase_2}{gas_phase_3}{gas_phase_4}{gas_phase_5}{gas_phase_6}'
        return body
    elif env == 'waterbox':
        body_vswitch = f'{liquid_phase_1}{VSWI}{liquid_phase_2}{liquid_phase_3}{liquid_phase_4}{liquid_phase_5}{liquid_phase_6}{liquid_phase_7}'
        body_vfswitch = f'{liquid_phase_1}{VFSW}{liquid_phase_2}{liquid_phase_3}{liquid_phase_4}{liquid_phase_5}{liquid_phase_6}{liquid_phase_7}'
        return body_vswitch, body_vfswitch


#testsuit
configuration = load_config_yaml(config='transformato/tests/config/test-toluene-methane-solvation-free-energy.yaml',
                                    input_dir='data/', output_dir='.')

charmm_factory(configuration,'structure1')