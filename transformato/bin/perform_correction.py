# general imports
import os
import glob
from endstate_correction.system import create_charmm_system, gen_box
from openmm.app import (
    CharmmParameterSet,
    CharmmPsfFile,
    PDBFile,
    CharmmCrdFile,
)
from endstate_correction.analysis import plot_endstate_correction_results
from endstate_correction.protocol import perform_endstate_correction, Protocol
import mdtraj
import pickle
from openmm import unit


# Variables will be generated by transformato
system_name = NAMEofSYSTEM
tlc = TLC

for env in ["waterbox", "vacuum"]:


    # load the charmm specific files (psf, rtf, crd files)
    psf_file = f"../{system_name}/intst1/lig_in_{env}.psf"
    psf = CharmmPsfFile(psf_file)
    pdb = PDBFile(f"../{system_name}/intst1/lig_in_{env}.pdb")
    crd = CharmmCrdFile(f"../{system_name}/intst1/lig_in_{env}.crd")

    # load forcefiled files (ligand.str and toppar files)
    parms = ()
    file = f"../{system_name}/intst1/{tlc.lower()}"
    if os.path.isfile(f"{file}.str"):
        parms += (f"{file}.str",)
    else:
        parms += (f"{file}_g.rtf",)
        parms += (f"{file}.prm",)

    parms += (f"../toppar/top_all36_cgenff.rtf",)
    parms += (f"../toppar/par_all36_cgenff.prm",)
    parms += (f"../toppar/toppar_water_ions.str",)

    params = CharmmParameterSet(*parms)

    # set up the treatment of the system for the specific environment
    if env == "waterbox":
        psf = gen_box(psf, crd)

    # define region that should be treated with the qml
    chains = list(psf.topology.chains())
    ml_atoms = [atom.index for atom in chains[0].atoms()]
    # define system
    sim = create_charmm_system(psf=psf, parameters=params, env=env, ml_atoms=ml_atoms)

    ########################################################
    ########################################################
    # ------------------- load samples ---------------------#

    files = glob.glob(f"../{system_name}/intst1/**/lig_in_{env}.dcd", recursive=True)
    traj = mdtraj.load(files, top=psf_file)

    if env == "waterbox":
        traj.image_molecules()

    mm_samples = []
    mm_samples.extend(traj.xyz * unit.nanometer)  # NOTE: this is in nanometer!

    ####################################################
    # ----------------------- FEP ----------------------
    ####################################################

    fep_protocoll = Protocol(
        method="NEQ",
        direction="unidirectional",
        sim=sim,
        trajectories=[mm_samples],
        nr_of_switches=500,  # 500
        neq_switching_length=5000,  # 5000
    )

    r = perform_endstate_correction(fep_protocoll)
    
    with open(f"results_{system_name}.pickle", "wb") as file:
        pickle.dump(r, file)
        
    plot_endstate_correction_results(
        system_name, r, f"results_neq_unidirectional_{env}.png"
    )
