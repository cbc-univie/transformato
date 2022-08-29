# general imports
from os import system
from endstate_correction.system import create_charmm_system, read_box
from openmm.app import (
    PME,
    CharmmParameterSet,
    CharmmPsfFile,
    PDBFile,
    CharmmCrdFile,
)
from endstate_correction.analysis import plot_endstate_correction_results
import endstate_correction
from endstate_correction.protocoll import perform_endstate_correction, Protocoll
import mdtraj
from openmm import unit
import sys

########################################################
########################################################
# ------------ set up the waterbox system --------------

system_name = sys.argv[1]

def gen_box(psf, crd):
    coords = crd.positions

    min_crds = [coords[0][0], coords[0][1], coords[0][2]]
    max_crds = [coords[0][0], coords[0][1], coords[0][2]]

    for coord in coords:
        min_crds[0] = min(min_crds[0], coord[0])
        min_crds[1] = min(min_crds[1], coord[1])
        min_crds[2] = min(min_crds[2], coord[2])
        max_crds[0] = max(max_crds[0], coord[0])
        max_crds[1] = max(max_crds[1], coord[1])
        max_crds[2] = max(max_crds[2], coord[2])

    boxlx = max_crds[0]-min_crds[0]
    boxly = max_crds[1]-min_crds[1]
    boxlz = max_crds[2]-min_crds[2]

    psf.setBox(boxlx, boxly, boxlz)
    return psf


for env in ["waterbox","vacuum"]:

    parameter_base = f"/site/raid3/johannes/free_solv_test/data"
    # load the charmm specific files (psf, pdb, rtf, prm and str files)
    psf_file = f"../{system_name}/intst1/lig_in_{env}.psf"
    psf = CharmmPsfFile(psf_file)
    pdb = PDBFile(f"../{system_name}/intst1/lig_in_{env}.pdb")
    crd = CharmmCrdFile(f"../{system_name}/intst1/lig_in_{env}.crd")
    params = CharmmParameterSet(
        f"../{system_name}/intst1/waterbox/unk/unk.str",
        f"../toppar/top_all36_cgenff.rtf",
        f"../toppar/par_all36_cgenff.prm",
        f"../toppar/toppar_water_ions.str",
    )
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

    mm_samples = []
    traj = mdtraj.load_dcd(
        f"../{system_name}/intst1/run_1/lig_in_{env}.dcd",
        top=psf_file,
    )
    if env == "waterbox":
        traj.image_molecules()
        
    mm_samples.extend(traj.xyz * unit.nanometer)  # NOTE: this is in nanometer!


    ####################################################
    # ----------------------- FEP ----------------------
    ####################################################

    fep_protocoll = Protocoll(
        method="NEQ",
        direction="unidirectional",
        sim=sim,
        trajectories=[mm_samples],
        nr_of_switches=500, # 500
        neq_switching_length=5000, # 5000
    )

    r = perform_endstate_correction(fep_protocoll)
    plot_endstate_correction_results(system_name, r, "results_neq_unidirectional.png")
