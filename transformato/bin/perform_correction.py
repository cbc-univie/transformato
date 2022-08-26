# general imports
from endstate_correction.system import create_charmm_system, read_box
from openmm.app import (
    PME,
    CharmmParameterSet,
    CharmmPsfFile,
    PDBFile,
)
from endstate_correction.analysis import plot_endstate_correction_results
import endstate_correction
from endstate_correction.protocoll import perform_endstate_correction, Protocoll
import mdtraj
from openmm import unit

########################################################
########################################################
# ------------ set up the waterbox system --------------


system_name = "ethane"
env = "vacuum"
# define the output directory
output_base = f"."
parameter_base = f"/site/raid3/johannes/free_solv_test/data"
# load the charmm specific files (psf, pdb, rtf, prm and str files)
psf_file = f"../methanol/intst1/lig_in_{env}.psf"
psf = CharmmPsfFile(psf_file)
pdb = PDBFile(f"../methanol/intst1/lig_in_{env}.pdb")
params = CharmmParameterSet(
    f"{parameter_base}/{system_name}/waterbox/unk/unk.str",
    f"{parameter_base}/{system_name}/waterbox/toppar/top_all36_cgenff.rtf",
    f"{parameter_base}/{system_name}/waterbox/toppar/par_all36_cgenff.prm",
    f"{parameter_base}/{system_name}/waterbox/toppar/toppar_water_ions.str",
)
# set up the treatment of the system for the specific environment
if env == "waterbox":
    psf = read_box(psf, f"{parameter_base}/{system_name}/waterbox/openmm/sysinfo.dat")

# define region that should be treated with the qml
chains = list(psf.topology.chains())
ml_atoms = [atom.index for atom in chains[0].atoms()]
# define system
sim = create_charmm_system(psf=psf, parameters=params, env=env, ml_atoms=ml_atoms)

########################################################
########################################################
# ------------------- load samples ---------------------#
n_samples = 500
n_steps_per_sample = 1250000
traj_base = f"../methanol/intst1/"
mm_samples = []
traj = mdtraj.load_dcd(
    f"{traj_base}/run_1/lig_in_{env}.dcd",
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
