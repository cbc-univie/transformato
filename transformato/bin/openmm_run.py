"""
File created in analogy to CHARMM-GUI (http://www.charmm-gui.org)
Last update: November, 2023
"""

from __future__ import print_function
import argparse
import sys
import os

from omm_readinputs import *
from omm_readparams import *
from omm_vfswitch import *

import openmm.unit as unit
from openmm.unit import *
from openmm import *
from openmm.app import *

parser = argparse.ArgumentParser()
parser.add_argument("-odcd", metavar="DCDFILE", dest="odcd")
parser.add_argument("-env", metavar="ENVIRONMENT", dest="env")
args = parser.parse_args()

# Load parameters
env = args.env
print(f"Loading parameters in this environment {env}")


inputs = read_inputs(f"lig_in_{env}.inp")

if os.path.isfile(f"lig_in_{env}.parm7"):
    top = AmberPrmtopFile(f"lig_in_{env}.parm7")
    crd = AmberInpcrdFile(f"lig_in_{env}.rst7")
    fftype = "amber"
else:
    fftype = "charmm"
    params = read_params("toppar.str")
    top = CharmmPsfFile(f"lig_in_{env}.psf")
    crd = read_crd(f"lig_in_{env}.crd")
    top = gen_box(top, crd)

# Build system
if env == "waterbox":
    nboptions = dict(
        nonbondedMethod=inputs.coulomb,
        nonbondedCutoff=inputs.r_off * unit.nanometers,
        constraints=inputs.cons,
        ewaldErrorTolerance=inputs.ewald_Tol,
    )
elif env == "vacuum":
    nboptions = dict(
        nonbondedMethod=NoCutoff,
        constraints=inputs.cons,
    )

if inputs.vdw == "Switch" and env != "vacuum":
    nboptions["switchDistance"] = inputs.r_on * unit.nanometers
if fftype == "amber":
    system = top.createSystem(**nboptions)
else:
    system = top.createSystem(params, **nboptions)


if inputs.vdw == "Force-switch" and fftype != "amber" and env != "vacuum":
    system = vfswitch(system, top, inputs)

if env != "vacuum":
    barostat = MonteCarloBarostat(inputs.p_ref * bar, inputs.temp * kelvin)
    system.addForce(barostat)

integrator = LangevinIntegrator(
    inputs.temp * kelvin, 1 / unit.picosecond, inputs.dt * unit.picoseconds
)

# Set platform
platform = Platform.getPlatformByName("CUDA")
prop = dict(CudaPrecision="mixed")

# Check if restraints.yaml exists - if it does, system uses restraints

# pdbpath = args.inpfile.replace(".inp", ".pdb")

# if os.path.exists("./restraints.yaml") and "complex" in pdbpath:
#     import transformato.restraints as tfrs
#     import yaml

#     print("Found restraints.yaml - applying restraints")
#     # Load tiny restraints config
#     with open("./restraints.yaml", "r") as stream:
#         try:
#             configuration = yaml.safe_load(stream)
#         except yaml.YAMLError as exc:
#             print(exc)

#     cc_names = configuration["system"]["structure"]["ccs"]

#     # Add forces via transformato.restraints

#     if not os.path.exists(pdbpath):
#         raise FileNotFoundError(
#             f"Couldnt find {pdbpath} necessary for Restraint Analysis"
#         )

#     restraintList = tfrs.create_restraints_from_config(configuration, pdbpath)

#     for restraint in restraintList:
#         restraint.createForce(cc_names)
#         restraint.applyForce(system)


# Build simulation context
simulation = Simulation(top.topology, system, integrator, platform, prop)
simulation.context.setPositions(crd.positions)
if os.path.isfile(f"lig_in_{env}.irst"):
    with open(f"lig_in_{env}.irst", "r") as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))


# Calculate initial system energy
print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Energy minimization
if inputs.mini_nstep > 0:
    print("\nEnergy minimization: %s steps" % inputs.mini_nstep)
    simulation.minimizeEnergy(maxIterations=inputs.mini_nstep)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Generate initial velocities
if inputs.gen_vel == "yes":
    print("\nGenerate initial velocities")
    if inputs.gen_seed:
        simulation.context.setVelocitiesToTemperature(inputs.gen_temp, inputs.gen_seed)
    else:
        simulation.context.setVelocitiesToTemperature(inputs.gen_temp)

# Production
print("\nMD run: %s steps" % inputs.nstep)
simulation.reporters.append(DCDReporter(args.odcd, inputs.nstdcd))
simulation.reporters.append(
    StateDataReporter(
        sys.stdout,
        inputs.nstout,
        step=True,
        time=True,
        potentialEnergy=True,
        temperature=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=inputs.nstep,
        separator="\t",
    )
)

simulation.step(inputs.nstep)

# needed for later analysis
file_name = f"lig_in_{env}"
print(file_name)
with open(file_name + "_integrator.xml", "w") as outfile:
    outfile.write(XmlSerializer.serialize(integrator))
with open(file_name + "_system.xml", "w") as outfile:
    outfile.write(XmlSerializer.serialize(system))
