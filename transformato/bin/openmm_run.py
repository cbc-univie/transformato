"""
Generated by CHARMM-GUI (http://www.charmm-gui.org)

openmm_run.py

This program is OpenMM running scripts written in python.

Correspondance: jul316@lehigh.edu or wonpil@lehigh.edu
Last update: March 29, 2017
"""

from __future__ import print_function
import argparse
import sys
import os

from omm_readinputs import *
from omm_readparams import *
from omm_vfswitch import *
from omm_barostat import *
from omm_restraints import *
from omm_rewrap import *

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="inpfile", help="Input parameter file", required=True)
parser.add_argument("-p", dest="psffile", help="Input CHARMM PSF file", required=True)
parser.add_argument("-c", dest="crdfile", help="Input CHARMM CRD file", required=True)
parser.add_argument(
    "-t", dest="toppar", help="Input CHARMM-GUI toppar stream file", required=True
)
parser.add_argument(
    "-b",
    dest="sysinfo",
    help="Input CHARMM-GUI sysinfo stream file (optional)",
    default=None,
)
parser.add_argument(
    "-icrst",
    metavar="RSTFILE",
    dest="icrst",
    help="Input CHARMM RST file (optional)",
    default=None,
)
parser.add_argument(
    "-irst",
    metavar="RSTFILE",
    dest="irst",
    help="Input restart file (optional)",
    default=None,
)
parser.add_argument(
    "-ichk",
    metavar="CHKFILE",
    dest="ichk",
    help="Input checkpoint file (optional)",
    default=None,
)
parser.add_argument(
    "-opdb",
    metavar="PDBFILE",
    dest="opdb",
    help="Output PDB file (optional)",
    default=None,
)
parser.add_argument(
    "-orst",
    metavar="RSTFILE",
    dest="orst",
    help="Output restart file (optional)",
    default=None,
)
parser.add_argument(
    "-ochk",
    metavar="CHKFILE",
    dest="ochk",
    help="Output checkpoint file (optional)",
    default=None,
)
parser.add_argument(
    "-odcd",
    metavar="DCDFILE",
    dest="odcd",
    help="Output trajectory file (optional)",
    default=None,
)
parser.add_argument(
    "-rewrap",
    dest="rewrap",
    help="Re-wrap the coordinates in a molecular basis (optional)",
    action="store_true",
    default=False,
)
args = parser.parse_args()

# Load parameters
print("Loading parameters")
inputs = read_inputs(args.inpfile)
params = read_params(args.toppar)
psf = CharmmPsfFile(args.psffile)
crd = read_crd(args.crdfile)
if args.sysinfo:
    psf = read_box(psf, args.sysinfo)
else:
    psf = gen_box(psf, crd)

# Build system
nboptions = dict(
    nonbondedMethod=inputs.coulomb,
    nonbondedCutoff=inputs.r_off * nanometers,
    constraints=inputs.cons,
    ewaldErrorTolerance=inputs.ewald_Tol,
)

if inputs.vdw == "Switch":
    nboptions["switchDistance"] = inputs.r_on * nanometers
system = psf.createSystem(params, **nboptions)
if inputs.vdw == "Force-switch":
    system = vfswitch(system, psf, inputs)

system = barostat(system, inputs)
if inputs.rest == "yes":
    system = restraints(system, crd, inputs)
integrator = LangevinIntegrator(
    inputs.temp * kelvin, inputs.fric_coeff / picosecond, inputs.dt * picoseconds
)

# Set platform
platform = Platform.getPlatformByName("CUDA")
prop = dict()

# Check if restraints.yaml exists - if it does, system uses restraints

pdbpath=args.inpfile.replace(".inp",".pdb")

if os.path.exists("./restraints.yaml") and "complex" in pdbpath:
    import transformato.restraints as tfrs
    import yaml

    print("Found restraints.yaml - applying restraints")
    # Load tiny restraints config
    with open("./restraints.yaml","r") as stream:
        try:
            configuration=yaml.safe_load(stream)
        except yaml.YAMLError as exc:
                print(exc)

    cc_names=configuration["system"]["structure"]["ccs"]

    # Add forces via transformato.restraints
    
    if not os.path.exists(pdbpath):
        raise FileNotFoundError(f"Couldnt find {pdbpath} necessary for Restraint Analysis")

    
    restraintList=tfrs.create_restraints_from_config(configuration,pdbpath)

    for restraint in restraintList:
        restraint.createForce(cc_names)
        restraint.applyForce(system)
        




# Build simulation context
simulation = Simulation(psf.topology, system, integrator, platform, prop)
simulation.context.setPositions(crd.positions)
if args.icrst:
    charmm_rst = read_charmm_rst(args.icrst)
    simulation.context.setPositions(charmm_rst.positions)
    simulation.context.setVelocities(charmm_rst.velocities)
    simulation.context.setPeriodicBoxVectors(
        charmm_rst.box[0], charmm_rst.box[1], charmm_rst.box[2]
    )
if args.irst:
    with open(args.irst, "r") as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))
if args.ichk:
    with open(args.ichk, "rb") as f:
        simulation.context.loadCheckpoint(f.read())

# Re-wrap
if args.rewrap:
    simulation = rewrap(simulation)

# Calculate initial system energy
print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Energy minimization
if inputs.mini_nstep > 0:
    print("\nEnergy minimization: %s steps" % inputs.mini_nstep)
    simulation.minimizeEnergy(
        tolerance=inputs.mini_Tol * kilojoule / mole, maxIterations=inputs.mini_nstep
    )
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Generate initial velocities
if inputs.gen_vel == "yes":
    print("\nGenerate initial velocities")
    if inputs.gen_seed:
        simulation.context.setVelocitiesToTemperature(inputs.gen_temp, inputs.gen_seed)
    else:
        simulation.context.setVelocitiesToTemperature(inputs.gen_temp)

# Production
if inputs.nstep > 0:
    print("\nMD run: %s steps" % inputs.nstep)
    if inputs.nstdcd > 0:
        if not args.odcd:
            args.odcd = "output.dcd"
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
    # Simulated annealing?
    if inputs.annealing == "yes":
        interval = inputs.interval
        temp = inputs.temp_init
        for i in range(inputs.nstep):
            integrator.setTemperature(temp * kelvin)
            simulation.step(1)
            temp += interval
    else:
        simulation.step(inputs.nstep)

# Write restart file
if not (args.orst or args.ochk):
    args.orst = "output.rst"
if args.orst:
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(args.orst, "w") as f:
        f.write(XmlSerializer.serialize(state))
if args.ochk:
    with open(args.ochk, "wb") as f:
        f.write(simulation.context.createCheckpoint())
if args.opdb:
    crd = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(psf.topology, crd, open(args.opdb, "w"))

