import argparse

from omm_readinputs import *
from omm_readparams import *
from openmm.unit import *
from openmm import *
from openmm.app import *

parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="inpfile", help="Input parameter file", required=True)
parser.add_argument("-p", dest="psffile", help="Input CHARMM PSF file", required=True)
parser.add_argument("-c", dest="crdfile", help="Input CHARMM CRD file", required=True)
parser.add_argument(
    "-t", dest="toppar", help="Input CHARMM-GUI toppar stream file", required=True
)
parser.add_argument(
    "-irst",
    metavar="RSTFILE",
    dest="irst",
    help="Input restart file (optional)",
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
    "-odcd",
    metavar="DCDFILE",
    dest="odcd",
    help="Output trajectory file (optional)",
    default=None,
)
args = parser.parse_args()

inputs = read_inputs(args.inpfile)
params = read_params(args.toppar)
psf = CharmmPsfFile(args.psffile)
crd = read_crd(args.crdfile)

# Build system
integrator = LangevinIntegrator(
    inputs.temp * kelvin, inputs.fric_coeff / picosecond, inputs.dt * picoseconds
)

nboptions = dict(
    nonbondedMethod=NoCutoff,
    constraints=inputs.cons,
)

system = psf.createSystem(params, **nboptions)

simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(crd.positions)

print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.minimizeEnergy(maxIterations=500)
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

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

if not (args.orst):
    args.orst = "output.rst"
if args.orst:
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(args.orst, "w") as f:
        f.write(XmlSerializer.serialize(state))
