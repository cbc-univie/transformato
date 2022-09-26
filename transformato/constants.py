import logging, sys
from openmm import unit

logger = logging.getLogger(__name__)

temperature = 303.15 * unit.kelvin
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
kT = kB * temperature
charmm_gpu = ""
charmm_gpu = "domdec-gpu"  # uncomment this if you want to use domdec-gpu
this = sys.modules[__name__]
# we can explicitly make assignments on it
this.NUM_PROC = 1

#####################################
# config for tets
test_platform_openMM = "GPU"
test_platform_CHARMM = "CPU"
test_platform_override = True
#####################################


def initialize_NUM_PROC(n_proc):
    if this.NUM_PROC == 1:
        # also in local function scope. no scope specifier like global is needed
        this.NUM_PROC = n_proc
    else:
        msg = "NUM_PROC is already initialized to {0}."
        # raise RuntimeError(msg.format(this.NUM_PROC))
        print(msg)


def change_platform_to_test_platform(configuration: dict, engine: str):

    if engine == "openMM":
        change_to = test_platform_openMM
    elif engine == "CHARMM":
        change_to = test_platform_CHARMM
    else:
        raise RecursionError()

    if change_to.upper() == "GPU":
        configuration["simulation"]["GPU"] = True
        print("Setting platform to GPU")
    elif change_to.upper() == "CPU":
        configuration["simulation"]["GPU"] = False
        print("Setting platform to CPU")
    else:
        raise RuntimeError("something went wrong")


def change_platform_to(configuration: dict, change_to: str):

    if change_to.upper() == "GPU":
        configuration["simulation"]["GPU"] = True
        print("Setting platform to GPU")
    elif change_to.upper() == "CPU":
        configuration["simulation"]["GPU"] = False
        print("Setting platform to CPU")
    else:
        raise RuntimeError("something went wrong")
