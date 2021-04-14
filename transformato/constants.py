import logging
from simtk import unit

logger = logging.getLogger(__name__)

platform = "CPU"
temperature = 303.15 * unit.kelvin
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
kT = kB * temperature


def check_platform(configuration: dict):

    if platform.upper() == "GPU":
        configuration["simulation"]["GPU"] = True
        print("Setting platform to GPU")
    elif platform.upper() == "CPU":
        configuration["simulation"]["GPU"] = False
        print("Setting platform to CPU")
    else:
        raise RuntimeError("something went wrong")