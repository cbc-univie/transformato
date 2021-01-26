import logging
from simtk import unit

logger = logging.getLogger(__name__)

platform = "cpu"
temperature = 303.15 * unit.kelvin
kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
kT = kB * temperature
