import logging
from simtk import unit
import parmed as pm
from simtk.openmm import XmlSerializer, System
from simtk.openmm.app import Simulation
import mdtraj
import numpy as np
from pymbar import mbar
from simtk.openmm.vec3 import Vec3


logger = logging.getLogger(__name__)

def return_reduced_potential(potential_energy:unit.Quantity, volume:unit.Quantity, temperature:unit.Quantity):
    """Retrieve the reduced potential for a given context.
    The reduced potential is defined as in Ref. [1]
    u = \beta [U(x) + p V(x)]
    where the thermodynamic parameters are
    \beta = 1/(kB T) is the inverse temperature
    p is the pressure
    and the configurational properties are
    x the atomic positions
    U(x) is the potential energy
    V(x) is the instantaneous box volume
    References
    ----------
    [1] Shirts MR and Chodera JD. Statistically optimal analysis of
    equilibrium states. J Chem Phys 129:124105, 2008.


    Parameters
    ----------
    potential_energy : simtk.unit of float
    context
    ensamble: NVT or NPT
    """

    assert(type(temperature) == unit.Quantity)
    assert(type(volume) == unit.Quantity)

    pressure = 1.0 * unit.atmosphere # atm      

    beta = 1.0 / (unit.BOLTZMANN_CONSTANT_kB * temperature)
    reduced_potential = potential_energy / unit.AVOGADRO_CONSTANT_NA
    if pressure is not None:
        reduced_potential += pressure * volume
    return beta * reduced_potential



def calculate_energies(env, state, structure_name, current_state=0, conf=None):
    # with modifications taken from:
    # https://github.com/pandegroup/openmm/blob/master/docs-source/usersguide/application.rst#computing-energies
    
    def _energy_over_trajectory(simulation, pos, bxl):
        a = Vec3(bxl.value_in_units_of(unit.nanometer), 0.0, 0.0)
        b = Vec3(0.0, bxl.value_in_units_of(unit.nanometer), 0.0)
        c = Vec3(0.0, 0.0, bxl.value_in_units_of(unit.nanometer))
        simulation.context.setPeriodicBoxVectors(a,b,c)
        simulation.context.setPositions((pos))
        state = simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy()

    def _setup_calculation(psf_file_path, traj_file_path, simulation):
        
        list_e = []
        traj  = mdtraj.load(traj_file_path, top=psf_file_path)
        for idx in range(traj.n_frames):
            bxl = traj.unitcell_lengths[idx][0] * (unit.nanometer)
            e = _energy_over_trajectory(simulation, traj.openmm_positions(idx), bxl) # energy in kJ/mol
            volumn = bxl ** 3
            red_e = return_reduced_potential(e, volumn, 300 * unit.kelvin)
            list_e.append(red_e)
        return list_e
    
    if conf['system']['structure1']['name'] == structure_name:
        structure = 'structure1'
    elif conf['system']['structure2']['name'] == structure_name:
        structure = 'structure2'
    else:
        raise RuntimeError(f"Could not finde structure entry for : {structure_name}")


    conf_sub = conf['system'][structure][env]
    base = f"{conf['system_dir']}/{structure_name}/intst{current_state}/"

    file_name = f"{base}/{conf_sub['intermediate-filename']}_system.xml"
    system  = XmlSerializer.deserialize(open(file_name).read())

    file_name = f"{base}/{conf_sub['intermediate-filename']}_integrator.xml"
    integrator  = XmlSerializer.deserialize(open(file_name).read())

    psf_file_path = f"{base}/{conf_sub['intermediate-filename']}.psf"
    psf = pm.charmm.CharmmPsfFile(psf_file_path)

    simulation = Simulation(psf.topology, system, integrator)
    simulation.context.setState(XmlSerializer.deserialize(open(f"{base}/{conf_sub['intermediate-filename']}/.rst", 'r').read()))
    logger.info('#############')
    logger.info('- Energy evaluation with potential from lambda: {}'.format(str(current_state)))

    traj_file_path = f"{base}/{conf_sub['intermediate-filename']}.dcd"

    logger.info('  - Looking at conformations from lambda: {}'.format(str(state)))
    energy = _setup_calculation(psf_file_path, traj_file_path, simulation)
    return energy
