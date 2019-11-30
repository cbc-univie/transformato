import logging
from simtk import unit
import parmed as pm
from simtk.openmm import XmlSerializer, System
from simtk.openmm.app import Simulation
import mdtraj
import numpy as np
from pymbar import mbar
from simtk.openmm.vec3 import Vec3
import json
from collections import defaultdict, namedtuple
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

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

    
class FreeEnergyCalculator(object):
    
    def __init__(self, configuration:dict, structure_name:str, thinning:int=10, mutation_route:list=[]):
        self.configuration = configuration
        self.structure_name = structure_name
        assert(type(thinning) == int)
        self.thinning = thinning

        # decide if the name of the system corresponds to structure1 or structure2
        if configuration['system']['structure1']['name'] == self.structure_name:
            structure = 'structure1'
        elif configuration['system']['structure2']['name'] == self.structure_name:
            structure = 'structure2'
        else:
            raise RuntimeError(f"Could not finde structure entry for : {self.structure_name}")

        self.base_path = f"{self.configuration['analysis_dir_base']}/{self.structure_name}/"

        self.structure = structure
        self.mutation_route = mutation_route
        self.waterbox_mbar = None
        self.complex_mbar = None
        self.snapshost, self.nr_of_states, self.N_k = self._merge_trajs()
        self._calculate_dG_to_common_core()

    def _merge_trajs(self)->(dict, int):


        #############
        # set all file paths for potential
        nr_of_states = len(next(os.walk(f"{self.base_path}"))[1])
        snapshost = {}
        for env in ['waterbox', 'complex']:
            confs = []
            conf_sub = self.configuration['system'][self.structure][env]
            N_k = []
            for i in range(1, nr_of_states+1):
                traj  = mdtraj.load(f"{self.base_path}/intst{i}/{conf_sub['intermediate-filename']}.dcd", 
                                    top=f"{self.base_path}/intst{i}/{conf_sub['intermediate-filename']}.psf")[50::self.thinning] 
                                    # NOTE: removing the first 50 confs and thinning
                
                if len(traj) < 10:
                    raise RuntimeError(f"Below 10 conformations per lambda -- decrease the thinning factor (currently: {self.thinning}).")
                
                confs.append(traj)
                logger.info(f"Nr of snapshots: {len(traj)}")
                N_k.append(len(traj))
            
            joined_trajs = mdtraj.join(confs, check_topology=True)
            logger.info(f"Combined nr of snapshots: {len(joined_trajs)}")
            snapshost[env] = joined_trajs
        
        return snapshost, nr_of_states, N_k

    def _analyse_results_using_mbar(self, env:str, snapshots:mdtraj.Trajectory, nr_of_states:int):

        def _energy_at_ts(simulation:Simulation, coordinates, bxl:unit.Quantity):
            """
            Calculates the potential energy with the correct periodic boundary conditions.
            """
            a = Vec3(bxl.value_in_unit(unit.nanometer), 0.0, 0.0)
            b = Vec3(0.0, bxl.value_in_unit(unit.nanometer), 0.0)
            c = Vec3(0.0, 0.0, bxl.value_in_unit(unit.nanometer))
            simulation.context.setPeriodicBoxVectors(a,b,c)
            simulation.context.setPositions((coordinates))
            state = simulation.context.getState(getEnergy=True)
            return state.getPotentialEnergy()

        def _evaluated_e_on_all_snapshots(snapshots, i:int, env:str):

            # read in necessary files
            conf_sub = self.configuration['system'][self.structure][env]
            file_name = f"{self.base_path}/intst{i}/{conf_sub['intermediate-filename']}_system.xml"
            system  = XmlSerializer.deserialize(open(file_name).read())
            file_name = f"{self.base_path}/intst{i}/{conf_sub['intermediate-filename']}_integrator.xml"
            integrator  = XmlSerializer.deserialize(open(file_name).read())
            psf_file_path = f"{self.base_path}/intst{i}/{conf_sub['intermediate-filename']}.psf"
            psf = pm.charmm.CharmmPsfFile(psf_file_path)

            # generate simulations object and set states
            simulation = Simulation(psf.topology, system, integrator)
            simulation.context.setState(XmlSerializer.deserialize(open(f"{self.base_path}/intst{i}/{conf_sub['intermediate-filename']}.rst", 'r').read()))      

            energies = []
            for ts in tqdm(range(snapshots.n_frames)):
                # extract the box size at the given ts
                bxl = snapshots.unitcell_lengths[ts][0] * (unit.nanometer)
                # calculate the potential energy 
                e = _energy_at_ts(simulation, snapshots.openmm_positions(ts), bxl)
                # obtain the reduced potential (for NpT)
                volumn = bxl ** 3
                red_e = return_reduced_potential(e, volumn, 300 * unit.kelvin)
                energies.append(red_e)
            return np.array(energies)

        u_kn = np.stack(
                [_evaluated_e_on_all_snapshots(snapshots, i, env) for i in range(1, self.nr_of_states+1)]
                )
        return mbar.MBAR(u_kn, self.N_k)

    def _calculate_dG_to_common_core(self):

        self.waterbox_mbar = self._analyse_results_using_mbar('waterbox', self.snapshost['waterbox'], self.nr_of_states)
        self.complex_mbar =  self._analyse_results_using_mbar('complex', self.snapshost['complex'], self.nr_of_states)


    @property
    def complex_free_energy_differences(self):
        """matrix of free energy differences"""
        return self.complex_mbar.getFreeEnergyDifferences()[0]
    
    @property
    def complex_free_energy_overlap(self):
        """overlap of lambda states"""
        return self.complex_mbar.computeOverlap()[-1]

    @property
    def complex_free_energy_difference_uncertainties(self):
        """matrix of asymptotic uncertainty-estimates accompanying free energy differences"""
        return self.complex_mbar.getFreeEnergyDifferences()[1]
    
    @property
    def waterbox_free_energy_differences(self):
        """matrix of free energy differences"""
        return self.waterbox_mbar.getFreeEnergyDifferences()[0]

    @property
    def waterbox_free_energy_overlap(self):
        """overlap of lambda states"""
        return self.waterbox_mbar.computeOverlap()[-1]
    
    @property
    def waterbox_free_energy_difference_uncertainties(self):
        """matrix of asymptotic uncertainty-estimates accompanying free energy differences"""
        return self.waterbox_mbar.getFreeEnergyDifferences()[1]

    @property
    def plot_complex_free_energy_overlap(self):
        plt.figure(figsize=[8,8], dpi=300)
        plt.imshow(self.complex_free_energy_overlap, cmap='Blues')
        plt.title('Overlap of lambda states for ligand in complex', fontsize=15)
        plt.xlabel('lambda state (0 to 1)', fontsize=15)
        plt.ylabel('lambda state (0 to 1)', fontsize=15)
        plt.legend()
        plt.colorbar()
        plt.show()
        plt.close()    

    @property
    def plot_waterbox_free_energy_overlap(self):
        plt.figure(figsize=[8,8], dpi=300)
        plt.imshow(self.waterbox_free_energy_overlap, cmap='Blues',)
        plt.title('Overlap of lambda states for ligand in waterbox', fontsize=15)
        plt.xlabel('lambda state (0 to 1)', fontsize=15)
        plt.ylabel('lambda state (0 to 1)', fontsize=15)       
        plt.legend()
        plt.colorbar()
        plt.show()
        plt.close()    

    @property
    def plot_complex_free_energy(self):
        x = [a for a in range(1, len(self.complex_free_energy_differences[0])+1)]
        y = self.complex_free_energy_differences[0]
        y_error = self.complex_free_energy_difference_uncertainties[0]

        plt.errorbar(x, y, yerr=y_error, label='ddG +- stddev [kT]')
        plt.title('Free energy estimate for ligand in complex', fontsize=15)
        plt.xlabel('Free energy estimate in kT')
        plt.ylabel('lambda state (0 to 1)', fontsize=15)       
        plt.show()
        plt.close()    

    @property
    def plot_waterbox_free_energy(self):
        x = [a for a in range(1, len(self.waterbox_free_energy_differences[0])+1)]
        y = self.waterbox_free_energy_differences[0]
        y_error = self.waterbox_free_energy_difference_uncertainties[0]

        plt.errorbar(x, y, yerr=y_error, label='ddG +- stddev [kT]')
        plt.title('Free energy estimate for ligand in waterbox', fontsize=15)
        plt.xlabel('Free energy estimate in kT')
        plt.ylabel('lambda state (0 to 1)', fontsize=15)

        plt.show()
        plt.close()    


    @property
    def end_state_free_energy_difference(self):
        """DeltaF[lambda=1 --> lambda=0]"""
        waterbox_DeltaF_ij, waterbox_dDeltaF_ij, _ = self.waterbox_mbar.getFreeEnergyDifferences()
        complex_DeltaF_ij, complex_dDeltaF_ij, _ = self.complex_mbar.getFreeEnergyDifferences()
        K = len(complex_DeltaF_ij)
        return complex_DeltaF_ij[0, K-1] - waterbox_DeltaF_ij[0, K-1], waterbox_dDeltaF_ij[0, K-1] + complex_dDeltaF_ij[0, K-1] 

    def show_summary(self):
        self.plot_complex_free_energy_overlap
        self.plot_waterbox_free_energy_overlap
        self.plot_complex_free_energy
        self.plot_waterbox_free_energy
        energy_estimate, uncertanty = self.end_state_free_energy_difference
        print(f"Free energy to common core: {energy_estimate} [kT] with uncertanty: {uncertanty} [kT].")


