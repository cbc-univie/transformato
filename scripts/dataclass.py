from dataclasses import dataclass

#dataclass
@dataclass
class solvation:
    steps_for_equilibration: int
@dataclass
class parameters:
    nstep: int
    nstdcd: int
    nstout: int
    cons: str
    dt: float
@dataclass
class simulation:
    parameters: parameters
    free_energy_type: str
@dataclass
class waterbox:
    dirname: str
    psf_file_name: str
    crd_file_name: str
    rst_file_name: str
    simulation_parameter: str
    intermediate_filename: str
@dataclass
class vacuum:
    dirname: str
    psf_file_name: str
    crd_file_name: str
    rst_file_name: str
    simulation_parameter: str
    intermediate_filename: str
@dataclass
class structure2:
    name: str
    tlc: str
    vacuum: vacuum
    waterbox: waterbox
    charmm_gui_dir: str
@dataclass
class waterbox:
    dirname: str
    psf_file_name: str
    crd_file_name: str
    rst_file_name: str
    simulation_parameter: str
    intermediate_filename: str
@dataclass
class vacuum:
    dirname: str
    psf_file_name: str
    crd_file_name: str
    rst_file_name: str
    simulation_parameter: str
    intermediate_filename: str
@dataclass
class structure1:
    name: str
    tlc: str
    vacuum: vacuum
    waterbox: waterbox
    charmm_gui_dir: str
@dataclass
class system:
    structure1: structure1
    structure2: structure2
    name: str
@dataclass
class input_dataclass:
    system: system
    simulation: simulation
    solvation: solvation
    bin_dir: str
    analysis_dir_base: str
    data_dir_base: str
    system_dir: str
    cluster_dir: str
