import os
from dataclasses import dataclass
import parmed as pm
import yaml
from dacite import from_dict


def get_bin_dir():
    """Returns the bin directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'bin'))


def get_toppar_dir():
    """Returns the toppar directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'toppar'))


def load_config_yaml(config, input_dir, output_dir):

    with open(f"{config}", 'r') as stream:
        try:
            settingsMap = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    if settingsMap['simulation']['parameters'].get('nstep') == None or settingsMap['simulation']['parameters'].get('nstdcd') == None:
        raise KeyError('nsteps or nstdcd is not defined in config file')
    else:
        if settingsMap['simulation']['parameters']['nstep']/settingsMap['simulation']['parameters']['nstdcd'] < 20:
            raise RuntimeError('nsteps size and nstdcd size in config file does not match')

    # set the bin, data and analysis dir
    settingsMap['bin_dir'] = get_bin_dir()
    settingsMap['analysis_dir_base'] = os.path.abspath(f"{output_dir}")
    settingsMap['data_dir_base'] = os.path.abspath(f"{input_dir}")
    system_name = f"{settingsMap['system']['structure1']['name']}-{settingsMap['system']['structure2']['name']}-{settingsMap['simulation']['free_energy_type']}"
    settingsMap['system_dir'] = f"{settingsMap['analysis_dir_base']}/{system_name}"
    settingsMap['cluster_dir'] = f"/data/local/{system_name}"

    settingsMap['system']['structure1']['charmm_gui_dir'] = f"{settingsMap['data_dir_base']}/{settingsMap['system']['structure1']['name']}/"
    settingsMap['system']['structure2']['charmm_gui_dir'] = f"{settingsMap['data_dir_base']}/{settingsMap['system']['structure2']['name']}/"
    settingsMap['system']['name'] = system_name
    return settingsMap

def fill_dataclass(input_dataclass,configuration):
    filled_class = from_dict(data_class=input_dataclass, data=configuration)
    return filled_class

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
class environment:
    dirname: str
    psf_file_name: str
    crd_file_name: str
    rst_file_name: str
    simulation_parameter: str
    intermediate_filename: str
@dataclass
class structure:
    name: str
    tlc: str
    vacuum: environment
    waterbox: environment
    charmm_gui_dir: str   
@dataclass
class system:
    structure1: structure
    structure2: structure
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

configuration = load_config_yaml(config='config/test-2oj9-solvation-free-energy_dataclass.yaml',
                                   input_dir='.', output_dir='data/')

filled_class = fill_dataclass(input_dataclass,configuration)

print(filled_class)