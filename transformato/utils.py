import os
import yaml
import parmed as pm

def get_data_dir():
    """Returns the data directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))

def get_config_dir():
    """Returns the data directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'config'))

def get_analysis_dir():
    """Returns the analysis directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'analysis'))

def get_bin_dir():
    """Returns the bin directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'bin'))

def get_toppar_dir():
    """Returns the bin directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'data/toppar'))

def load_config_yaml(system):

    with open(get_config_dir() + '/' + system + '.yaml', 'r') as stream:
        try:
            settingsMap = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    # set the bin, data and analysis dir
    settingsMap['bin_dir'] = get_bin_dir()
    settingsMap['analysis_dir_base'] = get_analysis_dir()
    settingsMap['data_dir_base'] = get_data_dir()
    system_name = settingsMap['system']['structure1']['name'] + '-' + settingsMap['system']['structure2']['name']
    system_path = settingsMap['analysis_dir_base'] + '/' + system_name
    settingsMap['system_dir'] = system_path
    settingsMap['cluster_dir'] = '/data/local/' + system_name

    settingsMap['system']['structure1']['charmm_gui_dir'] =  settingsMap['data_dir_base'] + '/' + settingsMap['system']['structure1']['name'] + '/'
    settingsMap['system']['structure2']['charmm_gui_dir'] =  settingsMap['data_dir_base'] + '/' + settingsMap['system']['structure2']['name'] + '/'
    settingsMap['system']['name'] = system_name
    return settingsMap

