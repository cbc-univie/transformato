import os
import yaml
import parmed as pm

def get_bin_dir():
    """Returns the bin directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'bin'))

def get_toppar_dir():
    """Returns the bin directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'toppar'))

def load_config_yaml(config, data_dir, output_dir):

    with open(f"{config}", 'r') as stream:
        try:
            settingsMap = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    # set the bin, data and analysis dir
    settingsMap['bin_dir'] = get_bin_dir()
    settingsMap['analysis_dir_base'] = f"{output_dir}"
    settingsMap['data_dir_base'] = f"{data_dir}"
    system_name = settingsMap['system']['structure1']['name'] + '-' + settingsMap['system']['structure2']['name']
    #system_path = settingsMap['analysis_dir_base'] + '/' + system_name
    #settingsMap['system_dir'] = system_path
    settingsMap['cluster_dir'] = '/data/local/' + system_name

    settingsMap['system']['structure1']['charmm_gui_dir'] =  settingsMap['data_dir_base'] + '/' + settingsMap['system']['structure1']['name'] + '/'
    settingsMap['system']['structure2']['charmm_gui_dir'] =  settingsMap['data_dir_base'] + '/' + settingsMap['system']['structure2']['name'] + '/'
    settingsMap['system']['name'] = system_name
    return settingsMap

