from numpy import arange 
import yaml
from urllib.request import urlretrieve
import os
from zipfile import ZipFile


def loader(input_location):
    try:
        os.makedirs(input_location)
        print ('Directory ' , input_location ,  ' Created')
    except FileExistsError:
        print('Directory ', input_location , ' already exists')
        pass
    
    print ("Downloading benchmark molecules...")
    url = 'https://github.com/wiederm/transformato-systems/archive/master.zip'
    dest_file = os.path.join(input_location, 'transformato-systems-master.zip')
    urlretrieve(url, dest_file)
    
    print ("Extracting...")
    with ZipFile(dest_file, 'r') as zipObj:
        zipObj.extractall(input_location)

    targets = next(os.walk(input_location + '/transformato-systems-master'))[1]

    try:
        os.makedirs(input_location + '/transformato-systems-master/output')
    except FileExistsError:
        pass
    
    print ("Cleaning...")
    os.remove(dest_file)
    print ("Finished!")
    return targets

def createList(r1,r2,d):
    return arange(r1,r2+1,d)

def createMap(nsteps_list,nstdcd_list):
    return [(x,y) for x in nsteps_list for y in nstdcd_list]

def createStepsMap(param1,param2):
    nsteps_list = createList(param1[0],param1[1],param1[2])
    nstdcd_list = createList(param2[0],param2[1],param2[2])
    steps_map = createMap(nsteps_list,nstdcd_list)
    return steps_map

def load_param_yaml(parameter):

    with open(f"{parameter}", 'r') as stream:
        try:
            settingsMap = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    if settingsMap['simulation']['parameters'].get('nstep') == None or settingsMap['simulation']['parameters'].get('nstdcd') == None:
        raise KeyError('nsteps or nstdcd is not defined in parameter file')

    return settingsMap

def check_folder(path):
    dir_list =  [x[0] for x in os.walk(path)]

    for dir in dir_list:
        if '/lig' in dir:
            tlc = 'LIG'
        elif '/unk' in dir:
            tlc = 'UNK'
        elif '/unl' in dir:
            tlc = 'UNL'

    return tlc


def create_input_yaml(path,pairs):

    for pair in pairs:

        new_path = path + 'output/' + pair[0] + '-' + pair[1] + '-solvation-free-energy'

        try:
            os.makedirs(new_path)       
        except FileExistsError:       
            pass

        tlc1 = check_folder(path + pair[0])
        tlc2 = check_folder(path + pair[1])

        input_dict = {'system': 
        {'structure1': 
        {'name': pair[0], 
        'tlc': tlc1, 
        'vacuum': 
        {'dirname': 'vacuum', 
        'psf_file_name': 'step3_charmm2omm', 
        'crd_file_name': 'step3_charmm2omm', 
        'rst_file_name': 'step4_equilibration',
        'simulation_parameter': 'step5_production.inp', 
        'intermediate-filename': 'lig_in_vacuum'}, 
        'waterbox':
        {'dirname': 'waterbox', 
        'psf_file_name': 'step3_charmm2omm', 
        'crd_file_name': 'step3_charmm2omm', 
        'rst_file_name': 'step4_equilibration', 
        'simulation_parameter': 'step5_production.inp', 
        'intermediate-filename': 'lig_in_waterbox'}}, 
        'structure2': 
        {'name': pair[1], 
        'tlc': tlc2, 
        'vacuum': 
        {'dirname': 'vacuum', 
        'psf_file_name': 'step3_charmm2omm', 
        'crd_file_name': 'step3_charmm2omm', 
        'rst_file_name': 'step4_equilibration', 
        'simulation_parameter': 'step5_production.inp', 
        'intermediate-filename': 'lig_in_vacuum'}, 
        'waterbox': 
        {'dirname': 'waterbox', 
        'psf_file_name': 'step3_charmm2omm', 
        'crd_file_name': 'step3_charmm2omm', 
        'rst_file_name': 'step4_equilibration', 
        'simulation_parameter': 'step5_production.inp', 
        'intermediate-filename': 'lig_in_waterbox'}}}, 
        'simulation': 
        {'parameters': 
        {'nstep': 500000, 
        'nstdcd': 10000, 
        'nstout': 10000, 
        'cons': 'None', 
        'dt': 0.001}, 
        'free-energy-type': 'solvation-free-energy'}, 
        'solvation': 
        {'steps_for_equilibration': 1000}}

        with open(f"{new_path + '/input.yaml'}", 'w') as file:
            try:
                yaml.dump(input_dict, file)
            except yaml.YAMLError as exc:
                print(exc)

def reverse(tuples): 
    new_tup = tuples[::-1] 
    return new_tup 


def targetList(input_location):
    targets = next(os.walk(input_location))[1]

    try:
        targets.remove('output')
    except ValueError:
        pass 
    except AttributeError:
        print("targets is not a list!")

    targets = createMap(targets,targets)

    for pair in targets:
        if pair[0] == pair[1]:    
            targets.remove(pair)

    for pair in targets:
        for other_pair in targets:
            if pair == reverse(other_pair):
                targets.remove(other_pair)

    return targets

