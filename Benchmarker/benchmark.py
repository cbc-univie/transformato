from utils import loader,load_param_yaml, createList, createMap, createStepsMap, targetList, create_input_yaml, check_folder
from simulation import benchmark_simulation
import os

user_input = input("Enter the path to destination folder: ")
#input = '/home/master/Master/benchmarker/test/input'
loader(user_input)

path = user_input + '/transformato-systems-master/'

parameter = load_param_yaml(parameter=path + 'parameter.yaml') 

steps_map = createStepsMap(parameter['simulation']['parameters']['nstep'],parameter['simulation']['parameters']['nstdcd'])

molecule_pairs = targetList(path)

create_input_yaml(path,molecule_pairs)

systems = next(os.walk(path + 'output'))[1]

for i in systems:
    benchmark_simulation(steps_map,path,i)








