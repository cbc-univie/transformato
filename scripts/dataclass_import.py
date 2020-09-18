import os
import yaml
from dataclasses import dataclass
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
    system_name = f"{settingsMap['system']['structure1']['name']}-{settingsMap['system']['structure2']['name']}-{settingsMap['simulation']['free-energy-type']}"
    settingsMap['system_dir'] = f"{settingsMap['analysis_dir_base']}/{system_name}"
    settingsMap['cluster_dir'] = f"/data/local/{system_name}"

    settingsMap['system']['structure1']['charmm_gui_dir'] = f"{settingsMap['data_dir_base']}/{settingsMap['system']['structure1']['name']}/"
    settingsMap['system']['structure2']['charmm_gui_dir'] = f"{settingsMap['data_dir_base']}/{settingsMap['system']['structure2']['name']}/"
    settingsMap['system']['name'] = system_name

    return settingsMap





def dict_depth(d):
    if isinstance(d, dict):
        return 1 + (max(map(dict_depth, d.values())) if d else 0)
    return 0


def dict_generator(indict, pre=None):   
    pre = pre[:] if pre else []
    if isinstance(indict, dict):
        for key, value in indict.items():
            if isinstance(value, dict):
                for d in dict_generator(value, pre + [key]):
                    yield d
            elif isinstance(value, list) or isinstance(value, tuple):
                for v in value:
                    for d in dict_generator(v, pre + [key]):
                        yield d
            else:
                yield pre + [key, value]
    else:
        yield pre + [indict]


class CodeBlock():

    def __init__(self, head, block):
        self.head = head
        self.block = block
        
    def __str__(self, indent=""):
        result = indent + self.head + ':\n'
        indent += '    '
        for block in self.block:
            if isinstance(block, CodeBlock):
                result += block.__str__(indent)
            else:
                result += indent + block + '\n'
        return result

class create_dataclass_file(object):

    def __init__(self, configuration):
        self.configuration = configuration
        self.generator = dict_generator(self.configuration)
        #self.__gen_test__(self.generator)
        self.tracker = self.__structure__(self.generator)
        print (self.tracker)
        #self.results = self.__dict_to_dataclass__(self.configuration)
        #self.__parser__(self.results)

    def __gen_test__(self, generator):
        typelst = []
        for i in self.generator:
            typelst.append(type(i[-1]))

        print (typelst)

        if typelst[0] is str:
            print ('Yolo')
        
    def __structure__(self, generator):
        tracker = []
        for i in self.generator:
            #print (i)
            depth_level = 0
            for j in i[:-1]:
                #print(j)
                tracker.append((j,depth_level))
                #print (tracker)
                depth_level += 1

        return tracker


    def __dict_to_dataclass__(self, configuration):
        self.setup = 'from dataclasses import dataclass' + '\n' + '\n' + '#dataclass' + '\n'


        self.code = str(CodeBlock('def print_success(x)', [configuration, 'print "Def finished"']))
        self.results = self.setup + self.code
        return self.results

    def __parser__(self, results):

        file_path = os.getcwd()
        file_name = '/dataclass.py'
        tmp_path = file_path + file_name

        try:
            with open(tmp_path, 'w') as f:
                f.write(self.results)
                f.close()

        except IOError:
            print(f"Data class could not be created: {file_name}")
            #logger.info(f"Data class could not be created: {file_name}")
            pass

    


configuration = load_config_yaml(config='config/ethane-ethanol-solvation-free-energy.yaml',
                                   input_dir='.', output_dir='data/')  



def iter_paths(d):
    def iter1(d, path):
        paths = []
        for k, v in d.items():
            if isinstance(v, dict):
                paths += iter1(v, path + [k])
            paths.append((path + [k], v))
        return paths
    return iter1(d, [])

#print(iter_paths(configuration))

#for x in configuration:
    #print(x)
#print (settingsMap)
#configuration= CodeBlock('if x>0', ['print x', 'print "Finished."'])
#configuration = CodeBlock('def print_success(x)', [ifblock, 'print "Def finished"'])
#print (configuration)
create_dataclass_file(configuration)
#print (configuration)
