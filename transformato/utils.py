import os
import logging

import parmed as pm
import yaml

from io import StringIO

logger = logging.getLogger(__name__)


def get_bin_dir():
    """Returns the bin directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "bin"))


def get_toppar_dir():
    """Returns the bin directory of this package"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "toppar"))


def map_lj_mutations_to_atom_idx(lj_mutations: list) -> dict:
    """
    map_lj_mutations_to_atom_idx
    returns a mapping between the mutations and the atom idx

    Parameters
    ----------
    lj_mutations : list of lj mutation objects
        [description]

    Returns
    -------
    dict
    """
    d = {}
    for lj in lj_mutations:
        key = lj.vdw_atom_idx
        d[tuple(key)] = lj

    return d


def print_mutations(mutation: dict):
    # TODO: this needs to be finished!
    if "charge" in mutation.keys():
        for charge_mutation in mutation["charge"]:
            print("Charge mutation")
            print(f"Atoms to be mutated: {charge_mutation.atoms_to_be_mutated}")
    if "hydrogen-lj" in mutation.keys():
        print("Hydrogen LJ mutation")
        print(f"Atoms to be mutated: {charge_mutation.atoms_to_be_mutated}")
        print(f"LJ mutation to default:")

    if "transform" in mutation.keys():
        pass


def load_config_yaml(config, input_dir, output_dir) -> dict:

    with open(f"{config}", "r") as stream:
        try:
            settingsMap = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    if (
        settingsMap["simulation"]["parameters"].get("nstep") == None
        or settingsMap["simulation"]["parameters"].get("nstdcd") == None
    ):
        raise KeyError("nsteps or nstdcd is not defined in config file")
    else:
        if (
            settingsMap["simulation"]["parameters"]["nstep"]
            / settingsMap["simulation"]["parameters"]["nstdcd"]
            < 20
        ):
            raise RuntimeError(
                "nsteps size and nstdcd size in config file does not match"
            )
    # making tlc caps always uppercase
    settingsMap["system"]["structure1"]["tlc"] = settingsMap["system"]["structure1"][
        "tlc"
    ].upper()
    settingsMap["system"]["structure2"]["tlc"] = settingsMap["system"]["structure2"][
        "tlc"
    ].upper()
    # set the bin, data and analysis dir
    settingsMap["bin_dir"] = get_bin_dir()
    settingsMap["analysis_dir_base"] = os.path.abspath(f"{output_dir}")
    settingsMap["data_dir_base"] = os.path.abspath(f"{input_dir}")
    system_name = f"{settingsMap['system']['structure1']['name']}-{settingsMap['system']['structure2']['name']}-{settingsMap['simulation']['free-energy-type']}"
    settingsMap["system_dir"] = f"{settingsMap['analysis_dir_base']}/{system_name}"
    settingsMap["cluster_dir"] = f"/data/local/{system_name}"

    settingsMap["system"]["structure1"][
        "charmm_gui_dir"
    ] = f"{settingsMap['data_dir_base']}/{settingsMap['system']['structure1']['name']}/"
    settingsMap["system"]["structure2"][
        "charmm_gui_dir"
    ] = f"{settingsMap['data_dir_base']}/{settingsMap['system']['structure2']['name']}/"
    settingsMap["system"]["name"] = system_name
    return settingsMap


def get_structure_name(configuration: dict, structure_name: str):
    if configuration["system"]["structure1"]["name"] == structure_name:
        return "structure1"
    elif configuration["system"]["structure2"]["name"] == structure_name:
        return "structure2"
    else:
        raise RuntimeError(f"Could not finde structure entry for : {structure_name}")


def psf_correction(str_object: StringIO):
    """Correcting the issue with 2 missing spaces in the waterbox psf files and replacing the !NGRP statement so its possible that CHARMM, CHARMM_OpenMM and OpenMM can handle the correspoing psf file"""
    str_object = str_object.getvalue()  # get the values as a long string
    new_str = ""
    i = 0
    second_line = -1
    correction_on = False
    for line in str_object.split("\n"):  # split on newline charactar
        i += 1
        if "!NATOM" in line:  # if !NATOM is found start correction mode
            new_str += f"{line}\n"
            correction_on = True
            continue

        if "!NBOND" in line:  # if !NBOND is found exit correction mode
            correction_on = False

        if (
            correction_on == True
        ):  # if in correction mode take the string, split on whitespace and put the values in a newly formated string
            values = line.split()
            if len(values) == 11:  # skip empty lines
                new_str += f"{values[0]:>10} {values[1]:8} {values[2]:8} {values[3]:8} {values[4]:8} {values[5]:6} {values[6]:>10}{values[7]:>14}{values[8]:>12}{values[9]:>10}{values[10]:>18}\n"
            elif len(values) == 8:
                values.extend(["0", "0.00000", "0.00000000000"])
                new_str += f"{values[0]:>10} {values[1]:8} {values[2]:8} {values[3]:8} {values[4]:8} {values[5]:6} {values[6]:>10}{values[7]:>14}{values[8]:>12}{values[9]:>10}{values[10]:>18}\n"
            elif len(values) == 9:
                values.extend(["0.00000", "0.00000000000"])
                new_str += f"{values[0]:>10} {values[1]:8} {values[2]:8} {values[3]:8} {values[4]:8} {values[5]:6} {values[6]:>10}{values[7]:>14}{values[8]:>12}{values[9]:>10}{values[10]:>18}\n"

            elif len(values) == 0:
                new_str += f"{line}\n"
            else:
                raise RuntimeError(f"Error with the psf file: {line}")

        else:  # otherwise add line to new_str

            if "!NGRP NST2" in line:
                second_line = i+1 # we want to remove the next line after !NGRP appears
                new_str += f"{line.replace('1','0')}\n"
            elif i == second_line:
                new_str += ' \n'
            else:
                new_str += f"{line}\n"



    return new_str
