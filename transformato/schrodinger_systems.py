import numpy as np
import logging
import transformato
logger = logging.getLogger(__name__)


from transformato.mutate import ProposeMutationRoute
from transformato.state import IntermediateStateFactory
from transformato.system import SystemStructure
from transformato.utils import load_config_yaml


def cdk2_structure1(conf_path: str, input_dir: str, output_dir: str):
    configuration = load_config_yaml(
        config=conf_path, input_dir=input_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    # s1_to_s2.calculate_common_core()
    s1_to_s2.propose_common_core()
    connected_dummy_regions_cc1 = [set(s1_to_s2.get_idx_not_in_common_core_for_mol1())]
    connected_dummy_regions_cc2 = [set(s1_to_s2.get_idx_not_in_common_core_for_mol2())]

    s1_to_s2.finish_common_core(
        connected_dummy_regions_cc1=connected_dummy_regions_cc1,
        connected_dummy_regions_cc2=connected_dummy_regions_cc2,
    )

    ###############################
    ####### Structure 1    ########
    ###############################

    # generate the mutation list for the original
    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())
    output_files = []
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    # write endpoint mutation
    output_file_base, intst = i.write_state(mutation_conf=[], intst_nr=1)
    output_files.append(output_file_base)

    # start with charges
    # turn off charges
    charges = mutation_list["charge"]

    nr_of_mutation_steps_charge = (
        1  # defines the number of mutation steps for charge mutation
    nr_of_mutation_steps_charge = (configuration['system']['structure1']['mutation']['steps_charge'])  # defines the number of mutation steps for charge mutation
    
    for lambda_value in np.linspace(1, 0, nr_of_mutation_steps_charge + 1)[1:]:
        output_file_base, intst = i.write_state(
        mutation_conf=charges,
        lambda_value_electrostatic=lambda_value,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # Turn off hydrogens
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base, intst = i.write_state(
    mutation_conf=hydrogen_lj_mutations,
    lambda_value_vdw=0.0,
    intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off lj of heavy atoms
    if type(configuration['system']['structure1']['mutation']['heavy_atoms']) == list:
        for h in configuration['system']['structure1']['mutation']['heavy_atoms']:
            logger.info(f'turning off lj of heavy atom: {h}')
            d = transformato.utils.map_lj_mutations_to_atom_idx(mutation_list["lj"])
            for m in [
                [d[(h,)]],
            ]:
                output_file_base, intst = i.write_state(
                    mutation_conf=m,
                    lambda_value_vdw=0.0,
                    intst_nr=intst,
                )
                output_files.append(output_file_base)
    else:
        logger.info(f'there are no heavy atoms indicatet in the yaml file')

    # generate terminal lj
    output_file_base, intst = i.write_state(
    mutation_conf=mutation_list["default-lj"],
    lambda_value_vdw=0.0,
    intst_nr=intst,

    )
    output_files.append(output_file_base)

    # change bonded parameters on common core
    m = mutation_list["transform"]
    for lambda_value in np.linspace(0.75, 0, configuration['system']['structure1']['mutation']['steps_common_core']):
    # interpolate between parameters
        output_file_base, intst = i.write_state(
        mutation_conf=m,
        common_core_transformation=lambda_value,
        intst_nr=intst,
    )
    output_files.append(output_file_base)
    return output_files, configuration # makes further processing easier

    #system2 is not necessary since we should use the same structure 2 for all cdk2 simulations


def cdk2_both_structures(conf_path: str, input_dir: str, output_dir: str):
    configuration = load_config_yaml(
        config=conf_path, input_dir=input_dir, output_dir=output_dir
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    #s1_to_s2.calculate_common_core()
    s1_to_s2.propose_common_core()
    connected_dummy_regions_cc1=[set(s1_to_s2.get_idx_not_in_common_core_for_mol1())]
    connected_dummy_regions_cc2=[set(s1_to_s2.get_idx_not_in_common_core_for_mol2())]

    s1_to_s2.finish_common_core(connected_dummy_regions_cc1=connected_dummy_regions_cc1, connected_dummy_regions_cc2=connected_dummy_regions_cc2)

    ###############################
    ####### Structure 1    ########
    ###############################
    
    # generate the mutation list for the original
    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())
    output_files = []
    i = IntermediateStateFactory(
    system=s1,
    configuration=configuration,
    )

    # write endpoint mutation
    output_file_base, intst = i.write_state(mutation_conf=[], intst_nr=1)
    output_files.append(output_file_base)

    # start with charges
    # turn off charges
    charges = mutation_list["charge"]
    nr_of_mutation_steps_charge = (configuration['system']['structure1']['mutation']['steps_charge'])  # defines the number of mutation steps for charge mutation
    
    for lambda_value in np.linspace(1, 0, nr_of_mutation_steps_charge + 1)[1:]:
        output_file_base, intst = i.write_state(
            mutation_conf=charges,
            lambda_value_electrostatic=lambda_value,
            intst_nr=intst,
        )
    output_files.append(output_file_base)

    # Turn off hydrogens
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base, intst = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # turn off lj of heavy atoms
    if type(configuration['system']['structure1']['mutation']['heavy_atoms']) == list:
        for h in configuration['system']['structure1']['mutation']['heavy_atoms']:
            logger.info(f'turning off lj of heavy atom: {h}')
            d = transformato.utils.map_lj_mutations_to_atom_idx(mutation_list["lj"])
            for m in [
                [d[(h,)]],
            ]:
                output_file_base, intst = i.write_state(
                    mutation_conf=m,
                    lambda_value_vdw=0.0,
                    intst_nr=intst,
                )
                output_files.append(output_file_base)
    else:
        logger.info(f'there are no heavy atoms indicatet in the yaml file')


    # generate terminal lj
    output_file_base, intst = i.write_state(
        mutation_conf=mutation_list["default-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)

    # change bonded parameters on common core
    m = mutation_list["transform"]

    for lambda_value in np.linspace(0.75, 0, 2):
        # interpolate between parameters
    for lambda_value in np.linspace(0.75, 0, configuration['system']['structure1']['mutation']['steps_common_core']):
    # interpolate between parameters

        output_file_base, intst = i.write_state(
            mutation_conf=m,
            common_core_transformation=lambda_value,
            intst_nr=intst,
        )
    output_files.append(output_file_base)

    ###############################
    ###### Structure 2     ########
    ###############################

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()

    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    output_files = []
    output_file_base, intst = i.write_state(mutation_conf=[], intst_nr=1)
    output_files.append(output_file_base)

    # start with charges
    # turn off charges
    charges = mutation_list["charge"]

    nr_of_mutation_steps_charge = (
        1  # defines the number of mutation steps for charge mutation
    )
    nr_of_mutation_steps_charge = (configuration['system']['structure2']['mutation']['steps_charge'])  # defines the number of mutation steps for charge mutation

    for lambda_value in np.linspace(1, 0, nr_of_mutation_steps_charge + 1)[1:]:
        output_file_base, intst = i.write_state(
            mutation_conf=charges,
            lambda_value_electrostatic=lambda_value,
            intst_nr=intst,
        )
    output_files.append(output_file_base)

    # Turn off hydrogens
    hydrogen_lj_mutations = mutation_list["hydrogen-lj"]
    output_file_base, intst = i.write_state(
        mutation_conf=hydrogen_lj_mutations,
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )
    output_files.append(output_file_base)
    d = transformato.utils.map_lj_mutations_to_atom_idx(mutation_list["lj"])

    # turn off lj of heavy atoms
    if type(configuration['system']['structure2']['mutation']['heavy_atoms']) == list:
        for h in configuration['system']['structure2']['mutation']['heavy_atoms']:
            logger.info(f'turning off lj of heavy atom: {h}')
            d = transformato.utils.map_lj_mutations_to_atom_idx(mutation_list["lj"])
            for m in [
                [d[(h,)]],
            ]:
                output_file_base, intst = i.write_state(
                    mutation_conf=m,
                    lambda_value_vdw=0.0,
                    intst_nr=intst,
                )
                output_files.append(output_file_base)
    else:
        logger.info(f'there are no heavy atoms indicatet in the yaml file')

    # generate terminal lj
    output_file_base, intst = i.write_state(
        mutation_conf=mutation_list["default-lj"],
        lambda_value_vdw=0.0,
        intst_nr=intst,
    )

    
