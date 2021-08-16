import numpy as np

import transformato
from transformato.mutate import ProposeMutationRoute
from transformato.state import IntermediateStateFactory
from transformato.system import SystemStructure
from transformato.utils import load_config_yaml
from transformato.constants import check_platform


def cdk2_withoutHeavyAtom(conf_path: str, input_dir: str, output_dir: str):
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
    nr_of_mutation_steps_charge = (
    1  # defines the number of mutation steps for charge mutation
    )
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

    # generate terminal lj
    output_file_base, intst = i.write_state(
    mutation_conf=mutation_list["default-lj"],
    lambda_value_vdw=0.0,
    intst_nr=intst,
    )
    output_files.append(output_file_base)