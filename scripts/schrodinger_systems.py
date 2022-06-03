import numpy as np
import logging
import transformato

logger = logging.getLogger(__name__)
from transformato.mutate import ProposeMutationRoute, perform_mutations
from transformato.state import IntermediateStateFactory
from transformato.system import SystemStructure
from transformato.utils import load_config_yaml


def cdk2_single_structure(conf_path: str, input_dir: str, output_dir: str):
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
    # generate the mutation list for structure 1
    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())

    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    output_files = perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
    )

    return output_files, configuration
