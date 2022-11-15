# Import package, test suite, and other packages as needed
import logging
import os
import warnings
from io import StringIO
import pdb
import parmed as pm
import pytest

# read in specific topology with parameters
# read in specific topology with parameters
from transformato import (
    SystemStructure,
    IntermediateStateFactory,
    ProposeMutationRoute,
    load_config_yaml,
    psf_correction,
)
from transformato.mutate import perform_mutations
from transformato.tests.paths import get_test_output_dir
from transformato_testsystems.testsystems import get_testsystems_dir

warnings.filterwarnings("ignore", module="parmed")


# @pytest.mark.rbfe
# @pytest.mark.requires_parmed_supporting_lp
# @pytest.mark.skipif(
#     os.getenv("CI") == "true",
#     reason="Skipping tests that cannot pass in github actions",
# )
def test_setting_up_point_mutation():

    configuration = load_config_yaml(
        config=f"/site/raid3/johannes/bioinfo/cano2_modif2.yaml",
        input_dir="/site/raid3/johannes/bioinfo",
        output_dir="/site/raid3/johannes/bioinfo",
    )

    # pdb.set_trace()
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(f"Die mutation list {mutation_list}")
    print(mutation_list["transform"][0].__dict__)
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )

    perform_mutations(
        configuration=configuration,
        nr_of_mutation_steps_charge=3,
        nr_of_mutation_steps_cc=3,
        i=i,
        mutation_list=mutation_list,
    )

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    print(f"Die mutation list {mutation_list}")
    i = IntermediateStateFactory(
        system=s2,
        configuration=configuration,
    )

    perform_mutations(
        configuration=configuration,
        nr_of_mutation_steps_charge=3,
        i=i,
        mutation_list=mutation_list,
    )
