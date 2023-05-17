# Import package, test suite, and other packages as needed
import logging
import os
import warnings
from io import StringIO
import pdb
import parmed as pm
import pytest
import pickle

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


@pytest.mark.point_mutation
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_point_mutation_2():
    configuration = load_config_yaml(
        config=f"/site/raid3/johannes/h2u/data/config/H2U_1_cano-H2U_1_mod.yaml",
        input_dir="/site/raid3/johannes/h2u/data",
        output_dir=f"/site/raid3/johannes",
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
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


@pytest.mark.point_mutation
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_setting_up_point_mutation():
    configuration = load_config_yaml(
        config=f"/site/raid3/johannes/bioinfo/data/config/cano10-psu10.yaml",
        input_dir="/site/raid3/johannes/bioinfo/data/psul",
        output_dir="/site/raid3/johannes/bioinfo/test",
    )

    # pdb.set_trace()
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
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


@pytest.mark.point_mutation
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_point_mutation():
    print("Setting up the system, should appear only once")
    configuration = load_config_yaml(
        config=f"/site/raid3/johannes/bioinfo/data/config/cano5-psu5.yaml",
        input_dir="/site/raid3/johannes/bioinfo/data/psu",
        output_dir="/site/raid3/johannes/bioinfo/test3",
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
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

    # mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol2()
    # i = IntermediateStateFactory(
    #     system=s2,
    #     configuration=configuration,
    # )

    # perform_mutations(
    #     configuration=configuration,
    #     nr_of_mutation_steps_charge=3,
    #     i=i,
    #     mutation_list=mutation_list,
    # )


# def expensive_computation():
#     print("called expensive_computation")
#     return "Hallo","wie","gehts"

# @pytest.fixture()
# def test_create(request):
#     print("Called test_create")
#     test = request.config.cache.get("test", None)
#     if test is None:
#         test = expensive_computation()
#         request.config.cache.set("test", test)
#     return test

# def test_use(test_create):
#     print("Run test_use")
#     s1,s2,s3 = test_create
#     print(s1,s2,s3)

# def test_use2(test_create):
#     print("Run test_use2")
#     assert test_create == "Hallo"
