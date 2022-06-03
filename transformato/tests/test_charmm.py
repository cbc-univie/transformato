import os
import pytest
import warnings
warnings.filterwarnings("ignore", module="parmed")


@pytest.mark.rbfe
@pytest.mark.charmm
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_28_1h1q_rbfe_production_with_CHARMM():
    # Generating output for a run of the CDK2 Ligand System
    from transformato import (
        IntermediateStateFactory,
        ProposeMutationRoute,
        SystemStructure,
        load_config_yaml,
    )

    from transformato.mutate import perform_mutations
    from transformato.utils import run_simulation

    configuration = load_config_yaml(
        config="transformato/tests/config/test-28_1h1q_rbfe.yaml",
        input_dir="data/",
        output_dir="/site/raid4/johannes/test",
    )

    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)

    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    perform_mutations(
        configuration=configuration,
        i=i,
        mutation_list=mutation_list,
        nr_of_mutation_steps_charge=1,
        nr_of_mutation_steps_cc=1,
    )
    # run_simulation(i.output_files, engine="charmm")


@pytest.mark.rbfe
@pytest.mark.charmm
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_28_1h1q_rsfe_analysis_with_CHARMM():
    # Running the analysis, first the productions script needs to be finisehd
    # Change output dir
    # Check transformato/bin/charmm_eval_energy.sh for the correct name of the charmm
    # executable!
    from transformato import load_config_yaml
    from transformato.utils import postprocessing
    from transformato.utils import postprocessing

    configuration = load_config_yaml(
        config="transformato/tests/config/test-28_1h1q_rsfe.yaml",
        input_dir="data/",
        output_dir="/site/raid4/johannes/test",
    )

    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration,
        name="cdk2-28",
        engine="CHARMM",
        max_snapshots=10,
        num_proc=1,
        analyze_traj_with="mdtraj",
    )
    print(f"Free energy difference: {ddG_openMM} +- {dddG} [kT")


@pytest.mark.rsfe
@pytest.mark.charmm
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_run_1a0q_1a07_rsfe_production_with_CHARMM(caplog):

    from transformato import (
        IntermediateStateFactory,
        ProposeMutationRoute,
        SystemStructure,
        load_config_yaml,
    )
    from transformato.mutate import perform_mutations

    workdir = "/site/raid4/johannes/test"
    conf = "transformato/tests/config/test-1a0q-1a07-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/test_systems_mutation", output_dir=workdir
    )
    s1 = SystemStructure(configuration, "structure1")
    s2 = SystemStructure(configuration, "structure2")
    s1_to_s2 = ProposeMutationRoute(s1, s2)
    s1_to_s2.propose_common_core()
    s1_to_s2.finish_common_core()

    # generate the mutation list for the original
    mutation_list = s1_to_s2.generate_mutations_to_common_core_for_mol1()
    print(mutation_list.keys())
    i = IntermediateStateFactory(
        system=s1,
        configuration=configuration,
    )
    perform_mutations(
        configuration=configuration,
        nr_of_mutation_steps_cc=2,
        nr_of_mutation_steps_charge=1,
        i=i,
        mutation_list=mutation_list,
    )
    # run_simulation(i.output_files, engine="CHARMM")

    # ddG_openMM, dddG = postprocessing(
    #     configuration,
    #     name="1a0q",
    #     engine="openMM",
    #     max_snapshots=5000,
    # )
    # print(ddG_openMM, dddG)


@pytest.mark.rsfe
@pytest.mark.charmm
@pytest.mark.skipif(
    os.getenv("CI") == "true",
    reason="Skipping tests that cannot pass in github actions",
)
def test_analyse_1a0q_1a07_rsfe_with_openMM(caplog):
    # Running the analysis, first the productions script needs to be finisehd
    # Change output dir
    # Check transformato/bin/charmm_eval_energy.sh for the correct name of the charmm
    # executable!
    from transformato import load_config_yaml
    from transformato.utils import postprocessing
    from transformato.utils import postprocessing

    workdir = "/site/raid4/johannes/test"
    conf = "transformato/tests/config/test-1a0q-1a07-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="data/test_systems_mutation", output_dir=workdir
    )

    ddG_openMM, dddG = postprocessing(
        configuration,
        name="1a0q",
        engine="CHARMM",
        max_snapshots=10,
        num_proc=1,
        analyze_traj_with="mdtraj",
    )

    print(ddG_openMM, dddG)
