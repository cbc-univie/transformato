# for two test systes dcd files are deposited --- this skript
# regenerates the test data needed
from transformato.tests.test_mutation import (
    setup_acetylacetone_tautomer_pair,
    setup_2OJ9_tautomer_pair_rsfe,
    setup_2OJ9_tautomer_pair_rbfe,
)
from transformato.tests.test_run_production import run_simulation
from transformato import load_config_yaml


def run_2OJ9_tautomer_pair_rsfe():
    conf_path = "config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf_path, input_dir="../../data/", output_dir=get_test_output_dir()
    )

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rsfe(
        configuration=configuration
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)


def run_acetylacetone_tautomer_pair_rsfe():
    conf_path = "config/test-acetylacetone-tautomer-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf_path, input_dir="../../data/", output_dir=get_test_output_dir()
    )

    (output_files_t1, output_files_t2), _, _ = setup_acetylacetone_tautomer_pair(
        configuration=configuration
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)


def run_2OJ9_tautomer_pair_rbfe():
    conf_path = "config/test-2oj9-tautomer-pair-rbfe.yaml"
    configuration = load_config_yaml(
        config=conf_path, input_dir="../../data/", output_dir=get_test_output_dir()
    )

    (output_files_t1, output_files_t2), _, _ = setup_2OJ9_tautomer_pair_rbfe(
        configuration=configuration
    )
    run_simulation(output_files_t1)
    run_simulation(output_files_t2)


if __name__ == "__main__":
    # 2OJ9_tautomer_pair_rsfe()
    # test_run_acetylacetone_tautomer_pair_rsfe()
    run_2OJ9_tautomer_pair_rbfe()
