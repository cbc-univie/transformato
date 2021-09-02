import numpy as np
from transformato.constants import initialize_NUM_PROC

# read in specific topology with parameters
from transformato import load_config_yaml, FreeEnergyCalculator


def postprocessing(
    configuration: dict,
    name: str = "methane",
    engine: str = "openMM",
    max_snapshots: int = 300,
    show_summary: bool = False,
    different_path_for_dcd: str = "",
    only_vacuum: bool = False,
):
    f = FreeEnergyCalculator(configuration, name)
    if only_vacuum:
        f.envs = ["vacuum"]

    if different_path_for_dcd:
        # this is needed if the trajectories are stored at a different location than the
        # potential definitions
        path = f.base_path
        f.base_path = different_path_for_dcd
        f.load_trajs(nr_of_max_snapshots=max_snapshots)
        f.base_path = path
    else:
        f.load_trajs(nr_of_max_snapshots=max_snapshots)

    f.calculate_dG_to_common_core(engine=engine)
    if only_vacuum:
        return -1, -1, f
    else:
        ddG, dddG = f.end_state_free_energy_difference
        print(f"Free energy difference: {ddG}")
        print(f"Uncertanty: {dddG}")
        if show_summary:
            f.show_summary()
        return ddG, dddG, f


def calculate_rsfe_mp():

    initialize_NUM_PROC(6)
    conf = "config/test-2oj9-tautomer-pair-rsfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="../../data/", output_dir="../../data"
    )  # NOTE: for preprocessing input_dir is the output dir

    # 2OJ9-original to tautomer common core
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="2OJ9-original", engine="openMM", max_snapshots=800
    )


def calculate_rbfe_mp():

    initialize_NUM_PROC(6)
    conf = "config/test-2oj9-tautomer-pair-rbfe.yaml"
    configuration = load_config_yaml(
        config=conf, input_dir="../../data/", output_dir="../../data"
    )  # NOTE: for preprocessing input_dir is the output dir

    # 2OJ9-original to tautomer common core
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="2OJ9-original", engine="openMM", max_snapshots=1000
    )


if __name__ == "__main__":
    # calculate_rsfe_mp()
    calculate_rbfe_mp()
