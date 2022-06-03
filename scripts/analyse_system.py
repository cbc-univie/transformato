from transformato import load_config_yaml, FreeEnergyCalculator


def postprocessing(
    configuration: dict,
    name: str = "methane",
    engine: str = "openMM",
    max_snapshots: int = 300,
    show_summary: bool = False,
):
    from transformato import FreeEnergyCalculator

    f = FreeEnergyCalculator(configuration, name)
    f.load_trajs(nr_of_max_snapshots=max_snapshots)

    f.calculate_dG_to_common_core(engine=engine)
    ddG, dddG = f.end_state_free_energy_difference
    print("######################################")
    print("Free energy to common core in kT")
    print("######################################")
    print(f"Free energy difference: {ddG} [kT]")
    print(f"Uncertanty: {dddG} [kT]")
    print("######################################")
    print("######################################")
    if show_summary:
        f.show_summary()
    return ddG, dddG, f


def calculate_rbfe():

    system_name = "1oi9"
    conf = f"/data/shared/projects/rbfe-tf/config/cdk2/{system_name}.yaml"
    configuration = load_config_yaml(
        config=conf,
        input_dir="../../data/",
        output_dir="/data/shared/projects/rbfe-tf/cdk2/run1/",
    )  # NOTE: for preprocessing input_dir is the output dir

    # 2OJ9-original to tautomer common core
    ddG_openMM, dddG, f_openMM = postprocessing(
        configuration, name="1oi9", engine="openMM", max_snapshots=2_000
    )


if __name__ == "__main__":
    calculate_rbfe()
