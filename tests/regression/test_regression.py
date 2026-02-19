import json
import os

import numpy as np
import pytest

from oqd_heisenberg_ion.common.driver.factory import DriverFactory
from oqd_heisenberg_ion.common.inputs.input_reader import InputReader
from oqd_heisenberg_ion.common.preprocessor.factory import PreprocessorFactory

with open(os.path.abspath("tests/regression/regression_results.json"), "r") as f:
    regression_data = json.load(f)


@pytest.mark.parametrize(
    ["hamiltonian_name", "alpha", "loop_type", "case_iter"],
    [
        ["fm_heisenberg_afm_Z", 1.0, "deterministic", 0],
        ["fm_heisenberg_fm_Z", 2.0, "deterministic", 1],
        ["XXZ", 1.5, "directed_loop", 2],
        ["XY", 1.0, "deterministic", 3],
        ["XY", 3.0, "deterministic", 4],
        ["XY", 10.0, "deterministic", 5],
    ],
)
def test_long_range(hamiltonian_name, alpha, loop_type, case_iter, tmp_path):

    # Long range QMC
    input_file = "tests/input_files/long_range.txt"
    inputs = InputReader(input_file_path=input_file)
    inputs.read_inputs_from_file()
    inputs.read_kwarg_inputs(hamiltonian_name=hamiltonian_name, alpha=alpha, loop_type=loop_type, root_folder=tmp_path)
    parameter_set_list = inputs.parameter_set_list
    preprocessor = PreprocessorFactory.create("long_range_qmc", parameter_set_list)
    driver_inputs = preprocessor.preprocess()
    driver = DriverFactory.create("long_range_qmc", preprocessor.simulation_folder, driver_inputs)
    driver.simulate()

    qmc_outputs_file = os.path.join(
        preprocessor.processed_configs[0]["misc"]["run_folder"], "qmc_output/estimators.csv"
    )
    estimator_outputs = compute_qmc_stats(qmc_outputs_file)

    regression_results = regression_data[case_iter]

    assert estimator_outputs[0][1] == pytest.approx(regression_results["energy_mean"])
    assert estimator_outputs[1][1] == pytest.approx(regression_results["energy_std"])
    assert estimator_outputs[0][2] == pytest.approx(regression_results["magnetization_mean"])
    assert estimator_outputs[1][2] == pytest.approx(regression_results["magnetization_std"])
    assert estimator_outputs[0][3] == pytest.approx(regression_results["stiffness_mean"])
    assert estimator_outputs[1][3] == pytest.approx(regression_results["stiffness_std"])


@pytest.mark.parametrize(
    ["hamiltonian_name", "J", "case_iter"],
    [
        ["XY", 1.0, 6],
        ["fm_heisenberg_fm_Z", 1.0, 7],
        ["afm_heisenberg_fm_Z", -1.0, 8],
    ],
)
def test_nearest_neighbor(hamiltonian_name, J, case_iter, tmp_path):

    # Nearest Neighbor QMC
    input_file = "tests/input_files/nearest_neighbor.txt"
    inputs = InputReader(input_file_path=input_file)
    inputs.read_inputs_from_file()
    inputs.read_kwarg_inputs(hamiltonian_name=hamiltonian_name, J=J, root_folder=tmp_path)
    parameter_set_list = inputs.parameter_set_list
    preprocessor = PreprocessorFactory.create("nearest_neighbor_qmc", parameter_set_list)
    driver_inputs = preprocessor.preprocess()
    driver = DriverFactory.create("nearest_neighbor_qmc", preprocessor.simulation_folder, driver_inputs)
    driver.simulate()
    qmc_outputs_file = os.path.join(
        preprocessor.processed_configs[0]["misc"]["run_folder"], "qmc_output/estimators.csv"
    )
    estimator_outputs = compute_qmc_stats(qmc_outputs_file)

    regression_results = regression_data[case_iter]

    assert estimator_outputs[0][1] == pytest.approx(regression_results["energy_mean"])
    assert estimator_outputs[1][1] == pytest.approx(regression_results["energy_std"])
    assert estimator_outputs[0][2] == pytest.approx(regression_results["magnetization_mean"])
    assert estimator_outputs[1][2] == pytest.approx(regression_results["magnetization_std"])


def compute_qmc_stats(qmc_outputs_file):
    estimator_outputs = np.loadtxt(qmc_outputs_file, delimiter=",", skiprows=2)
    return np.mean(estimator_outputs, axis=0), np.std(estimator_outputs, axis=0)
