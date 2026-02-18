import os
from pathlib import Path

import numpy as np
import pytest

from heisenberg_ion.common.postprocess import utils as post_proc
from heisenberg_ion.interface.main import main


@pytest.mark.parametrize(
    ["hamiltonian_name", "alpha"],
    [
        ["fm_heisenberg_afm_Z", "1.0"],
        ["fm_heisenberg_fm_Z", "2.0"],
        ["XXZ", "1.5"],
        ["XY", "1.0"],
        ["XY", "3.0"],
        ["XY", "10.0"],
    ],
)
def test_long_range_qmc(hamiltonian_name, alpha, tmp_path):

    test_file_dir = Path(__file__).parent
    T = "0.005"

    input_file = os.path.join(test_file_dir, "input_files/" + hamiltonian_name + "_long_range.txt")
    # out_dir = os.path.join(test_file_dir, "results/long_range_qmc")
    out_dir = str(tmp_path)
    uuid = "test_calculation"

    qmc_out_folder = "long_range_qmc_alpha_" + alpha + "_" + hamiltonian_name + "_results"
    return_code_qmc = simulate(
        input_file=input_file,
        root_folder=out_dir,
        output_folder_name=qmc_out_folder,
        hamiltonian_name=hamiltonian_name,
        alpha=alpha,
        simulator="long_range_qmc",
        uuid=uuid,
        T=T,
    )

    assert return_code_qmc == 0

    ed_out_folder = "ed_alpha_" + alpha + "_" + hamiltonian_name + "_results"
    return_code_ed = simulate(
        input_file=input_file,
        root_folder=out_dir,
        output_folder_name=ed_out_folder,
        hamiltonian_name=hamiltonian_name,
        alpha=alpha,
        simulator="exact_diagonalization",
        uuid=uuid,
        T=T,
    )

    assert return_code_ed == 0

    qmc_energy_mean, qmc_energy_err = compute_qmc_energy(out_dir, qmc_out_folder, uuid)
    ground_state_energy = compute_ed_energy(out_dir, ed_out_folder, uuid, T)

    print(qmc_energy_mean + qmc_energy_err)
    print(qmc_energy_mean - qmc_energy_err)

    print(ground_state_energy)

    assert qmc_energy_mean == pytest.approx(ground_state_energy, abs=0.05)


def test_long_range_qmc_exp_J_ij(tmp_path):

    test_file_dir = Path(__file__).parent
    # print(test_file_dir)
    hamiltonian_name = "XY"
    T = "0.005"
    interaction_matrix_file = (
        "/Users/shaeermoeed/Github/Heisenberg_Ion/tests/integration/input_files/Experimental_J_ij_mu_1.01_N_5.csv"
    )

    input_file = os.path.join(test_file_dir, "input_files/" + hamiltonian_name + "_long_range.txt")
    print(input_file)
    # out_dir = os.path.join(test_file_dir, "results/long_range_qmc")
    out_dir = str(tmp_path)
    uuid = "test_calculation"
    print(interaction_matrix_file)

    qmc_out_folder = "long_range_qmc" + "_" + hamiltonian_name + "_results"
    return_code_qmc = simulate(
        input_file=input_file,
        root_folder=out_dir,
        output_folder_name=qmc_out_folder,
        hamiltonian_name=hamiltonian_name,
        interaction_type="matrix_input",
        interaction_matrix_file=interaction_matrix_file,
        simulator="long_range_qmc",
        uuid=uuid,
        T=T,
    )

    assert return_code_qmc == 0

    ed_out_folder = "ed" + "_" + hamiltonian_name + "_results"
    return_code_ed = simulate(
        input_file=input_file,
        root_folder=out_dir,
        output_folder_name=ed_out_folder,
        hamiltonian_name=hamiltonian_name,
        interaction_type="matrix_input",
        interaction_matrix_file=interaction_matrix_file,
        simulator="exact_diagonalization",
        uuid=uuid,
        T=T,
    )

    assert return_code_ed == 0

    qmc_energy_mean, qmc_energy_err = compute_qmc_energy(out_dir, qmc_out_folder, uuid)
    ground_state_energy = compute_ed_energy(out_dir, ed_out_folder, uuid, T)

    print(qmc_energy_mean + qmc_energy_err)
    print(qmc_energy_mean - qmc_energy_err)

    print(ground_state_energy)

    assert qmc_energy_mean == pytest.approx(ground_state_energy, abs=0.05)


@pytest.mark.parametrize(
    ["hamiltonian_name"],
    [
        ["XY"],
        ["fm_heisenberg_fm_Z"],
        ["afm_heisenberg_fm_Z"],
    ],
)
def test_nearest_neighbor_qmc(hamiltonian_name, tmp_path):

    test_file_dir = Path(__file__).parent
    T = "0.1"

    if hamiltonian_name == "afm_heisenberg_fm_Z":
        J = "-1.0"
        N = "6"
        init_config = "0"
    else:
        J = "1.0"
        N = "5"
        init_config = "1"

    input_file = os.path.join(test_file_dir, "input_files/" + hamiltonian_name + "_nearest_neighbor.txt")
    # out_dir = os.path.join(test_file_dir, "results/nearest_neighbor_qmc")
    out_dir = str(tmp_path)
    uuid = "test_calculation"

    qmc_out_folder = "nearest_neighbor_qmc_results_" + hamiltonian_name + "_results"
    return_code_qmc = simulate(
        input_file=input_file,
        root_folder=out_dir,
        output_folder_name=qmc_out_folder,
        hamiltonian_name=hamiltonian_name,
        simulator="nearest_neighbor_qmc",
        uuid=uuid,
        T=T,
        J=J,
        N=N,
        initial_configuration_index=init_config,
        interaction_range="nearest_neighbor",
    )

    assert return_code_qmc == 0

    ed_out_folder = "ed_results_" + hamiltonian_name
    return_code_ed = simulate(
        input_file=input_file,
        root_folder=out_dir,
        output_folder_name=ed_out_folder,
        hamiltonian_name=hamiltonian_name,
        simulator="exact_diagonalization",
        uuid=uuid,
        T=T,
        J=J,
        N=N,
        initial_configuration_index=init_config,
        interaction_range="nearest_neighbor",
    )

    assert return_code_ed == 0

    qmc_energy_mean, qmc_energy_err = compute_qmc_energy(out_dir, qmc_out_folder, uuid)
    ground_state_energy = compute_ed_energy(out_dir, ed_out_folder, uuid, T)

    print(qmc_energy_mean + qmc_energy_err)
    print(qmc_energy_mean - qmc_energy_err)

    assert qmc_energy_mean == pytest.approx(ground_state_energy, abs=0.05)


def simulate(input_file, **kwargs):

    input_list = ["-i", input_file]
    for key, val in kwargs.items():
        input_list.append("-o")
        input_list.append(key)
        input_list.append(val)

    return_code = main(input_list)

    return return_code


def compute_qmc_energy(out_dir, qmc_out_folder, uuid):

    qmc_outputs_file = out_dir + "/" + qmc_out_folder + "/" + uuid + "/" + "qmc_output/estimators.csv"
    estimator_outputs = np.loadtxt(qmc_outputs_file, delimiter=",", skiprows=2)
    qmc_energy_estimates = estimator_outputs[:, 1]
    energy_stats = post_proc.statistics_binning(qmc_energy_estimates, 20, 500)

    return energy_stats[0], energy_stats[1]


def compute_ed_energy(out_dir, ed_out_folder, uuid, T):

    ed_outputs_file = out_dir + "/" + ed_out_folder + "/" + uuid + "/" + "ed_output/energies.csv"
    ed_outputs = np.loadtxt(ed_outputs_file, delimiter=",", skiprows=1)
    ground_state_energy = post_proc.ed_energy(ed_outputs, float(T))

    return ground_state_energy
