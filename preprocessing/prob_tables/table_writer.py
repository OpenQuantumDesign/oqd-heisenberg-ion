import os
import numpy as np
from coupling_strengths import get_J_ij_exp
from coupling_strengths import get_J_ij_power_law
from geometry import geometry
from heatbath_probs import compute_prob_tables_heat_bath
from directed_loop_probs import compute_prob_tables_directed_loops
from deterministic_probs import compute_probability_tables_deterministic

def write_prob_tables(gamma, Delta, h, N_1, N_2, lattice_type, alpha, ksi, J, loop_update_type, dist_dep_gamma, hamiltonian_type, boundary):

    out_dir = "ProbabilityDensities"
    sites, distances, geometry_table, N, num_bonds = geometry(N_1, N_2, lattice_type, boundary)
    num_bonds = int(N*(N-1)/2.0)
    h_B = h/(J*(N-1))

    J_ij_vector = get_J_ij_power_law(num_bonds, distances, alpha)

    if loop_update_type == "heat_bath":
        prob_tables = compute_prob_tables_heat_bath(num_bonds, sites, J_ij_vector, gamma, h_B, Delta)
        diag_prob_table = prob_tables[0]
        max_over_states = prob_tables[1]
        max_diag_norm = prob_tables[2]
        vertex_weights = prob_tables[3]
        spectrum_offset = prob_tables[4]
        loop_update_prob_table = prob_tables[5]

    elif loop_update_type == "directed_loops": 
        prob_tables = compute_prob_tables_directed_loops(num_bonds, sites, J_ij_vector, gamma, h_B, Delta, ksi, dist_dep_gamma)
        diag_prob_table = prob_tables[0]
        max_over_states = prob_tables[1]
        max_diag_norm = prob_tables[2]
        vertex_weights = prob_tables[3]
        spectrum_offset = prob_tables[4]
        loop_update_prob_table = prob_tables[5]

    else:
        raise Exception("Invalid key word argument for loop_update_type")
    
    if dist_dep_gamma:
        dist_dep_offset = 1
    else:
        dist_dep_offset = 0
    
    if (hamiltonian_type != 2):
        raise Exception("This function requires the hamiltonian type to 2\n")

    file_prefix = "N_{}_hamiltonian_type_{}_Delta_{}_h_{}_alpha_{}_gamma_{}_ksi_{}_J_{}_dist_dep_offset_{}_boundary_{}".format(N, hamiltonian_type, Delta, h, alpha, gamma, ksi, J, dist_dep_offset, boundary)

    geometry_file_name = os.path.join(out_dir, "N_{}_geometry.csv".format(N))
    diag_file_name = os.path.join(out_dir, file_prefix + "_diag_probs.csv")
    max_over_states_file_name = os.path.join(out_dir, file_prefix + "_max_over_states.csv")
    loop_update_table_file_name = os.path.join(out_dir, file_prefix + "_{}_off_diag_table.csv".format(loop_update_type))
    vertex_weights_file_name = os.path.join(out_dir, file_prefix + "_vertex_weights.csv".format(loop_update_type))

    np.savetxt(geometry_file_name, geometry_table, delimiter=",", fmt="%d", header="N={}, NumBonds={}".format(N, num_bonds))
    np.savetxt(diag_file_name, diag_prob_table, delimiter=",", header="N={},Delta={},h={},alpha={},gamma={},ksi={},J={}".format(N, Delta, h, alpha, gamma, ksi, J))
    np.savetxt(vertex_weights_file_name, vertex_weights, delimiter=",", header="N={},Delta={},h={},alpha={},gamma={},ksi={},J={}".format(N, Delta, h, alpha, gamma, ksi, J))
    np.savetxt(max_over_states_file_name, max_over_states, delimiter=",", header="N={},Delta={},h={},alpha={},gamma={},ksi={},J={},norm={},,spectrum_offset={}".format(N, Delta, h, alpha, gamma, ksi, J, max_diag_norm, spectrum_offset))
    np.savetxt(loop_update_table_file_name, loop_update_prob_table, delimiter=",", header="N={},Delta={},h={},alpha={},gamma={},ksi={},J={},spectrum_offset={},loop_update_type={}".format(N, Delta, h, alpha, gamma, ksi, J, spectrum_offset, loop_update_type))

    return 0

def write_prob_tables_deterministic(N_1, N_2, lattice_type, J_ij_type, J_ij_file_path, alpha, J, hamiltonian_type, boundary, mu):

    out_dir = "ProbabilityDensities"
    sites, distances, geometry_table, N, num_bonds = geometry(N_1, N_2, lattice_type, boundary)
    num_bonds = int(N*(N-1)/2.0)

    if J_ij_type == 1:
        if alpha != None or J != None:
            raise Exception("alpha and J should be None if J_ij_type is 1\n")
        if J_ij_file_path == None:
            raise Exception("No file path for J_ij matrix provided\n")
        alpha = "exp_mu_{}".format(mu)
        J = "exp"
        J_ij_vector = get_J_ij_exp(N, num_bonds, J_ij_file_path)
    elif J_ij_type == 0:
        if J_ij_file_path != None:
            raise Exception("No file path for J_ij matrix should be provided if J_ij_type is 0\n")
        if alpha == None or J == None:
            raise Exception("Need to provide alpha and J if J_ij_type is 0\n")
        J_ij_vector = get_J_ij_power_law(num_bonds, distances, alpha)

    prob_tables = compute_probability_tables_deterministic(num_bonds, J_ij_vector, hamiltonian_type)

    max_over_states = prob_tables[0]
    max_norm = prob_tables[1]
    spectrum_offset = prob_tables[2]
    Delta = prob_tables[3]

    file_prefix = "N_{}_hamiltonian_type_{}_Delta_{}_h_0.0_alpha_{}_gamma_0.0_ksi_0.0_J_{}_dist_dep_offset_0_boundary_{}".format(N, hamiltonian_type, Delta, alpha, J, boundary)

    geometry_file_name = os.path.join(out_dir, "N_{}_geometry_boundary_{}.csv".format(N, boundary))
    max_over_states_file_name = os.path.join(out_dir, file_prefix + "_max_over_states.csv")

    np.savetxt(geometry_file_name, geometry_table, delimiter=",", fmt="%d", header="N={}, NumBonds={}".format(N, num_bonds))
    np.savetxt(max_over_states_file_name, max_over_states, delimiter=",", header="N={},Delta={},alpha={},J={},norm={},spectrum_offset={}".format(N, Delta, alpha, J, max_norm, spectrum_offset))

    return 0