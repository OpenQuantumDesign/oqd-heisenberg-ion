import os

import numpy as np

from . import vertex_types as vt


class ConfigurationGenerator:
    def __init__(self, system_inputs, sampling_inputs, run_folder):

        self.run_folder = run_folder
        self.out_dir = run_folder + "/qmc_output/"
        os.mkdir(self.out_dir)

        # self.geometry = system_inputs.geometry
        self.hamiltonian = system_inputs.model_name
        self.Delta = system_inputs.hamiltonian_parameters.Delta

        self.init_M = sampling_inputs.initial_operator_list_size
        if self.init_M == -1:
            self.init_M = 50
        self.equilibration_steps = sampling_inputs.equilibration_steps
        self.mc_steps = sampling_inputs.mc_steps

        self.T = sampling_inputs.T
        self.J = system_inputs.hamiltonian_parameters.J
        self.beta = sampling_inputs.beta * np.abs(self.J)  # beta times |J| determines the energy scale
        self.a = sampling_inputs.operator_list_update_multiplier
        self.init_config_index = sampling_inputs.initial_configuration_index

        self.N = system_inputs.geometry.N
        self.num_bonds = system_inputs.geometry.num_bonds
        self.boundary = system_inputs.geometry.boundary

        self.sites = system_inputs.geometry.sites

        self.write_simulation_specs()

    def write_simulation_specs(self):

        specs_file_path = os.path.join(self.run_folder, "simulation_specs.txt")
        with open(specs_file_path, "w") as f:
            f.write(f"output_folder\t{self.run_folder}\n")
            f.write(f"N\t{self.N}\n")
            f.write(f"num_bonds\t{self.num_bonds}\n")
            f.write(f"hamiltonian\t{self.hamiltonian}\n")
            f.write(f"Delta\t{self.Delta}\n")
            f.write(f"boundary\t{self.boundary}\n")
            f.write(f"initial_M\t{self.init_M}\n")
            f.write(f"equilibration_steps\t{self.equilibration_steps}\n")
            f.write(f"mc_steps\t{self.mc_steps}\n")
            f.write(f"T\t{self.T}\n")
            f.write(f"J\t{self.J}\n")
            f.write(f"beta\t{self.beta}\n")
            f.write(f"a\t{self.a}\n")
            f.write(f"init_config_index\t{self.init_config_index}\n")

        return 0

    def initialize_spin_array(self):

        N = self.N

        if self.init_config_index == 1 or self.init_config_index == -1:
            spin_array = self.init_config_index * np.ones(N, dtype=int)

        elif self.init_config_index == 0:
            spin_array = np.zeros(N, dtype=int)
            for i in range(N):
                spin_array[i] = np.random.choice([1, -1])

        elif self.init_config_index > 1:
            spin_array = -1 * np.ones(N, dtype=int)
            for i in range(N):
                if i % self.init_config_index == 0:
                    spin_array[i] = 1

        else:
            raise Exception(f"Unrecognized initial spin configuration index: {self.init_config_index}")

        return spin_array

    def update_diagonal_estimators(self, spin_array, magnetization, step, N):

        M_z = 0.0
        for i in range(N):
            M_z += spin_array[i]

        M_z *= 0.5 / N
        magnetization[step] = M_z

        return magnetization

    def populate_linked_list(self, Sm_array, spin_array, M, N, n, vertex_map, total_num_legs):

        # temp arrays to store indices identifying the first and last vertices for each spin
        first_vertex_leg = -1 * np.ones(N, dtype=int)
        last_vertex_leg = -1 * np.ones(N, dtype=int)

        p = 0
        p_list = np.zeros(n, dtype=int)
        b_list = np.zeros(n, dtype=int)

        vertex_array = np.zeros(n, dtype=int)
        linked_list = np.zeros(total_num_legs, dtype=int)

        # TODO: optimization: can likely populate p_list and b_list in diagonal updates and iterate over n<M
        # here instead of M
        # For populating vertex_array and linked_list
        for t in range(M):
            if Sm_array[t] != 0:
                bond_num = Sm_array[t] // 2
                bond_index = bond_num - 1
                i_b = self.sites[bond_index, 0]
                j_b = self.sites[bond_index, 1]

                p_list[p] = t
                b_list[p] = bond_num

                # set first_vertex_leg if site never encountered before
                if last_vertex_leg[i_b] == -1:
                    first_vertex_leg[i_b] = vt.num_legs_per_vertex * p
                    last_vertex_leg[i_b] = vt.num_legs_per_vertex * p + 2  # exit leg identifier requires +2
                else:  # can link to previous vertex leg for this site if encountered before
                    linked_list[vt.num_legs_per_vertex * p] = last_vertex_leg[i_b]
                    linked_list[last_vertex_leg[i_b]] = vt.num_legs_per_vertex * p  # all links are bi-directional
                    last_vertex_leg[i_b] = vt.num_legs_per_vertex * p + 2

                # Need to follow this logic for both legs
                if last_vertex_leg[j_b] == -1:
                    first_vertex_leg[j_b] = vt.num_legs_per_vertex * p + 1
                    last_vertex_leg[j_b] = vt.num_legs_per_vertex * p + 3
                else:
                    linked_list[vt.num_legs_per_vertex * p + 1] = last_vertex_leg[j_b]
                    linked_list[last_vertex_leg[j_b]] = vt.num_legs_per_vertex * p + 1
                    last_vertex_leg[j_b] = vt.num_legs_per_vertex * p + 3

                # propagate for off-diagonal operators, set vertex in both cases
                if Sm_array[t] % 2 == 1:
                    vertex_array[p] = vertex_map[(spin_array[i_b], spin_array[j_b], -spin_array[i_b], -spin_array[j_b])]
                    spin_array[i_b] = -spin_array[i_b]
                    spin_array[j_b] = -spin_array[j_b]
                else:
                    vertex_array[p] = vertex_map[(spin_array[i_b], spin_array[j_b], spin_array[i_b], spin_array[j_b])]

                p += 1

        # Periodic due to trace - connect the first and last spins in expansion direction
        for i in range(N):
            if last_vertex_leg[i] != -1:
                linked_list[last_vertex_leg[i]] = first_vertex_leg[i]
                linked_list[first_vertex_leg[i]] = last_vertex_leg[i]

        return first_vertex_leg, linked_list, vertex_array, p_list, b_list

    def update_operators_and_spins(self, Sm_array, spin_array, M, N, leg_spin, vertex_array, first_vertex_leg, op_type):

        p = 0
        for t in range(M):
            if Sm_array[t] != 0:
                b = Sm_array[t] // 2
                Sm_array[t] = 2 * b + op_type[vertex_array[p]]
                p += 1

        for s in range(N):
            # Flip disconnected spins with probability 1/2
            if first_vertex_leg[s] == -1:
                spin_array[s] = np.random.choice([1, -1]) * spin_array[s]
            else:
                # Update connected spins according to the first vertex acting on them
                # consistency in (1+1)D space ensured by construction
                q = first_vertex_leg[s] // vt.num_legs_per_vertex
                l = first_vertex_leg[s] % vt.num_legs_per_vertex
                spin_array[s] = leg_spin[vertex_array[q]][l]

        return Sm_array, spin_array

    def off_diagonal_updates_isotropic(
        self, Sm_array, spin_array, M, N, n, isotropic_exit_leg_map, vertex_map, new_vertex_map, leg_spin
    ):

        # For a discussion of the algorithm, see Appendix A of: https://arxiv.org/pdf/cond-mat/0202316

        if n == 0:
            for s in range(N):
                spin_array[s] = np.random.choice([1, -1]) * spin_array[s]

            return Sm_array, spin_array

        else:
            total_num_legs = vt.num_legs_per_vertex * n

            first_vertex_leg, linked_list, vertex_array, p_list, b_list = self.populate_linked_list(
                Sm_array, spin_array, M, N, n, vertex_map, total_num_legs
            )

            # Construct deterministic loops
            disconnected_leg_flags = np.ones(total_num_legs, dtype=int)
            for j_0 in range(total_num_legs):
                if disconnected_leg_flags[j_0] == 1:
                    loop_closed = False
                    j = j_0
                    flip_spins = np.random.choice(2)

                    while not loop_closed:
                        disconnected_leg_flags[j] = 0

                        # Associated vertex and entrance leg
                        p = j // vt.num_legs_per_vertex
                        l_e = j % vt.num_legs_per_vertex

                        # only one exit leg (can only move horizontally along bars for the isotropic case)
                        l_x = isotropic_exit_leg_map[l_e]
                        vertex_type = vertex_array[p]

                        composite_leg_index = vt.num_legs_per_vertex * l_e

                        if flip_spins == 1:
                            # Update vertex by flipping the input leg and output leg spins
                            vertex_array[p] = new_vertex_map[vertex_type, composite_leg_index + l_x]

                        # Continue along linked list to get next j, close loop if required
                        j = vt.num_legs_per_vertex * p + l_x
                        disconnected_leg_flags[j] = 0
                        j = linked_list[j]
                        if j == j_0:
                            loop_closed = True

            Sm_array, spin_array = self.update_operators_and_spins(
                Sm_array, spin_array, M, N, leg_spin, vertex_array, first_vertex_leg, vt.isotropic_operator_type
            )

            return Sm_array, spin_array

    def off_diagonal_updates_XY(
        self, Sm_array, spin_array, M, N, n, xy_exit_leg_map, vertex_map, new_vertex_map, leg_spin
    ):

        # For a discussion of the algorithm, see Appendix A of: https://arxiv.org/pdf/cond-mat/0202316
        if n == 0:
            for s in range(N):
                spin_array[s] = np.random.choice([1, -1]) * spin_array[s]

            return Sm_array, spin_array

        else:
            total_num_legs = vt.num_legs_per_vertex * n

            first_vertex_leg, linked_list, vertex_array, p_list, b_list = self.populate_linked_list(
                Sm_array, spin_array, M, N, n, vertex_map, total_num_legs
            )

            # Construct deterministic loops
            disconnected_leg_flags = np.ones(total_num_legs, dtype=int)
            for j_0 in range(total_num_legs):
                if disconnected_leg_flags[j_0] == 1:
                    loop_closed = False
                    j = j_0
                    flip_spins = np.random.choice(2)

                    while not loop_closed:
                        disconnected_leg_flags[j] = 0

                        # Associated vertex and entrance leg
                        p = j // vt.num_legs_per_vertex
                        l_e = j % vt.num_legs_per_vertex

                        # only one exit leg (can only move horizontally along bars for the isotropic case)
                        vertex_type = vertex_array[p]
                        l_x = np.random.choice(xy_exit_leg_map[vertex_type * vt.num_legs_per_vertex + l_e, :])

                        composite_leg_index = vt.num_legs_per_vertex * l_e

                        if flip_spins == 1:
                            # Update vertex by flipping the input leg and output leg spins
                            vertex_array[p] = new_vertex_map[vertex_type, composite_leg_index + l_x]

                        # Continue along linked list to get next j, close loop if required
                        j = vt.num_legs_per_vertex * p + l_x
                        disconnected_leg_flags[j] = 0
                        j = linked_list[j]
                        if j == j_0:
                            loop_closed = True

            Sm_array, spin_array = self.update_operators_and_spins(
                Sm_array, spin_array, M, N, leg_spin, vertex_array, first_vertex_leg, vt.general_operator_type
            )

            return Sm_array, spin_array

    def diagonal_updates_isotropic(self, Sm_array, spin_array, M, n, num_bonds, diag_update_sign):

        sites = self.sites
        for t in range(M):
            if Sm_array[t] == 0:
                # propose adding a diagonal operator
                acceptance_prob = self.beta * 0.5 * num_bonds / (M - n)
                u = np.random.uniform(0.0, 1.0)
                if u <= acceptance_prob:
                    b_prime = np.random.choice(
                        num_bonds
                    )  # Gibbs sampling to determine which bond to potentially populate
                    i_b = sites[b_prime, 0]
                    j_b = sites[b_prime, 1]
                    if spin_array[i_b] == diag_update_sign * spin_array[j_b]:
                        bond_num = b_prime + 1
                        Sm_array[t] = 2 * bond_num
                        n += 1
            else:
                if Sm_array[t] % 2 == 1:
                    # propagate spin array for off-diagonal operators
                    b = Sm_array[t] // 2
                    bond_index = b - 1
                    i_b = sites[bond_index, 0]
                    j_b = sites[bond_index, 1]
                    spin_array[i_b] = -spin_array[i_b]
                    spin_array[j_b] = -spin_array[j_b]
                else:
                    # propose removing the diagonal operator
                    acceptance_prob = (M - n + 1) / (self.beta * 0.5 * num_bonds)
                    u = np.random.uniform(0.0, 1.0)
                    if u <= acceptance_prob:
                        Sm_array[t] = 0
                        n = n - 1

        return n, Sm_array, spin_array

    def diagonal_updates_XY(self, Sm_array, spin_array, M, n, num_bonds):

        sites = self.sites

        for t in range(M):
            if Sm_array[t] == 0:
                # propose adding a diagonal operator
                acceptance_prob = self.beta * 0.5 * num_bonds / (M - n)
                u = np.random.uniform(0.0, 1.0)
                if u <= acceptance_prob:
                    b_prime = np.random.choice(
                        num_bonds
                    )  # Gibbs sampling to determine which bond to potentially populate
                    i_b = sites[b_prime, 0]
                    j_b = sites[b_prime, 1]
                    bond_num = b_prime + 1
                    Sm_array[t] = 2 * bond_num
                    n += 1
            else:
                if Sm_array[t] % 2 == 1:
                    # propagate spin array for off-diagonal operators
                    b = Sm_array[t] // 2
                    bond_index = b - 1
                    i_b = sites[bond_index, 0]
                    j_b = sites[bond_index, 1]
                    spin_array[i_b] = -spin_array[i_b]
                    spin_array[j_b] = -spin_array[j_b]
                else:
                    # propose removing the diagonal operator
                    acceptance_prob = (M - n + 1) / (self.beta * 0.5 * num_bonds)
                    u = np.random.uniform(0.0, 1.0)
                    if u <= acceptance_prob:
                        Sm_array[t] = 0
                        n = n - 1

        return n, Sm_array, spin_array

    def simulate_isotropic(self):

        # Ferro => 0, Anti-ferro => 1
        interaction_type_index = vt.isotropic_interaction_type_map[self.hamiltonian]
        exit_leg_map = vt.isotropic_exit_leg_maps[interaction_type_index]
        diagonal_update_condition = (-1) ** interaction_type_index

        vertex_mapping = vt.isotropic_vertex_map[interaction_type_index]
        leg_spin = vt.isotropic_leg_spin[interaction_type_index]
        num_vertices = vt.num_vertices_isotropic
        new_vertex_types = vt.generate_new_vertex_type_array(vertex_mapping, num_vertices, leg_spin)

        num_bonds = self.num_bonds
        N = self.N
        spectrum_offset = num_bonds / 4.0
        M = self.init_M

        magnetization_array = np.zeros(self.mc_steps)

        step_spin = np.zeros((self.mc_steps, N), dtype=int)

        spin_array = self.initialize_spin_array()

        Sm_array = np.zeros(M, dtype=int)

        # Equilibration step
        n = 0
        n_array = []
        for step in range(self.equilibration_steps):
            n_new, Sm_array, spin_array = self.diagonal_updates_isotropic(
                Sm_array, spin_array, M, n, num_bonds, diagonal_update_condition
            )

            # track n to increase M as necessary
            n_array.append(n_new)
            M_new = max(int(self.a * np.max(n_array)), M)

            Sm_array_new = np.zeros(M_new, dtype=int)
            Sm_array_new[:M] = Sm_array
            Sm_array = Sm_array_new

            n = n_new
            M = M_new

            Sm_array, spin_array = self.off_diagonal_updates_isotropic(
                Sm_array, spin_array, M, N, n, exit_leg_map, vertex_mapping, new_vertex_types, leg_spin
            )

        n_array = np.zeros(self.mc_steps)
        for step in range(self.mc_steps):
            # diagonal updates first
            n, Sm_array, spin_array = self.diagonal_updates_isotropic(
                Sm_array, spin_array, M, n, num_bonds, diagonal_update_condition
            )
            # off-diagonal updates
            Sm_array, spin_array = self.off_diagonal_updates_isotropic(
                Sm_array, spin_array, M, N, n, exit_leg_map, vertex_mapping, new_vertex_types, leg_spin
            )
            self.update_diagonal_estimators(spin_array, magnetization_array, step, N)
            step_spin[step, :] = spin_array

            n_array[step] = n

        energy_array = -n_array / self.beta + spectrum_offset

        self.spin_array = spin_array
        self.energy_array = energy_array
        self.magnetization_array = magnetization_array

        return 0

    def simulate_XY(self):

        vertex_mapping = vt.general_vertex_map
        leg_spin = vt.general_leg_spin
        num_vertices = vt.num_vertices_general
        new_vertex_types = vt.generate_new_vertex_type_array(vertex_mapping, num_vertices, leg_spin)
        exit_leg_map = vt.generate_XY_exit_legs(new_vertex_types)

        num_bonds = self.num_bonds
        spectrum_offset = num_bonds / 2.0

        magnetization_array = np.zeros(self.mc_steps)

        N = self.N
        M = self.init_M

        spin_array = self.initialize_spin_array()

        Sm_array = np.zeros(M, dtype=int)
        step_spin = np.zeros((self.mc_steps, N), dtype=int)

        # Equilibration step
        n = 0
        n_array = []
        for step in range(self.equilibration_steps):
            n_new, Sm_array, spin_array = self.diagonal_updates_XY(Sm_array, spin_array, M, n, num_bonds)

            # track n to increase M as necessary
            n_array.append(n_new)
            M_new = max(int(self.a * np.max(n_array)), M)

            Sm_array_new = np.zeros(M_new, dtype=int)
            Sm_array_new[:M] = Sm_array
            Sm_array = Sm_array_new

            n = n_new
            M = M_new

            Sm_array, spin_array = self.off_diagonal_updates_XY(
                Sm_array, spin_array, M, N, n, exit_leg_map, vertex_mapping, new_vertex_types, leg_spin
            )

        n_array = np.zeros(self.mc_steps)
        for step in range(self.mc_steps):
            # diagonal updates first
            n, Sm_array, spin_array = self.diagonal_updates_XY(Sm_array, spin_array, M, n, num_bonds)
            # off-diagonal updates
            Sm_array, spin_array = self.off_diagonal_updates_XY(
                Sm_array, spin_array, M, N, n, exit_leg_map, vertex_mapping, new_vertex_types, leg_spin
            )
            self.update_diagonal_estimators(spin_array, magnetization_array, step, N)
            step_spin[step, :] = spin_array

            n_array[step] = n

        energy_array = -n_array / self.beta + spectrum_offset

        self.spin_array = spin_array
        self.energy_array = energy_array
        self.magnetization_array = magnetization_array

        return 0

    def simulate(self):

        if self.hamiltonian == "XY":
            self.simulate_XY()
        elif self.hamiltonian == "fm_heisenberg_fm_Z":
            self.simulate_isotropic()
        elif self.hamiltonian == "afm_heisenberg_fm_Z":
            self.simulate_isotropic()

    def write_outputs(self):

        spin_configs_file = os.path.join(self.out_dir, "final_spin_configurations.csv")
        step_outputs_file = os.path.join(self.out_dir, "estimators.csv")

        np.savetxt(spin_configs_file, self.spin_array, header="final_spin_configuration", delimiter=",", fmt="%d")

        with open(step_outputs_file, "w") as f:
            header = "MC Step Outputs\n"
            header += "MC Step, Energy, Magnetization"
            f.write(header)
            for i in range(self.mc_steps):
                f.write("\n")
                f.write(str(i + 1) + "," + str(self.energy_array[i]) + "," + str(self.magnetization_array[i]))
