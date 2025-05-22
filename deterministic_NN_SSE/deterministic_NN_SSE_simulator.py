import numpy as np
from deterministic_NN_SSE_initializer import *
import statistical_analysis as sa

def populate_linked_list(Sm_array, spin_array, M, N, n, sites, vertex_map, total_num_legs):

    # temp arrays to store indices identifying the first and last vertices for each spin
    first_vertex_leg = -1 * np.ones(N, dtype=int)
    last_vertex_leg = -1 * np.ones(N, dtype=int)

    p = 0
    p_list = np.zeros(n, dtype=int)
    b_list = np.zeros(n, dtype=int)

    vertex_array = np.zeros(n, dtype=int)
    linked_list = np.zeros(total_num_legs, dtype=int)

    # TODO: optimization: can likely populate p_list and b_list in diagonal updates and iterate over n<M here instead of M
    # For populating vertex_array and linked_list
    for t in range(M):
        if Sm_array[t] != 0:

            bond_num = Sm_array[t] // 2
            bond_index = bond_num-1
            i_b = sites[bond_index,0]
            j_b = sites[bond_index,1]

            p_list[p] = t
            b_list[p] = bond_num

            # set first_vertex_leg if site never encountered before
            if last_vertex_leg[i_b] == -1:
                first_vertex_leg[i_b] = num_legs_per_vertex*p
                last_vertex_leg[i_b] = num_legs_per_vertex*p+2 # exit leg identifier requires +2
            else: # can link to previous vertex leg for this site if encountered before
                linked_list[num_legs_per_vertex*p] = last_vertex_leg[i_b]
                linked_list[last_vertex_leg[i_b]] = num_legs_per_vertex*p # all links are bi-directional
                last_vertex_leg[i_b] = num_legs_per_vertex*p+2

            # Need to follow this logic for both legs
            if last_vertex_leg[j_b] == -1:
                first_vertex_leg[j_b] = num_legs_per_vertex*p+1
                last_vertex_leg[j_b] = num_legs_per_vertex*p+3
            else:
                linked_list[num_legs_per_vertex*p+1] = last_vertex_leg[j_b]
                linked_list[last_vertex_leg[j_b]] = num_legs_per_vertex*p+1
                last_vertex_leg[j_b] = num_legs_per_vertex*p+3

            # propagate for off-diagonal operators, set vertex in both cases 
            if Sm_array[t] % 2 == 1:
                vertex_array[p] = vertex_map[(spin_array[i_b], spin_array[j_b], -spin_array[i_b], -spin_array[j_b])]
                spin_array[i_b] = -spin_array[i_b]
                spin_array[j_b] = -spin_array[j_b]
            else:
                vertex_array[p] = vertex_map[(spin_array[i_b], spin_array[j_b], spin_array[i_b], spin_array[j_b])]
            
            p+=1

    # Periodic due to trace - connect the first and last spins in expansion direction
    for i in range(N):
        if last_vertex_leg[i] != -1:
            linked_list[last_vertex_leg[i]] = first_vertex_leg[i]
            linked_list[first_vertex_leg[i]] = last_vertex_leg[i]

    return first_vertex_leg, linked_list, vertex_array, p_list, b_list

def update_operators_and_spins(Sm_array, spin_array, M, N, leg_spin, vertex_array, first_vertex_leg):

    p=0
    for t in range(M):
        if Sm_array[t] != 0:
            b = Sm_array[t]//2
            Sm_array[t] = 2*b + isotropic_operator_type[vertex_array[p]]
            p += 1
        
    for s in range(N):
        # Flip disconnected spins with probability 1/2
        if first_vertex_leg[s] == -1:
            spin_array[s] = np.random.choice([1,-1]) * spin_array[s]
        else:
            # Update connected spins according to the first vertex acting on them - consistency in (1+1)D space ensured by construction
            q = first_vertex_leg[s] // num_legs_per_vertex
            l = first_vertex_leg[s] % num_legs_per_vertex
            spin_array[s] = leg_spin[vertex_array[q]][l]

    return Sm_array, spin_array

def off_diagonal_updates_isotropic(Sm_array, spin_array, M, N, n, sites, isotropic_exit_leg_map, vertex_map, new_vertex_map, leg_spin):

    # For a discussion of the algorithm, see Appendix A of: https://arxiv.org/pdf/cond-mat/0202316

    if n == 0:
        for s in range(N):
            spin_array[s] = np.random.choice([1,-1]) * spin_array[s]

        return Sm_array, spin_array
    
    else:
        total_num_legs = num_legs_per_vertex*n
        
        first_vertex_leg, linked_list, vertex_array = populate_linked_list(Sm_array, spin_array, M, N, n, sites, vertex_map, total_num_legs)

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
                    p = j // num_legs_per_vertex
                    l_e = j % num_legs_per_vertex

                    # only one exit leg (can only move horizontally along bars for the isotropic case)
                    l_x = isotropic_exit_leg_map[l_e]
                    vertex_type = vertex_array[p]
                
                    composite_leg_index = num_legs_per_vertex*l_e

                    if flip_spins == 1:
                        # Update vertex by flipping the input leg and output leg spins
                        vertex_array[p] = new_vertex_map[vertex_type, composite_leg_index + l_x]
                    
                    # Continue along linked list to get next j, close loop if required
                    j = num_legs_per_vertex*p + l_x
                    disconnected_leg_flags[j] = 0
                    if j == j_0:
                        loop_closed = True
                    else:
                        j = linked_list[j]
                        if j == j_0:
                            loop_closed = True

        Sm_array, spin_array = update_operators_and_spins(Sm_array, spin_array, M, N, leg_spin, vertex_array, first_vertex_leg)
        
        return Sm_array, spin_array
    
def off_diagonal_updates_XY(Sm_array, spin_array, M, N, n, sites, xy_exit_leg_map, vertex_map, new_vertex_map, leg_spin):

    # For a discussion of the algorithm, see Appendix A of: https://arxiv.org/pdf/cond-mat/0202316

    if n == 0:
        for s in range(N):
            spin_array[s] = np.random.choice([1,-1]) * spin_array[s]

        return Sm_array, spin_array
    
    else:
        total_num_legs = num_legs_per_vertex*n

        first_vertex_leg, linked_list, vertex_array = populate_linked_list(Sm_array, spin_array, M, N, n, sites, vertex_map, total_num_legs)

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
                    p = j // num_legs_per_vertex
                    l_e = j % num_legs_per_vertex

                    # only one exit leg (can only move horizontally along bars for the isotropic case)
                    vertex_type = vertex_array[p]
                    l_x = np.random.choice(xy_exit_leg_map[vertex_type * num_legs_per_vertex + l_e,:])
                
                    composite_leg_index = num_legs_per_vertex*l_e

                    if flip_spins == 1:
                        # Update vertex by flipping the input leg and output leg spins
                        vertex_array[p] = new_vertex_map[vertex_type, composite_leg_index + l_x]
                    
                    # Continue along linked list to get next j, close loop if required
                    j = num_legs_per_vertex*p + l_x
                    disconnected_leg_flags[j] = 0
                    if j == j_0:
                        loop_closed = True
                    else:
                        j = linked_list[j]
                        if j == j_0:
                            loop_closed = True

        Sm_array, spin_array = update_operators_and_spins(Sm_array, spin_array, M, N, leg_spin, vertex_array, first_vertex_leg)
        
        return Sm_array, spin_array

def diagonal_updates_isotropic(Sm_array, spin_array, sites, M, n, beta, num_bonds, diag_update_sign):

    for t in range(M):
        if Sm_array[t] == 0:
            # propose adding a diagonal operator
            acceptance_prob = beta * 0.5 * num_bonds/(M-n)
            u = np.random.uniform(0.0,1.0)
            if u <= acceptance_prob:
                b_prime = np.random.choice(num_bonds) # Gibbs sampling to determine which bond to potentially populate
                i_b = sites[b_prime,0]
                j_b = sites[b_prime,1]
                if spin_array[i_b] == diag_update_sign * spin_array[j_b]: 
                    bond_num = b_prime + 1
                    Sm_array[t] = 2*bond_num
                    n += 1
        else:
            if Sm_array[t] % 2 == 1:
                # propagate spin array for off-diagonal operators
                b = Sm_array[t] // 2
                bond_index = b-1
                i_b = sites[bond_index,0]
                j_b = sites[bond_index,1]
                spin_array[i_b] = -spin_array[i_b]
                spin_array[j_b] = -spin_array[j_b]
            else:
                # propose removing the diagonal operator
                acceptance_prob = (M-n+1) / (beta * 0.5 * num_bonds)
                u = np.random.uniform(0.0,1.0)
                if u <= acceptance_prob:
                    Sm_array[t] = 0
                    n = n - 1

    return n, Sm_array, spin_array

def diagonal_updates_XY(Sm_array, spin_array, sites, M, n, beta, num_bonds):

    for t in range(M):
        if Sm_array[t] == 0:
            # propose adding a diagonal operator
            acceptance_prob = beta * 0.5 * num_bonds/(M-n)
            u = np.random.uniform(0.0,1.0)
            if u <= acceptance_prob:
                b_prime = np.random.choice(num_bonds) # Gibbs sampling to determine which bond to potentially populate
                i_b = sites[b_prime,0]
                j_b = sites[b_prime,1]
                bond_num = b_prime + 1
                Sm_array[t] = 2*bond_num
                n += 1
        else:
            if Sm_array[t] % 2 == 1:
                # propagate spin array for off-diagonal operators
                b = Sm_array[t] // 2
                bond_index = b-1
                i_b = sites[bond_index,0]
                j_b = sites[bond_index,1]
                spin_array[i_b] = -spin_array[i_b]
                spin_array[j_b] = -spin_array[j_b]
            else:
                # propose removing the diagonal operator
                acceptance_prob = (M-n+1) / (beta * 0.5 * num_bonds)
                u = np.random.uniform(0.0,1.0)
                if u <= acceptance_prob:
                    Sm_array[t] = 0
                    n = n - 1

    return n, Sm_array, spin_array

def simulate_isotropic(N, M, equilibration_steps, mc_steps, sites, beta, interaction_type, boundary_conditions, a=1.25):

    # Ferro => 0, Anti-ferro => 1 
    interaction_type_index = isotropic_interaction_type_map[interaction_type]
    exit_leg_map = isotropic_exit_leg_maps[interaction_type_index]
    diagonal_update_condition = (-1)**interaction_type_index

    vertex_mapping = isotropic_vertex_map[interaction_type_index]
    leg_spin = isotropic_leg_spin[interaction_type_index]
    num_vertices = num_vertices_isotropic
    new_vertex_types = generate_new_vertex_type_array(vertex_mapping, num_vertices, leg_spin)

    num_bonds = get_number_of_bonds(N, boundary_conditions)
    spectrum_offset = num_bonds/4.0

    spin_array = np.zeros(N, dtype=int)
    for i in range(N):
        spin_array[i] = (-1)**(interaction_type_index * i)
    Sm_array = np.zeros(M, dtype=int)

    # Equilibration step
    n=0
    n_array = []
    for step in range(equilibration_steps):
        n_new, Sm_array, spin_array = diagonal_updates_isotropic(Sm_array, spin_array, sites, M, n, beta, num_bonds, diagonal_update_condition)

        # track n to increase M as necessary
        n_array.append(n_new)
        M_new = max(int(a * np.max(n_array)),M)

        Sm_array_new = np.zeros(M_new, dtype=int)
        Sm_array_new[:M] = Sm_array
        Sm_array = Sm_array_new

        n = n_new
        M = M_new

        Sm_array, spin_array = off_diagonal_updates_isotropic(Sm_array, spin_array, M, N, n, sites, exit_leg_map, vertex_mapping, new_vertex_types, leg_spin)

    n_array = np.zeros(mc_steps)
    for step in range(mc_steps):
        # diagonal updates first
        n, Sm_array, spin_array = diagonal_updates_isotropic(Sm_array, spin_array, sites, M, n, beta, num_bonds, diagonal_update_condition)
        # off-diagonal updates
        Sm_array, spin_array = off_diagonal_updates_isotropic(Sm_array, spin_array, M, N, n, sites, exit_leg_map, vertex_mapping, new_vertex_types, leg_spin)

        n_array[step] = n

    energy_array = -n_array/beta

    energy_mean, energy_error = sa.statistics_binning(energy_array) # use binning for errors etc., can likely improve this for more sensitive estimators

    return energy_array, energy_mean, energy_error, spectrum_offset

def simulate_XY(N, M, equilibration_steps, mc_steps, sites, beta, boundary_conditions, a=1.25):

    vertex_mapping = general_vertex_map
    leg_spin = general_leg_spin
    num_vertices = num_vertices_general
    new_vertex_types = generate_new_vertex_type_array(vertex_mapping, num_vertices, leg_spin)
    exit_leg_map = generate_XY_exit_legs(new_vertex_types)

    num_bonds = get_number_of_bonds(N, boundary_conditions)
    spectrum_offset = num_bonds/2.0

    spin_array = np.ones(N, dtype=int)
    
    for i in range(N):
        #spin_array[i] = np.random.choice([1,-1])
        spin_array[i] = (-1)**i
    
    Sm_array = np.zeros(M, dtype=int)

    # Equilibration step
    n=0
    n_array = []
    for step in range(equilibration_steps):
        n_new, Sm_array, spin_array = diagonal_updates_XY(Sm_array, spin_array, sites, M, n, beta, num_bonds)

        # track n to increase M as necessary
        n_array.append(n_new)
        M_new = max(int(a * np.max(n_array)),M)

        Sm_array_new = np.zeros(M_new, dtype=int)
        Sm_array_new[:M] = Sm_array
        Sm_array = Sm_array_new

        n = n_new
        M = M_new

        Sm_array, spin_array = off_diagonal_updates_XY(Sm_array, spin_array, M, N, n, sites, exit_leg_map, vertex_mapping, new_vertex_types, leg_spin)

    n_array = np.zeros(mc_steps)
    for step in range(mc_steps):
        # diagonal updates first
        n, Sm_array, spin_array = diagonal_updates_XY(Sm_array, spin_array, sites, M, n, beta, num_bonds)
        # off-diagonal updates
        Sm_array, spin_array = off_diagonal_updates_XY(Sm_array, spin_array, M, N, n, sites, exit_leg_map, vertex_mapping, new_vertex_types, leg_spin)

        n_array[step] = n

    energy_array = -n_array/beta

    energy_mean, energy_error = sa.statistics_binning(energy_array) # use binning for errors etc., can likely improve this for more sensitive estimators

    return energy_array, energy_mean, energy_error, spectrum_offset


