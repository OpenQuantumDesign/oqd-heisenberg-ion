import numpy as np
from SSE_initializer import *
import math
import estimators as est

def off_diagonal_updates(Sm_array, spin_array, M, N, N_l, n, max_loop_size, sites, off_diag_prob_table, num_free_flips, p_list, b_list):

    # For a discussion of the algorithm, see Appendix A of: https://arxiv.org/pdf/cond-mat/0202316

    if n == 0:

        for i in range(num_free_flips):
            for s in range(N):
                # Flip disconnected spins with probability 1/2
                spin_array[s] = np.random.choice([1,-1]) * spin_array[s]

        return Sm_array, spin_array, 0
    
    else:
        # temp arrays to store indices identifying the first and last vertices for each spin
        first_vertex_leg = -1 * np.ones(N, dtype=int)
        last_vertex_leg = -1 * np.ones(N, dtype=int)

        '''
        p = 0
        p_list = np.zeros(n, dtype=int)
        b_list = np.zeros(n, dtype=int)
        '''

        vertex_array = np.zeros(n, dtype=int)
        linked_list = np.zeros(num_legs_per_vertex*n, dtype=int)

        num_composite_leg_indices = num_legs_per_vertex**2

        # TODO: optimization: can likely populate p_list and b_list in diagonal updates and iterate over n<M here instead of M
        # For populating vertex_array and linked_list
        '''
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
        '''

        for p in range(len(p_list)):
            t = p_list[p]

            bond_num = b_list[p]
            bond_index = bond_num-1
            i_b = sites[bond_index,0]
            j_b = sites[bond_index,1]

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
            

        # Periodic due to trace - connect the first and last spins in expansion direction
        for i in range(N):
            if last_vertex_leg[i] != -1:
                linked_list[last_vertex_leg[i]] = first_vertex_leg[i]
                linked_list[first_vertex_leg[i]] = last_vertex_leg[i]

        skip_loop_update = False
        loop_size=0
        cumulative_loop_size = 0

        # Construct probabilistic loops
        for loop_num in range(N_l):

            # Cancel off-diagonal update if loop size exceeds allowed maximum
            if loop_size >= max_loop_size:
                skip_loop_update = True
                break

            cumulative_loop_size += loop_size # For tracking number of vertices visited in off-diagonal update
            loop_size = 0
            
            # pick a leg randomly to start
            j_0 = np.random.randint(num_legs_per_vertex*n)
            j = j_0

            while loop_size < max_loop_size:

                # Associated vertex and entrance leg
                p = j // num_legs_per_vertex
                l_e = j % num_legs_per_vertex
                vertex_type = vertex_array[p]
                
                # Index probability table to get exit leg given input leg and vertex
                composite_leg_index = num_legs_per_vertex*l_e
                row_index = num_composite_leg_indices*vertex_type + composite_leg_index

                out_leg_probs = off_diag_prob_table[row_index:row_index+4,b_list[p]-1]
                l_x = np.random.choice(num_legs_per_vertex,p=out_leg_probs)

                # Update vertex by flipping the input leg and output leg spins
                vertex_array[p] = new_vertex_map[vertex_type, composite_leg_index + l_x]
                
                # Continue along linked list to get next j, close loop if required
                j = num_legs_per_vertex*p + l_x
                if l_e != l_x:
                    loop_size += 1
                if j == j_0:
                    break
                else:
                    j = linked_list[j]
                    if j == j_0:
                        break
        
        # Update Sm_array and spin_array using the updated vertex array if loops of allowed sizes were successfully constructed
        if not skip_loop_update:
            # TODO: optimization: can iterate over p_list here instead of M
            for p in range(len(p_list)):
                t = p_list[p]
                b = b_list[p]
                Sm_array[t] = 2*b + operator_type[vertex_array[p]]
            
            free_spins = False
            free_spin_indices = []
            for s in range(N):
                # Flip disconnected spins with probability 1/2
                if first_vertex_leg[s] == -1:
                    spin_array[s] = np.random.choice([1,-1]) * spin_array[s]
                    free_spins = True
                    free_spin_indices.append(s)
                else:
                    # Update connected spins according to the first vertex acting on them - consistency in (1+1)D space ensured by construction
                    q = first_vertex_leg[s] // num_legs_per_vertex
                    l = first_vertex_leg[s] % num_legs_per_vertex
                    spin_array[s] = leg_spin[vertex_array[q]][l]

            if free_spins:
                for i in range(num_free_flips-1):
                    for s in free_spin_indices:
                        spin_array[s] = np.random.choice([1,-1]) * spin_array[s]

            return Sm_array, spin_array, cumulative_loop_size

        free_spin_indices = []
        free_spins = False
        for s in range(N):
            if first_vertex_leg[s] == -1:
                spin_array[s] = np.random.choice([1,-1]) * spin_array[s]
                free_spins = True
                free_spin_indices.append(s)

        for i in range(num_free_flips-1):
            for s in free_spin_indices:
                spin_array[s] = np.random.choice([1,-1]) * spin_array[s]

        return Sm_array, spin_array, cumulative_loop_size
        

def diagonal_updates(Sm_array, spin_array, sites, diagonal_norm, max_norm_diag_probs, diag_probs, M, n, beta, num_bonds):

    p_list = []
    b_list = []
    for t in range(M):
        if Sm_array[t] == 0:
            # propose adding a diagonal operator
            acceptance_prob = beta * diagonal_norm/(M-n)
            u = np.random.uniform(0.0,1.0)
            if u <= acceptance_prob:
                b_prime = np.random.choice(num_bonds, p=max_norm_diag_probs) # Gibbs sampling to determine which bond to potentially populate
                i_b = sites[b_prime,0]
                j_b = sites[b_prime,1] 
                vertex_type = vertex_map[(spin_array[i_b], spin_array[j_b], spin_array[i_b], spin_array[j_b])]
                u = np.random.uniform(0.0,1.0) # For accepting gibbs proposed update 
                if u <= diag_probs[vertex_type, b_prime]:
                    bond_num = b_prime + 1
                    Sm_array[t] = 2*bond_num
                    n += 1
                    p_list.append(t)
                    b_list.append(bond_num)
        else:
            if Sm_array[t] % 2 == 1:
                p_list.append(t)
                b = Sm_array[t] // 2
                b_list.append(b)
                # propagate spin array for off-diagonal operators
                bond_index = b-1
                i_b = sites[bond_index,0]
                j_b = sites[bond_index,1]
                spin_array[i_b] = -spin_array[i_b]
                spin_array[j_b] = -spin_array[j_b]
            else:
                # propose removing the diagonal operator
                acceptance_prob = (M-n+1) / (beta * diagonal_norm)
                u = np.random.uniform(0.0,1.0)
                if u <= acceptance_prob:
                    Sm_array[t] = 0
                    n = n - 1
                else:
                    p_list.append(t)
                    b = Sm_array[t] // 2
                    b_list.append(b)

    p_list = np.array(p_list)
    b_list = np.array(b_list)

    return n, Sm_array, spin_array, p_list, b_list

def simulate_XXZ(N, M, num_bonds, equilibration_steps, mc_steps, sites, beta, diag_probs, max_norm_diag_probs, diagonal_norm, off_diag_probs, starting_config, a=1.25):

    #spin_array = np.zeros(N, dtype=int)
    '''
    for i in range(N):
        spin_array[i] = np.random.choice([1,-1])
    '''
    spin_array = starting_config
    Sm_array = np.zeros(M, dtype=int)

    magnetization_array = np.zeros(mc_steps)

    # Equilibration step
    n=0
    n_array = np.zeros(equilibration_steps)
    cumulative_loop_size_list = np.zeros(equilibration_steps, dtype=int)
    for step in range(equilibration_steps):
        n_new, Sm_array, spin_array, p_list, b_list = diagonal_updates(Sm_array, spin_array, sites, diagonal_norm, max_norm_diag_probs, diag_probs, M, n, beta, num_bonds)

        # track n to increase M as necessary
        n_array[step] = n_new
        M_new = max(int(a * np.max(n_array[:step+1])),M)
        '''
        if M - n_new < n_new/3:
            M_new = max(int(n_new + n_new/3), M)
        '''
        Sm_array_new = np.zeros(M_new, dtype=int)
        Sm_array_new[:M] = Sm_array
        Sm_array = Sm_array_new
        M = M_new

        n = n_new
        
        max_loop_size = int(100*np.mean(n_array[:step+1]))
        N_l = 2*M
        Sm_array, spin_array, cumulative_loop_size = off_diagonal_updates(Sm_array, spin_array, M, N, N_l, n, max_loop_size, sites, off_diag_probs, N_l, p_list, b_list)
        cumulative_loop_size_list[step] = cumulative_loop_size # track cumulative loop size to determine N_l such that cumulative loop size is approximately 2M

    # Use the last 1/10 equilibration steps to compute N_l
    average_cumulative_loop_size = np.mean(cumulative_loop_size_list[-int(equilibration_steps/10):])
    N_l = max(int(np.round(2*M/average_cumulative_loop_size)),10)

    n_array = np.zeros(mc_steps, dtype=int)
    #spin_array = starting_config
    #Sm_array = np.zeros(M, dtype=int)
    #n=0
    for step in range(mc_steps):
        # diagonal updates first
        n, Sm_array, spin_array, p_list, b_list = diagonal_updates(Sm_array, spin_array, sites, diagonal_norm, max_norm_diag_probs, diag_probs, M, n, beta, num_bonds)
        # off-diagonal updates
        Sm_array, spin_array, cumulative_loop_size = off_diagonal_updates(Sm_array, spin_array, M, N, N_l, n, max_loop_size, sites, off_diag_probs, N_l, p_list, b_list)
        '''
        # pick random time-slice for recording structural properties
        time_slice = np.random.choice(M)
        for t in range(time_slice):
            if Sm_array[t] % 2 == 1:
                # propagate spin array for off-diagonal operators
                b = Sm_array[t] // 2
                bond_index = b-1
                i_b = sites[bond_index,0]
                j_b = sites[bond_index,1]
                spin_array[i_b] = -spin_array[i_b]
                spin_array[j_b] = -spin_array[j_b]
        '''
        magnetization_array = est.update_diagonal_estimators(spin_array, magnetization_array, step, N)
        '''
        for t in range(time_slice, M):
            if Sm_array[t] % 2 == 1:
                # propagate spin array for off-diagonal operators
                b = Sm_array[t] // 2
                bond_index = b-1
                i_b = sites[bond_index,0]
                j_b = sites[bond_index,1]
                spin_array[i_b] = -spin_array[i_b]
                spin_array[j_b] = -spin_array[j_b]
        '''

        n_array[step] = n

    energy_array = -n_array/beta

    '''
    energy_mean, energy_error = statistics_binning(energy_array) # use binning for errors etc., can likely improve this for more sensitive estimators
    magnetization_mean, magnetization_error = statistics_binning(magnetization_array)
    '''

    return energy_array, magnetization_array

