import numpy as np
from SSE_initializer import *
import math

def off_diagonal_updates(Sm_array, spin_array, M, N, N_l, n, max_loop_size, sites, off_diag_prob_table):

    # For a discussion of the algorithm, see Appendix A of: https://arxiv.org/pdf/cond-mat/0202316

    if n == 0:
        return Sm_array, spin_array, 0
    else:
        # temp arrays to store indices identifying the first and last vertices for each spin
        first_vertex_leg = -1 * np.ones(N, dtype=int)
        last_vertex_leg = -1 * np.ones(N, dtype=int)

        p = 0
        p_list = np.zeros(n, dtype=int)
        b_list = np.zeros(n, dtype=int)

        vertex_array = np.zeros(n, dtype=int)
        linked_list = np.zeros(num_legs_per_vertex*n, dtype=int)

        num_composite_leg_indices = num_legs_per_vertex**2

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

        skip_loop_update = False
        loop_size=0
        cumulative_loop_size = 0

        # Construct probabilistic loops
        for loop_num in range(N_l):

            # Cancel off-diagonal update if loop size exceeds allowed maximum
            if loop_size == max_loop_size:
                skip_loop_update = True
                break

            cumulative_loop_size += loop_size # For tracking number of vertices visited in off-diagonal update
            loop_size = 0
            
            # pick a pick randomly to start
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
            p=0
            for t in range(M):
                if Sm_array[t] != 0:
                    b = Sm_array[t]//2
                    Sm_array[t] = 2*b + operator_type[vertex_array[p]]
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
        
        return Sm_array, spin_array, cumulative_loop_size

def diagonal_updates(Sm_array, spin_array, sites, diagonal_norm, max_norm_diag_probs, diag_probs, M, n, beta, num_bonds):

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
                acceptance_prob = (M-n+1) / (beta * diagonal_norm)
                u = np.random.uniform(0.0,1.0)
                if u <= acceptance_prob:
                    Sm_array[t] = 0
                    n = n - 1

    return n, Sm_array, spin_array

def simulate_XXZ(N, M, num_bonds, equilibration_steps, mc_steps, sites, beta, diag_probs, max_norm_diag_probs, diagonal_norm, off_diag_probs, a=1.25):

    spin_array = np.ones(N, dtype=int)
    Sm_array = np.zeros(M, dtype=int)

    # Equilibration step
    n=0
    n_array = []
    cumulative_loop_size_list = []
    for step in range(equilibration_steps):
        n_new, Sm_array, spin_array = diagonal_updates(Sm_array, spin_array, sites, diagonal_norm, max_norm_diag_probs, diag_probs, M, n, beta, num_bonds)

        # track n to increase M as necessary
        n_array.append(n_new)
        M_new = max(int(a * np.max(n_array)),M)

        Sm_array_new = np.zeros(M_new, dtype=int)
        Sm_array_new[:M] = Sm_array
        Sm_array = Sm_array_new

        n = n_new
        M = M_new
        
        max_loop_size = 100*np.mean(n_array)

        Sm_array, spin_array, cumulative_loop_size = off_diagonal_updates(Sm_array, spin_array, M, N, 2*M, n, max_loop_size, sites, off_diag_probs)
        cumulative_loop_size_list.append(cumulative_loop_size) # track cumulative loop size to determine N_l such that cumulative loop size is approximately 2M

    # Use the last 1/10 equilibration steps to compute N_l 
    average_cumulative_loop_size = np.mean(cumulative_loop_size_list[-int(equilibration_steps/10):])
    N_l = max(int(np.round(2*M/average_cumulative_loop_size)),10)

    #print("Equilibration Finished")

    n_array = np.zeros(mc_steps)
    for step in range(mc_steps):
        # diagonal updates first
        n, Sm_array, spin_array = diagonal_updates(Sm_array, spin_array, sites, diagonal_norm, max_norm_diag_probs, diag_probs, M, n, beta, num_bonds)
        # off-diagonal updates
        Sm_array, spin_array, cumulative_loop_size = off_diagonal_updates(Sm_array, spin_array, M, N, N_l, n, max_loop_size, sites, off_diag_probs)

        n_array[step] = n

    energy_array = -n_array/beta

    energy_mean, energy_error = statistics_binning(energy_array) # use binning for errors etc., can likely improve this for more sensitive estimators

    return energy_array, energy_mean, energy_error

def statistics_binning(arr):
    # Average and standard error using the binning method
    workingNdim  = int(math.log(len(arr))/math.log(2))
    trunc = int(len(arr)-2**workingNdim)
    mean = np.mean(arr[trunc:])
    standardError = max_error_binning(arr[trunc:], workingNdim-6)
    return mean, standardError

def error_propagation(data):
    ndim   = len(data)
    error = np.std(data,ddof=0)/np.sqrt(ndim)
    return error

def max_error_binning(data, workingNdim):
    if(workingNdim<=1):
        raise Exception('Not enough points MC steps were used for the binning method, please increase the number of MC steps')
    error = np.zeros(workingNdim)
    i = 0
    error[0] = error_propagation(data)

    for i in range(1,workingNdim):
        ndim = int(len(data)/2)
        data1 = np.zeros(ndim)

        for j in range(ndim):
            data1[j] = 0.5*(data[2*j]+data[2*j+1])
        data = data1
        error[i] = error_propagation(data)
    return np.max(error)

