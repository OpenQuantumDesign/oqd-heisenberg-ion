if __name__=="__main__":
    
    #Delta_list = [-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
    #Delta_list = [-3.0,-4.0,-5.0,-6.0,-7.0,-8.0,-9.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
    Delta_list = [-20.0,-19.0,-18.0,-17.0,-11.0, -12.0, -13.0, -14.0, -15.0, -16.0, -3.0,-4.0,-5.0,-6.0,-7.0,-8.0,-9.0,-10.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
    #Delta_list = [-3.0]
    Delta_list = [-1.0, 0.0, 1.0]

    h_list = [0.0]
    N_1 = 20
    N_2 = 0
    alpha = 1.0
    ksi = 0.0
    J = 1.0
    loop_update_type = "directed_loops"
    lattice_type = "1d"
    dist_dep_gamma = False
    alpha_list=[4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0]

    '''
    for alpha in alpha_list:
        for Delta in Delta_list:
            gamma = 0.1
            for h in h_list:
                write_prob_tables(gamma, Delta, h, N_1, N_2, lattice_type, alpha, ksi, J, loop_update_type, dist_dep_gamma, 2, 1, 0.0)
    '''

    N_list = [3]
    #N_list=[11]
    N_2 = 0
    J = 1.0
    #J = None
    #alpha_list = [2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00]
    alpha_list = [10.0]
    #alpha_list = [None]
    lattice_type = "1d"
    boundary = 0
    J_ij_type = 0
    #mu_list = [1.41]

    for N_1 in N_list:
        for alpha in alpha_list:
            #write_prob_tables_deterministic(N_1, N_2, lattice_type, alpha, J, 1, boundary)
            #write_prob_tables_deterministic(N_1, N_2, lattice_type, alpha, J, -1, boundary)
            write_prob_tables_deterministic(N_1, N_2, lattice_type, J_ij_type, None, alpha, J, 0, boundary, None)

    '''
    for N_1 in N_list:
        for mu in mu_list:
            J_ij_file_path = "/Users/shaeermoeed/Github/oqd-trical/Experimental_J_ij_mu_{}.csv".format(mu)
            write_prob_tables_deterministic(N_1, N_2, lattice_type, J_ij_type, J_ij_file_path, None, None, 0, boundary, mu)
    '''