import numpy as np

def geometry(N_1, N_2, geometry_type, boundary_conditions):
    
    if geometry_type == "1d":
        if N_2 == 0:
            if boundary_conditions == 0:
                sites, distances, geometry_table, N, num_bonds = geometry_1d_OBC(N_1)
            elif boundary_conditions == 1: 
                sites, distances, geometry_table, N, num_bonds = geometry_1d_PBC(N_1)
            else: 
                raise Exception("Boundary conditions for 1d need to be either 0 (OBC) or 1 (PBC).")
        else: 
            raise Exception("Geometry type is 1d but N_2 is non-zero.")
    elif geometry_type == "2d_triangular":
        if boundary_conditions == 0:
            sites, distances, geometry_table, N, num_bonds = geometry_2d_triangular_OBC(N_1, N_2)
        else:
            raise Exception("Boundary conditions for 2d need to be 0 (OBC).")
    else:
        raise Exception("Geometry type: {} is not implemented. Valid possible choices are: {} and {}.".format(geometry_type, "1d", "2d_triangular"))

    return sites, distances, geometry_table, N, num_bonds

def geometry_1d_OBC(N):

    num_bonds = int(N*(N-1)/2)
    sites = np.zeros((num_bonds,2), dtype=int)
    geometry_table = np.zeros((num_bonds,3))
    distances = np.zeros(num_bonds)
    b = 0
    for i in range(N):
        for j in range(i+1,N):
            sites[b,0] = i
            sites[b,1] = j
            distances[b] = j-i
            geometry_table[b,0] = i
            geometry_table[b,1] = j
            geometry_table[b,2] = distances[b]
            b += 1

    return sites, distances, geometry_table, N, num_bonds

def geometry_1d_PBC(N):

    num_bonds = int(N*(N-1)/2)
    sites = np.zeros((num_bonds,2), dtype=int)
    geometry_table = np.zeros((num_bonds,3))
    distances = np.zeros(num_bonds)
    b = 0
    for i in range(N):
        for j in range(i+1,N):
            sites[b,0] = i
            sites[b,1] = j
            geometry_table[b,0] = i
            geometry_table[b,1] = j
            if (j-i) <= N - (j-i):
                distances[b] = j-i
                geometry_table[b,2] = distances[b]
            else:
                distances[b] = N - (j-i)
                geometry_table[b,2] = -distances[b]
            b += 1

    return sites, distances, geometry_table, N, num_bonds

def geometry_2d_triangular_OBC(N_1, N_2):

    N = N_1 * N_2
    num_bonds = int(N*(N-1)/2)
    sites = np.zeros((num_bonds,2), dtype=int)
    distances = np.zeros(num_bonds)
    geometry_table = np.zeros(num_bonds,3)
    b = 0
    a_1 = np.array([1.0,0.0])
    a_2 = np.array([-0.5, np.sqrt(3.0)/2.0])
    for i1 in range(N_1):
        for i2 in range(N_2):
            i = i1*N_2 + i2
            for j1 in range(N_1):
                for j2 in range(N_2):
                    j = j1*N_2 + j2
                    if j > i:
                        sites[b,0] = i
                        sites[b,1] = j
                        distances[b] = np.norm((j1-i1)*a_1 + (j2-i2)*a_2)
                        geometry_table[b,0] = i
                        geometry_table[b,1] = j
                        geometry_table[b,2] = distances[b]
                        b += 1

    return sites, distances, geometry_table, N, num_bonds

# TODO: Implement 2d triangular PBC geometry