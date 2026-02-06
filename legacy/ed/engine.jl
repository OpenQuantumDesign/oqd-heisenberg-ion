#XXZ with delta
using LinearAlgebra
using CSV, DataFrames
using DelimitedFiles

#function for the distance 1/(r_ij ^alpha)
function r_distance_obc(index_i, index_j, alpha)
    dist = abs(index_i - index_j)^alpha
    return (1 / dist)
end

function r_distance_pbc(index_i, index_j, alpha, N)
    if index_i == index_j
        return 0.0
    end
    pbc = min(abs(index_i - index_j), N - abs(index_i - index_j))
    if alpha != 0.0
        dist = pbc^alpha
    else
        dist = 1.0
    end
    return (1 / dist)
end

function r_distance_NN_pbc(index_i, index_j, N)
    if index_i == index_j
        return 0
    end
    pbc = min(abs(index_i - index_j), N - abs(index_i - index_j))
    if (pbc > 1)
        return 0
    else
        return 1
    end
end

function r_distance_NN_obc(index_i, index_j)
    if index_i == index_j
        return 0
    end
    pbc = abs(index_i - index_j)
    if (pbc > 1)
        return 0
    else
        return 1
    end
end

function hamiltonian_pbc_long_range_twist_even_N(N, delta, J_Y, h, B, alpha, theta)

    Dim = 2^N
    J_Z = delta

    Hamiltonian = zeros(ComplexF64, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1 
            for r = 1:round(Int32, N/2)-1
                j = i + r
                if j > (N-1)
                    j = j - N
                end

                #S^Z part
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                one_over_r_b_alpha = r_distance_pbc(i, j, alpha, N)
                
                Diagonal += -1.0 * one_over_r_b_alpha * J_Z * 0.25 * Spin1 * Spin2

                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j
                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                #off-diagonal part

                sign = (si == sj) ? 0 : 1
                i_minus_j = min(abs(i - j), N - abs(i - j))
                Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * one_over_r_b_alpha * sign * exp((si-sj) * i_minus_j * theta * 1im)

            end
        end

        for i=0:round(Int64, N/2)-1
            r = round(Int64, N/2)
            j = (i + r)
            if j > (N-1)
                j = j - N
            end
            
            Spin1 = 2 * ((Ket >> i) & 1) - 1
            Spin2 = 2 * ((Ket >> j) & 1) - 1
            one_over_r_b_alpha = r_distance_pbc(i, j, alpha, N)
            Diagonal += -1.0 * one_over_r_b_alpha * J_Z * 0.25 * Spin1 * Spin2

            bit_i = 2^i
            bit_j = 2^j
            Bra = Ket ⊻ bit_i ⊻ bit_j
            si = (Ket >> i) & 1
            sj = (Ket >> j) & 1

            #off-diagonal part

            sign = (si == sj) ? 0 : 1
            i_minus_j = min(abs(i - j), N - abs(i - j))
            Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * one_over_r_b_alpha * sign * exp((si-sj) * i_minus_j * theta * 1im)
        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian 
end

function hamiltonian_pbc_long_range_twist_odd_N(N, delta, J_Y, h, B, alpha, theta)

    Dim = 2^N
    J_Z = delta

    Hamiltonian = zeros(ComplexF64, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1 
            for r = 1:round(Int32, (N-1)/2)
                j = i + r
                if j > (N-1)
                    j = j - N
                end

                #S^Z part
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                one_over_r_b_alpha = r_distance_pbc(i, j, alpha, N)
                
                Diagonal += -1.0 * one_over_r_b_alpha * J_Z * 0.25 * Spin1 * Spin2

                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j
                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                #off-diagonal part

                sign = (si == sj) ? 0 : 1
                i_minus_j = min(abs(i - j), N - abs(i - j))
                Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * one_over_r_b_alpha * sign * exp((si-sj) * i_minus_j * theta * 1im)

            end
        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian 
end

function hamiltonian_pbc_long_range(N, delta, J_Y, h, B, alpha)

    Dim = 2^N
    J_Z = delta

    Hamiltonian = zeros(Float64, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1
            for j = i+1:N-1

                #S^Z part
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                one_over_r_b_alpha = r_distance_pbc(i, j, alpha, N)
                Diagonal += -1.0 * one_over_r_b_alpha * J_Z * 0.25 * Spin1 * Spin2

                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j

                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                sign = (si == sj) ? 0 : 1
                Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * one_over_r_b_alpha * sign

            end

        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        #field along X
        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian 
end

function hamiltonian_obc_long_range_twist(N, delta, J_Y, h, B, alpha, theta)

    Dim = 2^N
    J_Z = delta

    Hamiltonian = zeros(ComplexF64, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1 
            for j = i+1:N-1

                #S^Z part, similar to before 
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                one_over_r_b_alpha = r_distance_obc(i, j, alpha)
                Diagonal += -1.0 * one_over_r_b_alpha * J_Z * 0.25 * Spin1 * Spin2

                #S^X part 
                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j

                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                sign = (si == sj) ? 0 : 1
                i_minus_j = min(abs(i - j), N - abs(i - j))
                Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * one_over_r_b_alpha * sign * exp((si-sj) * i_minus_j * theta * 1im)

            end

        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian
end

function hamiltonian_obc_long_range(N, delta, J_Y, h, B, alpha)

    Dim = 2^N
    J_Z = delta

    cutoff = 1e-3
    # initialise Hamiltonian
    Hamiltonian = zeros(Float64, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1 
            for j = i+1:N-1

                #S^Z part
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                one_over_r_b_alpha = r_distance_obc(i, j, alpha)
                Diagonal += -1.0 * one_over_r_b_alpha * J_Z * 0.25 * Spin1 * Spin2

                #S^X part 
                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j

                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                sign = (si == sj) ? 0 : 1
                Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * one_over_r_b_alpha * sign
            end

        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian 
end

function hamiltonian_obc_long_range_exp_J_ij(N, delta, exp_J_ij_file)

    Dim = 2^N
    J_Z = delta

    cutoff = 1e-3
    # initialise Hamiltonian
    Hamiltonian = zeros(Float64, Dim, Dim)

    J_ij_matrix = readdlm(exp_J_ij_file, ',', skipstart=1)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1 
            for j = i+1:N-1

                #S^Z part
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                #one_over_r_b_alpha = r_distance_obc(i, j, alpha)
                Diagonal += -1.0 * J_ij_matrix[i+1,j+1] * J_Z * 0.25 * Spin1 * Spin2

                #S^X part 
                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j

                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                sign = (si == sj) ? 0 : 1
                Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_ij_matrix[i+1,j+1] * sign
            end

        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian 
end

function hamiltonian_pbc_nearest_neighbour_twist(N, delta, J_Y, h, B, theta)

    Dim = 2^N
    J_Z = delta

    Hamiltonian = zeros(ComplexF64, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-2
            j = i+1

            #S^Z part, similar to before 
            Spin1 = 2 * ((Ket >> i) & 1) - 1
            Spin2 = 2 * ((Ket >> j) & 1) - 1
            int_strength = r_distance_NN_pbc(i, j, N)
            if (int_strength != 0)
                Diagonal += -1.0 * int_strength * J_Z * 0.25 * Spin1 * Spin2
            end

            bit_i = 2^i
            bit_j = 2^j
            Bra = Ket ⊻ bit_i ⊻ bit_j

            si = (Ket >> i) & 1
            sj = (Ket >> j) & 1

            sign = (si == sj) ? 0 : 1
            if (int_strength != 0)
                Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * int_strength * sign * exp((si-sj) * theta * 1im)
            end

        end

        if N != 2
            i = 0
            j = N-1

            #S^Z part, similar to before 
            Spin1 = 2 * ((Ket >> i) & 1) - 1
            Spin2 = 2 * ((Ket >> j) & 1) - 1
            int_strength = r_distance_NN_pbc(i, j, N)
            Diagonal += -1.0 * int_strength * J_Z * 0.25 * Spin1 * Spin2

            bit_i = 2^i
            bit_j = 2^j
            Bra = Ket ⊻ bit_i ⊻ bit_j

            si = (Ket >> i) & 1
            sj = (Ket >> j) & 1

            sign = (si == sj) ? 0 : 1
            Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * int_strength * sign * exp(-(si-sj) * theta * 1im)
        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian 
end

function hamiltonian_pbc_nearest_neighbour(N, delta, J_Y, h, B)

    Dim = 2^N
    J_Z = delta

    Hamiltonian = zeros(Float64, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1 
            for j = i+1:N-1

                #S^Z part
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                int_strength = r_distance_NN_pbc(i, j, N)
                if (int_strength != 0)
                    Diagonal += -1.0 * int_strength * J_Z * 0.25 * Spin1 * Spin2
                end

                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j

                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                sign = (si == sj) ? 0 : 1
                if (int_strength != 0)
                    Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * int_strength * sign
                end

            end

        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian
end

function hamiltonian_obc_nearest_neighbour_twist(N, delta, J_Y, h, B, theta)

    Dim = 2^N
    J_Z = delta

    Hamiltonian = zeros(ComplexF64, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1
            for j = i+1:N-1

                #S^Z part
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                int_strength = r_distance_NN_obc(i, j)
                if (int_strength != 0)
                    Diagonal += -1.0 * int_strength * J_Z * 0.25 * Spin1 * Spin2
                end

                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j

                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                sign = (si == sj) ? 0 : 1
                if (int_strength != 0)
                    Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * int_strength * sign * exp((si-sj) * theta * 1im)
                end

            end

        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian 
end

function hamiltonian_obc_nearest_neighbour(N, delta, J_Y, h, B)

    Dim = 2^N
    J_Z = delta

    Hamiltonian = zeros(Float64, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1
            for j = i+1:N-1

                #S^Z part
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                int_strength = r_distance_NN_obc(i, j)
                if (int_strength != 0)
                    Diagonal += -1.0 * int_strength * J_Z * 0.25 * Spin1 * Spin2
                end

                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j
                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                sign = (si == sj) ? 0 : 1
                if (int_strength != 0)
                    Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * int_strength * sign
                end

            end

        end

        #field along Z
        for k in 0:N-1
            spin_k = 2 * ((Ket >> k) & 1) - 1
            Diagonal += -0.5 * h * spin_k
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            bit = 2^SpinIndex   #The "label" of the bit to be flipped
            Bra = Ket ⊻ bit    #Binary XOR flips the bit
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian 

end

#=
function exact_diagonalization(N, delta, J_X, alpha, theta, hamiltonian_type, exp_J_ij_filepath)
    
    Dim = 2^N
    h = 0.0
    B = 0.0
    J_Y = J_X
    if hamiltonian_type == 1
        if N % 2 == 0
            Hamiltonian = hamiltonian_pbc_long_range_twist_even_N(N, delta, J_X, alpha, theta)
        else
            Hamiltonian = hamiltonian_pbc_long_range_twist_odd_N(N, delta, J_X, alpha, theta)
        end
        filename = "Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(round(theta, digits=3))_Heisenberg_PBC.csv"
    elseif hamiltonian_type == 2
        Hamiltonian = hamiltonian_obc_long_range_twist(N, delta, J_X, alpha, theta)
        filename = "Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(round(theta, digits=3))_Heisenberg_OBC.csv"
    elseif hamiltonian_type == 3
        if alpha != "inf"
            throw("theta should be inf with hamiltonian type 3")
        end
        Hamiltonian = hamiltonian_pbc_nearest_neighbour_twist(N, delta, J_X, theta)
        filename = "Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(round(theta, digits=3))_Heisenberg_PBC.csv"
    elseif hamiltonian_type == 4
        if alpha != "inf"
            throw("theta should be inf with hamiltonian type 4")
        end
        Hamiltonian = hamiltonian_obc_nearest_neighbour_twist(N, delta, J_X, theta)
        filename = "Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(round(theta, digits=3))_Heisenberg_OBC.csv"
    elseif hamiltonian_type == 5
        if theta != 0.0
            throw("theta should be 0.0 with hamiltonian type 5")
        end
        Hamiltonian = hamiltonian_pbc_long_range(N, delta, J_X, alpha)
        filename = "Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(round(theta, digits=3))_Heisenberg_PBC.csv"
    elseif hamiltonian_type == 6
        if theta != 0.0
            throw("theta should be 0.0 with hamiltonian type 6")
        end
        Hamiltonian = hamiltonian_obc_long_range(N, delta, J_X, alpha)
        filename = "Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(round(theta, digits=3))_Heisenberg_OBC.csv"
    elseif hamiltonian_type == 7
        if theta != 0.0
            throw("theta should be 0.0 with hamiltonian type 7")
        end
        if alpha != "inf"
            throw("theta should be inf with hamiltonian type 7")
        end
        Hamiltonian = hamiltonian_pbc_nearest_neighbour(N, delta, J_X)
        filename = "Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(round(theta, digits=3))_Heisenberg_PBC.csv"
    elseif hamiltonian_type == 8
        if theta != 0.0
            throw("theta should be 0.0 with hamiltonian type 8")
        end
        if alpha != "inf"
            throw("theta should be inf with hamiltonian type 8")
        end
        Hamiltonian = hamiltonian_obc_nearest_neighbour(N, delta, J_X)
        filename = "Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(round(theta, digits=3))_Heisenberg_OBC.csv"
    elseif hamiltonian_type == 9
        if (SubString(alpha, 1,3) != "exp") || (J_X != "exp")
            throw("alpha and J_X should be exp with hamiltonian_type 9")
        end
        Hamiltonian = hamiltonian_obc_long_range_exp_J_ij(N, delta, exp_J_ij_filepath)
        filename = "Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_exp_Jy_exp_alpha_$(alpha)_B_$(B)_theta_$(round(theta, digits=3))_Heisenberg_OBC.csv"
    else
        throw("Hamiltonian type not implemented")
    end

    Diag = eigen(Hamiltonian)
    EigenVectors = Diag.vectors  
    EigenEnergies = Diag.values

    #magnetization per site 
    #magnetization along Z
    mag_Z_all = zeros(Float32, Dim)
    mag_X_all = zeros(Float32, Dim)

    for ev in 1:Dim

        mag_Z = 0.0
        for Ket = 0:Dim-1
            spin_sum = 0.0
            for i = 0:N-1
                spin_sum += (2 * ((Ket >> i) & 1) - 1) * 0.5
            end
            mag_Z += spin_sum * abs2(EigenVectors[Ket+1, ev])
        end

        #magnetization along X
        mag_X = 0.0

        for Ket = 0:Dim-1
            amp_ket = EigenVectors[Ket+1, ev]
            for i = 0:N-1

                Bra = Ket ⊻ (2^i)
                amp_bra = EigenVectors[Bra+1, ev]
                mag_X += 0.5 * conj(amp_ket) * amp_bra
            end
        end

        mag_Z_all[ev] = mag_Z / N
        mag_X_all[ev] = real(mag_X) / N
    end

    results_df = DataFrame(
        Index=1:Dim,
        Eigenvalue=EigenEnergies,
        MagZ=mag_Z_all,
        MagX=mag_X_all,
    )

    CSV.write(filename, results_df)
    println(EigenEnergies[1])

end
=#