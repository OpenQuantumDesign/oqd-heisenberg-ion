#XXZ with delta
using LinearAlgebra
using CSV, DataFrames
using DelimitedFiles


function get_zz_term(Ket, i, j)

    Spin1 = 2 * ((Ket >> i) & 1) - 1
    Spin2 = 2 * ((Ket >> j) & 1) - 1

    return -0.25 * Spin1 * Spin2
end


function get_non_zero_row_index(Ket, i, j)

    bit_i = 2^i
    bit_j = 2^j
    return Ket ⊻ bit_i ⊻ bit_j
end


function get_off_diag_term(Ket, i, j, r, theta)

    si = (Ket >> i) & 1
    sj = (Ket >> j) & 1
    mask = (si == sj) ? 0 : 1
    return -0.5 * mask * exp((si-sj) * r * theta * 1im)
end


function get_longitudinal_field_term(Ket, k)

    spin_k = 2 * ((Ket >> k) & 1) - 1
    return -0.5 * spin_k
end


function get_transverse_field_index(Ket, SpinIndex)

    bit = 2^SpinIndex   # The "label" of the bit to be flipped
    return Ket ⊻ bit    # Binary XOR flips the bit
end


function hamiltonian_long_range(N, Delta, h, B, J_ij, J, Hamiltonian, theta)

    Dim = 2^N

    even_num_spins = (N % 2 == 0)
    r_max = even_num_spins ? round(Int32, N/2)-1 : round(Int32, (N-1)/2)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1
            for r = 1:r_max
                j = (i + r) % N

                #S^Z part
                Diagonal += Delta * J_ij[i+1,j+1] * get_zz_term(Ket, i, j)

                #S^X part
                Bra = get_non_zero_row_index(Ket, i, j)
                Hamiltonian[Bra+1, Ket+1] += J * J_ij[i+1,j+1] * get_off_diag_term(Ket, i, j, r, theta)
            end
        end
        
        if even_num_spins
            for i=0:round(Int64, N/2)-1
                r = round(Int64, N/2)
                j = (i + r) % N
                
                #S^Z part
                Diagonal += Delta * J_ij[i+1,j+1] * get_zz_term(Ket, i, j)

                #S^X part
                Bra = get_non_zero_row_index(Ket, i, j)
                Hamiltonian[Bra+1, Ket+1] += J * J_ij[i+1,j+1] * get_off_diag_term(Ket, i, j, r, theta)
            end
        end

        #field along Z
        for k in 0:N-1
            Diagonal += h * get_longitudinal_field_term(Ket, k)
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            Bra = get_transverse_field_index(Ket, SpinIndex)
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end
    end

    return Hamiltonian
end

function hamiltonian_nearest_neighbour(N, Delta, h, B, J, Hamiltonian, boundary, theta)

    Dim = 2^N

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-2
            j=i+1

            #S^Z part
            Diagonal += Delta * get_zz_term(Ket, i, j)

            Bra = get_non_zero_row_index(Ket, i, j)
            si = (Ket >> i) & 1
            sj = (Ket >> j) & 1
            mask = (si == sj) ? 0 : 1
            Hamiltonian[Bra+1, Ket+1] += J * get_off_diag_term(Ket, i, j, 1, theta)
        end

        if boundary == "periodic" && N != 2
            i = N-1
            j = 0

            #S^Z part
            Diagonal += Delta * get_zz_term(Ket, i, j)

            Bra = get_non_zero_row_index(Ket, i, j)
            si = (Ket >> i) & 1
            sj = (Ket >> j) & 1
            mask = (si == sj) ? 0 : 1
            Hamiltonian[Bra+1, Ket+1] += J * get_off_diag_term(Ket, i, j, 1, theta)
        end

        #field along Z
        for k in 0:N-1
            Diagonal += h * get_longitudinal_field_term(Ket, k)
        end

        Hamiltonian[Ket+1, Ket+1] = Diagonal

        for SpinIndex = 0:N-1
            Bra = get_transverse_field_index(Ket, SpinIndex)
            Hamiltonian[Bra+1, Ket+1] += -0.5 * B
        end

    end

    return Hamiltonian

end

function initialize_hamiltonian(N, theta)

    Dim = 2^N
    if theta == 0
        Hamiltonian = zeros(FloatF64, Dim, Dim)
    else
        Hamiltonian = zeros(ComplexF64, Dim, Dim)
    end

    return Hamiltonian
end

function diagonalize_hamiltonian(; N::Int, Delta::Float64, h::Float64, B::Float64, J::Float64, 
    theta::Int, interaction_type::String, kwargs...)

    Hamiltonian = initialize_hamiltonian(N, theta)

    if interaction_type == "long_range"

        J_ij_file = kwargs.J_ij_file
        J_ij = readdlm(J_ij_file, ',', skipstart=1)

        Hamiltonian = hamiltonian_long_range(N, Delta, h, B, J_ij, J, Hamiltonian, theta)

    elseif interaction_type == "nearest_neighbors"

        boundary = kwargs.boundary
        Hamiltonian = hamiltonian_nearest_neighbour(N, Delta, h, B, J, Hamiltonian, boundary, theta)

    else
        throw("interaction_type can be long_range or nearest_neighbors")
    end
    
    Diag = eigen(Hamiltonian)
    evecs = Diag.vectors
    evals = Diag.values

    return evals, evecs
end