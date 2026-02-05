project_root = joinpath(@__DIR__, "..", "..")
include(joinpath(project_root, "legacy/ed", "engine.jl"))
include(joinpath(project_root, "src/heisenberg_ion/simulators/ed", "engine.jl"))

function check_hamiltonian_entries_obc_nearest_neighbors(N, Delta, J, h, B)

    ham_matrix = zeros(Float64, 2^N, 2^N)
    ham_1 = hamiltonian_obc_nearest_neighbour(N, Delta, J, h, B)
    ham_2 = hamiltonian_nearest_neighbour(N, Delta, h, B, J, ham_matrix, "open", 0)

    if !isapprox(ham_1, ham_2; rtol=1e-12, atol=1e-12)
        println("N = ", N)
        println("Delta = ", Delta)
        println("J = ", J)
        println("h = ", h)
        println("B = ", B)

        println(ham_1)
        println(ham_2)

        throw("Hamiltonians do not match")
    end

end

function check_hamiltonian_entries_pbc_nearest_neighbors(N, Delta, J, h, B)

    ham_matrix = zeros(Float64, 2^N, 2^N)
    ham_1 = hamiltonian_pbc_nearest_neighbour(N, Delta, J, h, B)
    ham_2 = hamiltonian_nearest_neighbour(N, Delta, h, B, J, ham_matrix, "periodic", 0)

    if !isapprox(ham_1, ham_2; rtol=1e-12, atol=1e-12)
        println("N = ", N)
        println("Delta = ", Delta)
        println("J = ", J)
        println("h = ", h)
        println("B = ", B)
        println(ham_1)
        println(ham_2)

        throw("Hamiltonians do not match")
    end

end

function check_hamiltonian_entries_pbc_nearest_neighbor_twist(N, Delta, J, h, B, theta)

    ham_matrix = zeros(ComplexF64, 2^N, 2^N)
    ham_1 = hamiltonian_pbc_nearest_neighbour_twist(N, Delta, J, h, B, theta)
    ham_2 = hamiltonian_nearest_neighbour(N, Delta, h, B, J, ham_matrix, "periodic", theta)

    if !isapprox(ham_1, ham_2; rtol=1e-12, atol=1e-12)
        println("N = ", N)
        println("Delta = ", Delta)
        println("J = ", J)
        println("h = ", h)
        println("B = ", B)
        println("theta = ", theta)
        println(ham_1)
        println(ham_2)

        throw("Hamiltonians do not match")
    end

end

function check_hamiltonian_entries_obc_long_range(N, Delta, J, h, B, alpha)

    ham_matrix = zeros(Float64, 2^N, 2^N)
    ham_1 = hamiltonian_obc_long_range(N, Delta, J, h, B, alpha)
    J_ij = zeros(Float64, N, N)
    for i=1:N
        for j = 1:N
            if i != j
                J_ij[i,j] = 1.0/((abs(i - j))^alpha)
            end
        end
    end
    ham_2 = hamiltonian_long_range(N, Delta, h, B, J_ij, J, ham_matrix, 0)

    if !isapprox(ham_1, ham_2; rtol=1e-12, atol=1e-12)
        println("N = ", N)
        println("Delta = ", Delta)
        println("J = ", J)
        println("h = ", h)
        println("B = ", B)
        println("alpha = ", alpha)

        println(ham_1)
        println(ham_2)

        throw("Hamiltonians do not match")
    end
end

function check_hamiltonian_entries_pbc_long_range(N, Delta, J, h, B, alpha)

    ham_matrix = zeros(Float64, 2^N, 2^N)
    ham_1 = hamiltonian_pbc_long_range(N, Delta, J, h, B, alpha)
    J_ij = zeros(Float64, N, N)
    for i=1:N
        for j = 1:N
            if i != j
                r = min(abs(i - j), N - abs(i - j))
                J_ij[i,j] = 1.0/((r)^alpha)
            end
        end
    end
    ham_2 = hamiltonian_long_range(N, Delta, h, B, J_ij, J, ham_matrix, 0)

    if !isapprox(ham_1, ham_2; rtol=1e-12, atol=1e-12)
        println("N = ", N)
        println("Delta = ", Delta)
        println("J = ", J)
        println("h = ", h)
        println("B = ", B)
        println("alpha = ", alpha)

        println(ham_1)
        println(ham_2)

        throw("Hamiltonians do not match")
    end
end

function check_hamiltonian_entries_pbc_long_range(N, Delta, J, h, B, alpha)

    ham_matrix = zeros(Float64, 2^N, 2^N)
    ham_1 = hamiltonian_pbc_long_range(N, Delta, J, h, B, alpha)
    J_ij = zeros(Float64, N, N)
    for i=1:N
        for j = 1:N
            if i != j
                r = min(abs(i - j), N - abs(i - j))
                J_ij[i,j] = 1.0/((r)^alpha)
            end
        end
    end
    ham_2 = hamiltonian_long_range(N, Delta, h, B, J_ij, J, ham_matrix, 0)

    if !isapprox(ham_1, ham_2; rtol=1e-12, atol=1e-12)
        println("N = ", N)
        println("Delta = ", Delta)
        println("J = ", J)
        println("h = ", h)
        println("B = ", B)
        println("alpha = ", alpha)

        println(ham_1)
        println(ham_2)

        throw("Hamiltonians do not match")
    end
end

function check_hamiltonian_entries_pbc_long_range_twist(N, Delta, J, h, B, alpha, theta)

    ham_matrix = zeros(ComplexF64, 2^N, 2^N)
    if (N % 2 == 0)
        ham_1 = hamiltonian_pbc_long_range_twist_even_N(N, Delta, J, h, B, alpha, theta)
    else
        ham_1 = hamiltonian_pbc_long_range_twist_odd_N(N, Delta, J, h, B, alpha, theta)
    end
    J_ij = zeros(Float64, N, N)
    for i=1:N
        for j = 1:N
            if i != j
                r = min(abs(i - j), N - abs(i - j))
                J_ij[i,j] = 1.0/((r)^alpha)
            end
        end
    end
    ham_2 = hamiltonian_long_range(N, Delta, h, B, J_ij, J, ham_matrix, theta)

    if !isapprox(ham_1, ham_2; rtol=1e-12, atol=1e-12)
        println("N = ", N)
        println("Delta = ", Delta)
        println("J = ", J)
        println("h = ", h)
        println("B = ", B)
        println("alpha = ", alpha)

        println(ham_1)
        println(ham_2)

        throw("Hamiltonians do not match")
    end
end

N_list = [2,3,4,5]
Delta_list = [0.0, 1.0, 2.0, -1.0, -2.0]
J_list = [0.0, 1.0, 2.0, -1.0, -2.0]
h_list = [0.0, 1.0, 2.0, -1.0, -2.0]
B_list = [0.0, 1.0, 2.0, -1.0, -2.0]
theta_list = [0.0, 1.0, 2.0, -1.0, -2.0]
alpha_list = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
for N in N_list
    for Delta in Delta_list
        for J in J_list
            for h in h_list
                for B in B_list
                    for theta in theta_list
                        check_hamiltonian_entries_pbc_nearest_neighbor_twist(N, Delta, J, h, B, theta)
                        for alpha in alpha_list
                            check_hamiltonian_entries_pbc_long_range_twist(N, Delta, J, h, B, alpha, theta)
                        end
                    end
                    check_hamiltonian_entries_obc_nearest_neighbors(N, Delta, J, h, B)
                    check_hamiltonian_entries_pbc_nearest_neighbors(N, Delta, J, h, B)
                    for alpha in alpha_list
                        check_hamiltonian_entries_obc_long_range(N, Delta, J, h, B, alpha)
                        check_hamiltonian_entries_pbc_long_range(N, Delta, J, h, B, alpha)
                    end
                end
            end
        end
    end
end
println("All Tests Passed")