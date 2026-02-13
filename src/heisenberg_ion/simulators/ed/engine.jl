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
    return -0.5 * mask * exp((si - sj) * r * theta * 1im)
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
    r_max = even_num_spins ? round(Int32, N / 2) - 1 : round(Int32, (N - 1) / 2)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1
            for r = 1:r_max
                j = (i + r) % N

                #S^Z part
                Diagonal += J * Delta * J_ij[i+1, j+1] * get_zz_term(Ket, i, j)

                #S^X part
                Bra = get_non_zero_row_index(Ket, i, j)
                Hamiltonian[Bra+1, Ket+1] += J * J_ij[i+1, j+1] * get_off_diag_term(Ket, i, j, r, theta)
            end
        end

        if even_num_spins
            for i = 0:round(Int64, N / 2)-1
                r = round(Int64, N / 2)
                j = (i + r) % N

                #S^Z part
                Diagonal += J * Delta * J_ij[i+1, j+1] * get_zz_term(Ket, i, j)

                #S^X part
                Bra = get_non_zero_row_index(Ket, i, j)
                Hamiltonian[Bra+1, Ket+1] += J * J_ij[i+1, j+1] * get_off_diag_term(Ket, i, j, r, theta)
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
            j = i + 1

            #S^Z part
            Diagonal += J * Delta * get_zz_term(Ket, i, j)

            Bra = get_non_zero_row_index(Ket, i, j)
            si = (Ket >> i) & 1
            sj = (Ket >> j) & 1
            mask = (si == sj) ? 0 : 1
            Hamiltonian[Bra+1, Ket+1] += J * get_off_diag_term(Ket, i, j, 1, theta)
        end

        if boundary == "periodic" && N != 2
            i = N - 1
            j = 0

            #S^Z part
            Diagonal += J * Delta * get_zz_term(Ket, i, j)

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
        Hamiltonian = zeros(Float64, Dim, Dim)
    else
        Hamiltonian = zeros(ComplexF64, Dim, Dim)
    end

    return Hamiltonian
end

function diagonalize_hamiltonian(; kwargs...)

    N = parse(Int, kwargs[:N])
    theta = parse(Float64, kwargs[:theta])
    Delta = parse(Float64, kwargs[:Delta])
    h = parse(Float64, kwargs[:h])
    B = parse(Float64, kwargs[:B])
    J = parse(Float64, kwargs[:J])
    interaction_name = kwargs[:interaction_range]

    Hamiltonian = initialize_hamiltonian(N, theta)

    if interaction_name == "long_range"

        J_ij_file = kwargs[:J_ij_file]
        J_ij = readdlm(J_ij_file, ',', skipstart=1)

        ed_parameters = (output_folder=kwargs[:run_folder], N=N, theta=theta, Delta=Delta, h=h, B=B, J=J,
            interaction_name=interaction_name, J_ij_file=J_ij_file)

        Hamiltonian = hamiltonian_long_range(N, Delta, h, B, J_ij, J, Hamiltonian, theta)

    elseif interaction_name == "nearest_neighbor"

        boundary = kwargs[:boundary]

        ed_parameters = (output_folder=kwargs[:run_folder], N=N, theta=theta, Delta=Delta, h=h, B=B, J=J,
            interaction_name=interaction_name, boundary=boundary)

        Hamiltonian = hamiltonian_nearest_neighbour(N, Delta, h, B, J, Hamiltonian, boundary, theta)

    else
        throw("interaction_type can be long_range or nearest_neighbor")
    end

    Diag = eigen(Hamiltonian)
    evecs = Diag.vectors
    evals = Diag.values

    return evals, evecs
end


function write_simulation_specs(parameter_set)

    simulation_spec_file_path = parameter_set.output_folder + "/simulation_specs.txt"
    open(simulation_spec_file_path) do f
        for (key, val) in pairs(parameter_set)
            println(key + "\t" + val + "\n")
        end
    end
end


function write_results(run_folder, evals, evecs)

    ed_out_dir = run_folder * "/ed_output"
    mkdir(ed_out_dir)

    evals_df = DataFrame(evals=evals)
    CSV.write(ed_out_dir * "/energies.csv", evals_df)

    column_names = [Symbol("|E_$i>") for i in 0:length(evals)-1]
    evecs_df = DataFrame(evecs, column_names)
    CSV.write(ed_out_dir * "/eigenvectors.csv", evecs_df)

end


function extract_inputs(input_file_path)

    input_parameters = Dict{Symbol,Vector}()

    num_parameter_sets = 0
    simulation_folder = 0
    simulation_folder_found = false
    open(input_file_path) do f
        for line in eachline(f)
            key_val = split(line, "\t")
            key = key_val[1]
            val = key_val[2]
            if key == "simulation_folder"
                simulation_folder_found = true
                simulation_folder = val
            else
                val_entries = split(val, "\t")
                num_parameter_sets = length(val_entries)
                input_parameters[Symbol(key)] = val_entries
            end
        end
    end

    if !simulation_folder_found
        throw("Unable to find simulation_folder in input file\n")
    end
    input_parameters[Symbol("simulation_folder")] = fill(simulation_folder, num_parameter_sets)

    parameter_sets = []
    for i = 1:num_parameter_sets
        single_parameter_set = Dict{Symbol,String}()
        for (key, val_list) in input_parameters
            single_parameter_set[key] = val_list[i]
        end
        push!(parameter_sets, single_parameter_set)
    end

    return parameter_sets

end


function execute_simulations(parameter_sets)

    num_parameter_sets = length(parameter_sets)

    for i = 1:num_parameter_sets
        kwargs = parameter_sets[i]
        evals, evecs = diagonalize_hamiltonian(; kwargs...)
        write_results(kwargs[:run_folder], evals, evecs)
    end
end


function main()

    input_file_path = ARGS[1]
    #input_file_path = "/Users/shaeermoeed/Github/Heisenberg_Ion/tests/integration/results/2026_02_05_20_07_40/inputs.txt"
    parameter_sets = extract_inputs(input_file_path)
    execute_simulations(parameter_sets)

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end