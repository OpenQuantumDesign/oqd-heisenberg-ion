#XXZ with delta
using LinearAlgebra
using CSV, DataFrames

#function for the distance 1/(r_ij ^alpha)
function r_distance(index_i, index_j, alpha)
    dist = abs(index_i - index_j)^alpha
    return (1 / dist)
end

function exact_diagonalization_theta(N, delta, h, J_X, J_Y, alpha, B, theta)

    println(delta)
    println("\n")

    Dim = 2^N
    J_Z = delta

    cutoff = 1e-3
    # initialise Hamiltonian
    Hamiltonian = zeros(ComplexF32, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1 #we're doing i!=j so here looping first for i<j (effectively half a sum). Then we multiply by 2 every Hamiltonian contribution to account for all pairs. 
            for j = i+1:N-1

                #if cutoff not wanted, comment the following lines

                #if var_ij < cutoff
                #    continue
                #end

                #S^Z part, similar to before 
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                Diagonal += -1.0 * r_distance(i, j, alpha) * J_Z * 0.25 * Spin1 * Spin2


                #S^X part 
                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j
                #Hamiltonian[Bra+1, Ket+1] += -1.0 * J_X * 0.25 * r_distance(i, j, alpha) * cos(theta)

                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                #S^Y part- not all terms have the same sign from the fact that S^Y=1/2i (S^+-S^-). I expanded S^Y*S^Y to check which terms would have an overall minus sign. 

                sign = (si == sj) ? 0 : 1
                Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * r_distance(i, j, alpha) * sign * exp((si-sj) * theta * 1im)


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

    #for row_index=1:8
    #    println(Hamiltonian[row_index,:])
    #end

    Diag = eigen(Hamiltonian)
    EigenVectors = Diag.vectors  #this gives the groundstate eigenvector
    EigenEnergies = Diag.values

    #index = findmin(EigenEnergies)

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

        #println("Magnetization along Z per site: ", mag_Z / N)

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

        #println("Magnetization along X per site: ", mag_X / N)

        mag_Z_all[ev] = mag_Z / N
        mag_X_all[ev] = real(mag_X) / N
    end

    for row_index=0:7
        if imag(EigenEnergies[row_index+1]) < 1e-12
            EigenEnergies[row_index+1] = real(EigenEnergies[row_index+1])
        else
            println("Large imaginary part in energy")
            println(EigenEnergies[row_index + 1])
        end
    end

    #EigenEnergies = real.(EigenEnergies)
    println(EigenEnergies)
    println("\n")

    results_df = DataFrame(
        Index=1:Dim,
        Eigenvalue=EigenEnergies,
        MagZ=mag_Z_all,
        MagX=mag_X_all,
    )

    CSV.write("Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(theta)_Heisenberg_OBC.csv", results_df)

end

function exact_diagonalization(N, delta, h, J_X, J_Y, alpha, B)

    println(delta)
    println("\n")

    Dim = 2^N
    J_Z = delta

    cutoff = 1e-3
    # initialise Hamiltonian
    Hamiltonian = zeros(Float32, Dim, Dim)

    for Ket = 0:Dim-1
        Diagonal = 0.0

        for i = 0:N-1 #we're doing i!=j so here looping first for i<j (effectively half a sum). Then we multiply by 2 every Hamiltonian contribution to account for all pairs. 
            for j = i+1:N-1

                #if cutoff not wanted, comment the following lines

                #if var_ij < cutoff
                #    continue
                #end

                #S^Z part, similar to before 
                Spin1 = 2 * ((Ket >> i) & 1) - 1
                Spin2 = 2 * ((Ket >> j) & 1) - 1
                Diagonal += -1.0 * r_distance(i, j, alpha) * J_Z * 0.25 * Spin1 * Spin2


                #S^X part 
                bit_i = 2^i
                bit_j = 2^j
                Bra = Ket ⊻ bit_i ⊻ bit_j
                #Hamiltonian[Bra+1, Ket+1] += -1.0 * J_X * 0.25 * r_distance(i, j, alpha) * cos(theta)

                si = (Ket >> i) & 1
                sj = (Ket >> j) & 1

                #S^Y part- not all terms have the same sign from the fact that S^Y=1/2i (S^+-S^-). I expanded S^Y*S^Y to check which terms would have an overall minus sign. 

                sign = (si == sj) ? 0 : 1
                Hamiltonian[Bra+1, Ket+1] += -0.5 * 1 * J_Y * r_distance(i, j, alpha) * sign


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

    #for row_index=1:8
    #    println(Hamiltonian[row_index,:])
    #end

    Diag = eigen(Hamiltonian)
    EigenVectors = Diag.vectors  #this gives the groundstate eigenvector
    EigenEnergies = Diag.values

    #index = findmin(EigenEnergies)

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

        #println("Magnetization along Z per site: ", mag_Z / N)

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

        #println("Magnetization along X per site: ", mag_X / N)

        mag_Z_all[ev] = mag_Z / N
        mag_X_all[ev] = real(mag_X) / N
    end

    for row_index=0:7
        if imag(EigenEnergies[row_index+1]) < 1e-12
            EigenEnergies[row_index+1] = real(EigenEnergies[row_index+1])
        else
            println("Large imaginary part in energy")
            println(EigenEnergies[row_index + 1])
        end
    end

    #EigenEnergies = real.(EigenEnergies)
    println(EigenEnergies)
    println("\n")

    results_df = DataFrame(
        Index=1:Dim,
        Eigenvalue=EigenEnergies,
        MagZ=mag_Z_all,
        MagX=mag_X_all,
    )

    CSV.write("Results/Exact_Diagonalization/ED_N_$(N)_Delta_$(delta)_h_$(h)_Jx_$(J_X)_Jy_$(J_Y)_alpha_$(alpha)_B_$(B)_theta_$(0.0)_Heisenberg_OBC.csv", results_df)

end

N = 3
J=1.0
J_X = J
J_Y = J
delta_list = [-11.0, -12.0, -13.0, -14.0, -15.0, -16.0, -17.0, -18.0, -19.0, -20.0, -3.0,-4.0,-5.0,-6.0,-7.0,-8.0,-9.0,-10.0,-2.0,-1.0, 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0]
#delta_list = [-17.0, -18.0, -19.0, -20.0]
#delta_list = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
#delta_list = [0.0]
#h_list = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
h_list = [0.0]
delta_list = [-1.0, 0.0, 1.0]
theta = 0.1
alpha = 5.0
B=0.0
for delta in delta_list
    for h in h_list
        exact_diagonalization_theta(N, delta, h, J_X, J_Y, alpha, B, theta)
        exact_diagonalization(N, delta, h, J_X, J_Y, alpha, B)
        exact_diagonalization_theta(N, delta, h, J_X, J_Y, alpha, B, -theta)
    end
end
