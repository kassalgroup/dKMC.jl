module dKMC_hopping_rates

#Loading the required packages.
include("../shared_functions/package_loading.jl")
include("constants.jl")


"""
    kdelta(a::Number,b::Number)

The Kronecker delta function.

# Arguments:
- `a,b`: Two numbers 

# Output:
- `1`: if a == b
- `0`: otherwise
"""
function kdelta(a::Number,b::Number)

    if (a == b)
        return 1
    else
        return 0
    end

end


"""
    lambda(a::Number, b::Number, c::Number, d::Number)

The lambda function in the sPTRE rate, a series of delta functions.

# Arguments:
- `a,b,c,d`: Four numbers.

# Output:
- `lambda`: The lamba value required in sPTRE rate.

"""
function lambda(a::Number, b::Number, c::Number, d::Number)

    lambda = kdelta(a,c) - kdelta(a,d) + kdelta(b,d) - kdelta(b,c)

    return lambda 

end


"""
    charge_transport_dKMC_rate(initial_state::Integer,final_state::Integer,
    transformed_coupling::Matrix{<:AbstractFloat},evals::Vector{<:AbstractFloat},evecs::Matrix{<:AbstractFloat},
    K_tot::Matrix{ComplexF64},E_step::AbstractFloat,E_limit::Integer,initial_state_relevant_sites::Vector{<:Integer},
    final_state_relevant_sites::Vector{<:Integer})

Calculates the dKMC charge transport delocalised hopping rate from one state to another. The rate is based on sPTRE, 
but only uses a subset of sites that are relevant to each state.

# Arguments:
- `initial_state`: Index of the initial state.
- `final_state`: Index of the final state.
- `transformed_coupling`: Matrix containing only polaron-transformed coupling terms.
- `evals`: Energy eigenvalues of the polaron trasformed system hamiltonian.
- `evecs`: Energy eigenvectors of the polaron trasformed system hamiltonian.
- `K_tot`: Precaclulated K values required for dKMC rate calculations.
- `E_step`: Energy step used in K value calculation.
- `E_limit`: Energy limit used in K value calculation.
- `initial_state_relevant_sites`: Sites relevant to the initial state, those that contain signicant contributions to the states amplitudes.
- `final_state_relevant_sites`: Sites relevant to the final state, those that contain signicant contributions to the states amplitudes.

# Output:
- `rate`: dKMC rate for delocalised hopping.

"""
function charge_transport_dKMC_rate(initial_state::Integer,final_state::Integer,transformed_coupling::Matrix{<:AbstractFloat},evals::Vector{<:AbstractFloat},evecs::Matrix{<:AbstractFloat},K_tot::Matrix{ComplexF64},E_step::AbstractFloat,E_limit::Integer,initial_state_relevant_sites::Vector{<:Integer},final_state_relevant_sites::Vector{<:Integer})

    #Start the damping rate at zero.
    gamma = 0

    #Check whether energy lies outside the range of precalculated K values.
    E = evals[initial_state] - evals[final_state]
    if abs(E) < E_limit

        #Find the K value (at each of the 5 lambda values) at an energy that is closest to the precomputed values stored in K_tot.
        K = real.(K_tot[Int(round((E+E_limit)/E_step))+1,:])

        #Find pairs of sites which have non zero couplings between them. 
        site_pairs_1 = [CartesianIndex(i,j) for i in initial_state_relevant_sites, j in final_state_relevant_sites if transformed_coupling[i,j]!=0]
        site_pairs_2 = [CartesianIndex(i,j) for i in final_state_relevant_sites, j in initial_state_relevant_sites if transformed_coupling[i,j]!=0]
        
        #Charge transfer terms.
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(i[1],i[2],i2[1],i2[2])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * K[cur_lamb+3]
                end
            end
        end
    end

    rate = (2/(constants.hbar^2))*real(gamma)

    return rate

end


"""
    exciton_transport_dKMC_rate(initial_state::Integer,final_state::Integer,
    transformed_coupling::Matrix{<:AbstractFloat},evals::Vector{<:AbstractFloat},evecs::Matrix{<:AbstractFloat},
    K_tot::Matrix{ComplexF64},E_step::AbstractFloat,E_limit::Integer,initial_state_relevant_sites::Vector{<:Integer},
    final_state_relevant_sites::Vector{<:Integer})

Calculates the dKMC exciton transport delocalised hopping rate from one state to another. The rate is based on sPTRE, 
but only uses a subset of sites that are relevant to each state.

# Arguments:
- `initial_state`: Index of the initial state.
- `final_state`: Index of the final state.
- `transformed_coupling`: Matrix containing only polaron-transformed coupling terms.
- `evals`: Energy eigenvalues of the polaron trasformed system hamiltonian.
- `evecs`: Energy eigenvectors of the polaron trasformed system hamiltonian.
- `K_tot`: Precaclulated K values required for dKMC rate calculations.
- `E_step`: Energy step used in K value calculation.
- `E_limit`: Energy limit used in K value calculation.
- `initial_state_relevant_sites`: Sites relevant to the initial state, those that contain signicant contributions to the states amplitudes.
- `final_state_relevant_sites`: Sites relevant to the final state, those that contain signicant contributions to the states amplitudes.

# Output:
- `rate`: dKMC rate for delocalised hopping.
    
"""
function exciton_transport_dKMC_rate(initial_state::Integer,final_state::Integer,transformed_coupling::Matrix{<:AbstractFloat},evals::Vector{<:AbstractFloat},evecs::Matrix{<:AbstractFloat},K_tot::Matrix{ComplexF64},E_step::AbstractFloat,E_limit::Integer,initial_state_relevant_sites::Vector{<:Integer},final_state_relevant_sites::Vector{<:Integer})

    #Start the damping rate at zero.
    gamma = 0

    #Check whether energy lies outside the range of precalculated K values.
    E = evals[initial_state] - evals[final_state]
    if abs(E) < E_limit
        
        #Find the K value (at each of the 5 lambda values) at an energy that is closest to the precomputed values stored in K_tot.
        K=real.(K_tot[Int(round((E+E_limit)/E_step))+1,:])

        #Exciton transfer terms.
        for n in final_state_relevant_sites, m in initial_state_relevant_sites
            gamma_1 = transformed_coupling[m,n] * conj(evecs[m,initial_state]) * evecs[n,final_state]
            for np in initial_state_relevant_sites, mp in final_state_relevant_sites
                cur_lamb = lambda(m,n,mp,np)
                if cur_lamb != 0
                    gamma += transformed_coupling[mp,np] * gamma_1 * conj(evecs[mp,final_state]) * evecs[np,initial_state] * K[cur_lamb+3]
                end
            end 
        end

    end 

    rate = (2/(constants.hbar^2))*real(gamma)

    return rate

end


"""
    charge_separation_dKMC_rate(initial_state::Integer,final_state::Integer,
    transformed_coupling::Matrix{<:AbstractFloat},evals::Vector{<:AbstractFloat},evecs::Matrix{<:AbstractFloat},
    K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:AbstractFloat},E_limits::Vector{<:Integer},
    initial_state_relevant_sites::Vector{<:Integer},final_state_relevant_sites::Vector{<:Integer},
    electron_index::Vector{<:Integer},hole_index::Vector{<:Integer},bath_index::Matrix{Integer})

Calculates the dKMC charge separation delocalised hopping rate from one state to another. The rate is based on sPTRE, 
but only uses a subset of sites that are relevant to each state.

# Arguments:
- `initial_state`: Index of the initial state.
- `final_state`: Index of the final state.
- `transformed_coupling`: Matrix containing only polaron-transformed coupling terms.
- `evals`: Energy eigenvalues of the polaron trasformed system hamiltonian.
- `evecs`: Energy eigenvectors of the polaron trasformed system hamiltonian.
- `K_tots`: Vector of precaclulated K_tot values for electrons and holes.
- `E_steps`: Vector of the energy step used for values that K_tot is calculated at for electrons and holes.
- `E_limits`: Vector of the energy limits used for values that K_tot is calculated at for electrons and holes.
- `initial_state_relevant_sites`: Sites relevant to the initial state, those that contain signicant contributions to the states amplitudes.
- `final_state_relevant_sites`: Sites relevant to the final state, those that contain signicant contributions to the states amplitudes.
- `electron_index`: Vector of indices to only the electron position of every site-pair.
- `hole_index`: Vector of indices to only the hole position of every site-pair.
- `bath_index`: Matrix of references to which particle type, and therefore which bath, corresponds to the coupling elements in Hamiltonian (1 for electrons, 2 for holes).

# Output:
- `rate`: dKMC rate for delocalised hopping.

"""
function charge_separation_dKMC_rate(initial_state::Integer,final_state::Integer,transformed_coupling::Matrix{<:AbstractFloat},evals::Vector{<:AbstractFloat},evecs::Matrix{<:AbstractFloat},K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:AbstractFloat},E_limits::Vector{<:Integer},initial_state_relevant_sites::Vector{<:Integer},final_state_relevant_sites::Vector{<:Integer},electron_index::Vector{<:Integer},hole_index::Vector{<:Integer},bath_index::Matrix{<:Integer})
    
    #Start the damping rate at zero.
    gamma = 0
    
    #Check whether energy lies outside the range of precalculated K values.
    E = evals[initial_state] - evals[final_state]
    if abs(E) < E_limits[1] && abs(E) < E_limits[2]

        #Find the K values (at each of the 5 lambda values) for electrons and holes at an energy that is closest to the precomputed values stored in K_tots.
        electron_K = real.(K_tots[1][Int(round((E+E_limits[1])/E_steps[1]))+1,:])
        hole_K = real.(K_tots[2][Int(round((E+E_limits[2])/E_steps[2]))+1,:])

        #Find pairs of sites which have non zero couplings between them. 
        pairs_1 = [CartesianIndex(i,j) for i in initial_state_relevant_sites, j in final_state_relevant_sites if transformed_coupling[i,j]!=0]
        pairs_2 = [CartesianIndex(i,j) for i in final_state_relevant_sites, j in initial_state_relevant_sites if transformed_coupling[i,j]!=0]
        
        #Electron transfer terms.		
        site_pairs_1 = [i for i in pairs_1 if bath_index[i]==1]
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==1]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(electron_index[i[1]],electron_index[i[2]],electron_index[i2[1]],electron_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * electron_K[cur_lamb+3]
                end 
            end 
        end 

        #Hole transfer terms.
        site_pairs_1 = [i for i in pairs_1 if bath_index[i]==2]
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==2]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(hole_index[i[1]],hole_index[i[2]],hole_index[i2[1]],hole_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * hole_K[cur_lamb+3]
                end 
            end 
        end
    end 

    rate = (2/(constants.hbar^2))*real(gamma)
    
    return rate

end


"""
    charge_generation_dKMC_rate(initial_state::Integer,final_state::Integer,
    transformed_coupling::Matrix{<:AbstractFloat},evals::Vector{<:AbstractFloat},evecs::Matrix{<:AbstractFloat},
    K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:AbstractFloat},E_limits::Vector{<:Integer},
    initial_state_relevant_sites::Vector{<:Integer},final_state_relevant_sites::Vector{<:Integer},
    electron_index::Vector{<:Integer},hole_index::Vector{<:Integer},exciton_index::Vector{<:Integer},
    bath_index::Matrix{Integer})

Calculates the dKMC charge generation delocalised hopping rate from one state to another. The rate is based on sPTRE, 
but only uses a subset of sites that are relevant to each state.

# Arguments:
- `initial_state`: Index of the initial state.
- `final_state`: Index of the final state.
- `transformed_coupling`: Matrix containing only polaron-transformed coupling terms.
- `evals`: Energy eigenvalues of the polaron trasformed system hamiltonian.
- `evecs`: Energy eigenvectors of the polaron trasformed system hamiltonian.
- `K_tots`: Vector of precaclulated K_tot values for electrons, holes, and excitons.
- `E_steps`: Vector of the energy step used for values that K_tot is calculated at for electrons, holes, and excitons.
- `E_limits`: Vector of the energy limits used for values that K_tot is calculated at for electrons, holes, and excitons.
- `initial_state_relevant_sites`: Sites relevant to the initial state, those that contain signicant contributions to the states amplitudes.
- `final_state_relevant_sites`: Sites relevant to the final state, those that contain signicant contributions to the states amplitudes.
- `electron_index`: Vector of indices to only the electron position of every site-pair.
- `hole_index`: Vector of indices to only the hole position of every site-pair.
- `exciton_index`: Vector of indices to only the exciton position of every site-pair.
- `bath_index`: Matrix of references to which particle type, and therefore which bath, corresponds to the coupling elements in Hamiltonian (1 for electrons, 2 for holes, 3 for excitons).

# Output:
- `rate`: dKMC rate for delocalised hopping.

"""
function charge_generation_dKMC_rate(initial_state::Integer,final_state::Integer,transformed_coupling::Matrix{<:AbstractFloat},evals::Vector{<:AbstractFloat},evecs::Matrix{<:AbstractFloat},K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:AbstractFloat},E_limits::Vector{<:Integer},initial_state_relevant_sites::Vector{<:Integer},final_state_relevant_sites::Vector{<:Integer},electron_index::Vector{<:Integer},hole_index::Vector{<:Integer},exciton_index::Vector{<:Integer},bath_index::Matrix{<:Integer})

    #Start the damping rate at zero.
    gamma = 0

    #Check whether energy lies outside the range of precalculated K values.
    E = evals[initial_state] - evals[final_state]
    if abs(E) < E_limits[1] && abs(E) < E_limits[2] && abs(E) < E_limits[3]

        #Find the K values (at each of the 5 lambda values) for electrons, holes, and excitons at an energy that is closest to the precomputed values stored in K_tots.
        electron_K = real.(K_tots[1][Int(round((E+E_limits[1])/E_steps[1]))+1,:])
        hole_K = real.(K_tots[2][Int(round((E+E_limits[2])/E_steps[2]))+1,:])
        exciton_K = real.(K_tots[3][Int(round((E+E_limits[3])/E_steps[3]))+1,:])
        electron_and_hole_K = real.(K_tots[4][Int(round((E+max(E_limits[1],E_limits[2]))/min(E_steps[1],E_steps[2])))+1,:])
        electron_and_exciton_K = real.(K_tots[5][Int(round((E+max(E_limits[1],E_limits[3]))/min(E_steps[1],E_steps[3])))+1,:])
        hole_and_exciton_K = real.(K_tots[6][Int(round((E+max(E_limits[2],E_limits[3]))/min(E_steps[2],E_steps[3])))+1,:])

        #Find pairs of sites which have non zero couplings between them. 
        pairs_1 = [CartesianIndex(i,j) for i in initial_state_relevant_sites, j in final_state_relevant_sites if transformed_coupling[i,j]!=0]
        pairs_2 = [CartesianIndex(i,j) for i in final_state_relevant_sites, j in initial_state_relevant_sites if transformed_coupling[i,j]!=0]
        
        #Electron transfer terms.		
        site_pairs_1 = [i for i in pairs_1 if bath_index[i]==1]
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==1]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(electron_index[i[1]],electron_index[i[2]],electron_index[i2[1]],electron_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * electron_K[cur_lamb+3]
                end 
            end 
        end 

        #Electron and hole mixed transfer terms.		
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==2]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(electron_index[i[1]],electron_index[i[2]],hole_index[i2[1]],hole_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * electron_and_hole_K[cur_lamb+3]
                end 
            end 
        end 

        #Electron and exciton mixed transfer terms.		
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==3]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(electron_index[i[1]],electron_index[i[2]],exciton_index[i2[1]],exciton_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * electron_and_exciton_K[cur_lamb+3]
                end 
            end 
        end 


        #Hole transfer terms.
        site_pairs_1 = [i for i in pairs_1 if bath_index[i]==2]
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==2]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(hole_index[i[1]],hole_index[i[2]],hole_index[i2[1]],hole_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * hole_K[cur_lamb+3]
                end 
            end 
        end 

        #Hole and electron mixed transfer terms.		
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==1]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(hole_index[i[1]],hole_index[i[2]],electron_index[i2[1]],electron_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * electron_and_hole_K[cur_lamb+3]
                end 
            end 
        end 

        #Hole and exciton mixed transfer terms.		
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==3]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(hole_index[i[1]],hole_index[i[2]],exciton_index[i2[1]],exciton_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * hole_and_exciton_K[cur_lamb+3]
                end 
            end 
        end 


        #Exciton transfer terms.
        site_pairs_1 = [i for i in pairs_1 if bath_index[i]==3]
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==3]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(exciton_index[i[1]],exciton_index[i[2]],exciton_index[i2[1]],exciton_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * exciton_K[cur_lamb+3]
                end 
            end 
        end 

        #Exciton and electron mixed transfer terms.		
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==1]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(exciton_index[i[1]],exciton_index[i[2]],electron_index[i2[1]],electron_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * electron_and_exciton_K[cur_lamb+3]
                end 
            end 
        end 

        #Exciton and hole mixed transfer terms.		
        site_pairs_2 = [i for i in pairs_2 if bath_index[i]==2]
        for i in site_pairs_1
            gamma_1 = transformed_coupling[i] * conj(evecs[i[1],initial_state]) * evecs[i[2],final_state]
            for i2 in site_pairs_2
                cur_lamb = lambda(exciton_index[i[1]],exciton_index[i[2]],hole_index[i2[1]],hole_index[i2[2]])
                if cur_lamb != 0
                    gamma += transformed_coupling[i2] * gamma_1 * conj(evecs[i2[1],final_state]) * evecs[i2[2],initial_state] * hole_and_exciton_K[cur_lamb+3]
                end 
            end 
        end

    end 

    rate = (2/(constants.hbar^2))*real(gamma)

    return rate

end

end
