module approximating_radii_functions

#Loading the required packages.
include("../shared_functions/package_loading.jl")
include("setup_hamiltonian.jl")
include("dKMC_hopping_rates.jl")


"""
    collect_charge_approximating_radii(dimension::Integer,N::Integer,disorder::Number,electronic_coupling::Number,
    bath_reorganisation_energy::Number,bath_cutoff_energy::Number,T::Number,kappa::AbstractFloat,
    accuracy::AbstractFloat,landscape_iterations::Integer,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)

Collects the approximating hopping and hamiltonian radii. If the radii have already been calculated, then it collects 
them from their saved location. Otherwise, it will calculate them by running charge_approximating_radii() on many 
realisations of disorder, splitting the work among available processes. The mean approximating radii are then saved 
for future use.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorder`: Energetic disorder (in meV).
- `electronic_coupling`: Nearest neighbour electronic coupling (in meV).
- `bath_reorganisation_energy`: Reorganisation energy of the bath (in meV).
- `bath_cutoff_energy`: Cutoff energy of the bath (in meV).
- `T`: Temperature (in K).
- `kappa`: Renormalisation constant for electronic couplings by the polaron transformation.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `landscape_iterations`: Number of realisations of disorder that approximating radii are calculated and averaged across.
- `K_tot`: Precaclulated K values required for dKMC rate calculations.
- `E_step`: Energy step used for values that K_tot is calculated at.
- `E_limit`: Energy limit used for values that K_tot is calculated at.

# Output:
- `hamiltonian_radius`: Mean hamiltonian radius cutoff required for outgoing rate sum to converge to a chosen accuracy.
- `hopping_radius`: Mean hopping radius cutoff required for outgoing rate sum to converge to a chosen accuracy.

"""
function collect_charge_approximating_radii(dimension::Integer,N::Integer,disorder::Number,electronic_coupling::Number,bath_reorganisation_energy::Number,bath_cutoff_energy::Number,T::Number,kappa::AbstractFloat,accuracy::AbstractFloat,landscape_iterations::Integer,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)
    
    #Create path for data to be stored.
    DATA_PATH = string(dirname(dirname(dirname(@__DIR__))),"/raw_data")
    if !isdir(DATA_PATH)
        mkdir(DATA_PATH)
    end    

    #Checking if the approximating hamiltonian and hopping radii have already been calculated and saved.
    if isfile("$DATA_PATH/hopping_radii/$dimension"*"D_σ"*string(disorder)*"_V$electronic_coupling"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$bath_reorganisation_energy" * "_ωc$bath_cutoff_energy" *"_T$T"*".dat")
                
        #If they are already calculated, read them from their saved location.
        hamiltonian_radius = readdlm("$DATA_PATH/hamiltonian_radii/$dimension"*"D_σ"*string(disorder)*"_V$electronic_coupling"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$bath_reorganisation_energy" * "_ωc$bath_cutoff_energy" *"_T$T"*".dat")[1]
        hopping_radius = readdlm("$DATA_PATH/hopping_radii/$dimension"*"D_σ"*string(disorder)*"_V$electronic_coupling"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$bath_reorganisation_energy" * "_ωc$bath_cutoff_energy" *"_T$T"*".dat")[1]

    else
        
        #If they are not already calculated, calculate them.
        results = pmap((disorders)->charge_approximating_radii(dimension,N,disorders,electronic_coupling,bath_reorganisation_energy,kappa,accuracy,K_tot,E_step,E_limit),Int.(ones(landscape_iterations)).*disorder)
        hamiltonian_radius,hopping_radius = mean(results[findall(x->x!=[0,0],results)])
        
        #Save the approximating hamiltonian and hopping radii to file for future use.
        if !isdir("$DATA_PATH/hamiltonian_radii")
            mkdir("$DATA_PATH/hamiltonian_radii")
        end
        writedlm("$DATA_PATH/hamiltonian_radii/$dimension"*"D_σ"*string(disorder)*"_V$electronic_coupling"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$bath_reorganisation_energy" * "_ωc$bath_cutoff_energy" *"_T$T"*".dat",hamiltonian_radius)
        
        if !isdir("$DATA_PATH/hopping_radii")
            mkdir("$DATA_PATH/hopping_radii")
        end
        writedlm("$DATA_PATH/hopping_radii/$dimension"*"D_σ"*string(disorder)*"_V$electronic_coupling"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$bath_reorganisation_energy" * "_ωc$bath_cutoff_energy" *"_T$T"*".dat",hopping_radius)
        
    end

    return hamiltonian_radius, hopping_radius

end


"""
    charge_approximating_radii(dimension::Integer,N::Integer,disorder::Number,electronic_coupling::Number,
    bath_reorganisation_energy::Number,kappa::AbstractFloat,accuracy::AbstractFloat,K_tot::Matrix{ComplexF64},
    E_step::Number,E_limit::Number)

Calculates the approximating hopping and hamiltonian radii on a single realisation of disorder. It does so by 
periodically increasing the radii until the total outgoing rate sum converges to a chosen accuracy.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorder`: Energetic disorder (in meV).
- `electronic_coupling`: Nearest neighbour electronic coupling (in meV).
- `bath_reorganisation_energy`: Reorganisation energy of the bath (in meV).
- `kappa`: Renormalisation constant for electronic couplings by the polaron transformation.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `K_tot`: Precaclulated K values required for dKMC rate calculations.
- `E_step`: Energy step used for values that K_tot is calculated at.
- `E_limit`: Energy limit used for values that K_tot is calculated at.

# Output:
- `hamiltonian_radius`: Hamiltonian radius cutoff required for outgoing rate sum to converge to a chosen accuracy.
- `hopping_radius`: Hopping radius cutoff required for outgoing rate sum to converge to a chosen accuracy.

"""
function charge_approximating_radii(dimension::Integer,N::Integer,disorder::Number,electronic_coupling::Number,bath_reorganisation_energy::Number,kappa::AbstractFloat,accuracy::AbstractFloat,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)

    #Assign energies to every site from a Gaussian distribution with standard deviation given by the disorder.
    site_energies = disorder.*randn(N^dimension)

    #Creating lists to track rate sums and hopping_radius at each increase of hopping_radius.
    total_hopping_rate_sum = [0.0]
    hopping_radii = [0]

    #Start hopping and hamiltonian radius at 1.
    hopping_radius = 1
    hamiltonian_radius = 1

    #Defining the initial position of the charge, in the middle of the lattice.	
    current_location = ones(dimension) .* (N/2)

    #Calculating the original system Hamiltonian, the polaron-trandformed system Hamiltonian, postion vectors, and coupling vector.
    H,Ht,r,transformed_coupling,site_indexes = setup_hamiltonian.current_charge_transport_hamiltonian(dimension,N,site_energies,electronic_coupling,bath_reorganisation_energy,kappa,hamiltonian_radius,current_location)

    #Calculating eigenvectors and eigenvalues in the energy eigenbasis.
    evals,evecs = eigen(Ht)

    #Calculating the expectation value of positions of energy eigenstates.
    centres = setup_hamiltonian.compute_centres(dimension,evecs,r)

    #Choose initial state closest to middle of system.
    current_state = argmin(sum([(centres[:,i] .- current_location[i]).^2 for i in 1:dimension]))
    current_location = centres[current_state,:]

    #Calculating distance from every state to the current state.
    distances_squared = sum([(centres[:,i] .-  current_location[i]).^2 for i in 1:dimension])

    #Create an empty vector to store hamiltonian information at each increase in hamiltonian radius, and then store the current hamiltonian information.
    stored_hamiltonian_information = Vector{Any}(zeros(Int(round(N/2))))
    stored_hamiltonian_information[hamiltonian_radius] = [H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared]
    
    #Increase hopping_radius until the rate sum converges to a chosen accuracy.
    for a = 1:Int(round(N/2))

        #Increase hopping radius.
        hopping_radius += 1
        
        #Save the total hopping rate sum from previous increment.
        hopping_rate_sum = [total_hopping_rate_sum[end]]

        #Creating vector to track hamiltonian radii at each increase of hopping_radius.
        hamiltonian_radii = [hamiltonian_radius]

        #Increase hamiltonian_radius until the rate sum to states within current hopping_radius converges to a chosen accuracy.
        for b = 1:Int(round(N/2))

            #Increase exciton hamiltonian radius.
            if b != 1
                hamiltonian_radius += 1
            elseif b == 1 && hopping_radius > hamiltonian_radius
                hamiltonian_radius += 1
            end

            #Checking whether the Hamiltonian got bigger than the predetermined system size N.
            if hamiltonian_radius >= N/2 - 1 
                if hopping_rate_sum[end] != 0
                    @warn "Convergence not reached while calculating approximating radii. Increase the length of the system (N) to avoid this error."
                end 
                return [0, 0]
            end

            #If the hamiltonian information at current hamiltonian_radius hasn't already been calculated and stored, calculate and store it.
            if stored_hamiltonian_information[hamiltonian_radius] == 0
                previous_eigenstate = evecs[:,current_state]
                previous_site_indexes = copy(site_indexes)
                H,Ht,r,transformed_coupling,site_indexes = setup_hamiltonian.current_charge_transport_hamiltonian(dimension,N,site_energies,electronic_coupling,bath_reorganisation_energy,kappa,hamiltonian_radius,current_location)
                evals,evecs = eigen(Ht)
                centres = setup_hamiltonian.compute_centres(dimension,evecs,r)
                current_state = argmax(abs2.((previous_eigenstate[findall(in(site_indexes),previous_site_indexes)]' * evecs[filter(!iszero,something.(indexin(previous_site_indexes,site_indexes),0)),:])[:]))
                current_location = centres[current_state,:]
                distances_squared = sum([(centres[:,i] .-  current_location[i]).^2 for i in 1:dimension])
                stored_hamiltonian_information[hamiltonian_radius] = [H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared]
            else
                H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared = stored_hamiltonian_information[hamiltonian_radius]
            end

            #Finding states that are within the new hopping_radius, but not the previous hopping_radius.
            accessible_states = findall(x-> 0 < x <= hopping_radius^2, distances_squared)
            
            #If there are no accessible_states within current hamiltonian_radius, increase the radii and set up and rediagonalise a larger subset of the Hamiltonian.
            if isempty(accessible_states)
                hopping_radius += 1
                if hopping_radius > hamiltonian_radius
                    hamiltonian_radius += 1
                    if stored_hamiltonian_information[hamiltonian_radius] == 0
                        previous_eigenstate = evecs[:,current_state]
                        previous_site_indexes = copy(site_indexes)
                        H,Ht,r,transformed_coupling,site_indexes = setup_hamiltonian.current_charge_transport_hamiltonian(dimension,N,site_energies,electronic_coupling,bath_reorganisation_energy,kappa,hamiltonian_radius,current_location)
                        evals,evecs = eigen(Ht)
                        centres = setup_hamiltonian.compute_centres(dimension,evecs,r)
                        current_state = argmax(abs2.((previous_eigenstate[findall(in(site_indexes),previous_site_indexes)]' * evecs[filter(!iszero,something.(indexin(previous_site_indexes,site_indexes),0)),:])[:]))
                        current_location = centres[current_state,:]
                        distances_squared = sum([(centres[:,i] .-  current_location[i]).^2 for i in 1:dimension])
                        stored_hamiltonian_information[hamiltonian_radius] = [H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared]
                    else
                        H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared = stored_hamiltonian_information[hamiltonian_radius]
                    end        
                end
                accessible_states = findall(x-> 0 < x <= hopping_radius^2, distances_squared)
            end

            #Calculating hopping rates to all states in accessible_states.
            hopping_rates = zeros(length(accessible_states))
            current_state_relevant_sites = dKMC_hopping_rates.relevant_sites(evecs[:,current_state],accuracy)
            for (f,destination_state) in enumerate(accessible_states)
                destination_state_relevant_sites = dKMC_hopping_rates.relevant_sites(evecs[:,destination_state],accuracy)
                hopping_rate = dKMC_hopping_rates.charge_transport_dKMC_rate(current_state,destination_state,transformed_coupling,evals,evecs,K_tot,E_step,E_limit,current_state_relevant_sites,destination_state_relevant_sites)
                if hopping_rate > 0
                    hopping_rates[f] = hopping_rate
                end
            end
            push!(hopping_rate_sum, sum(hopping_rates))

            #Checking whether the hopping_rate_sum at hopping_radius has converged at hamiltonian_radius to a given accuracy.
            if hopping_rate_sum[b]/hopping_rate_sum[b+1] > accuracy && hopping_rate_sum[b+1]/hopping_rate_sum[b] > accuracy
                push!(total_hopping_rate_sum,hopping_rate_sum[b])
                hamiltonian_radius = hamiltonian_radii[b]
                push!(hopping_radii,hopping_radius)
                break
            end

            #Store the hamiltonian_radius at this iteration.
            push!(hamiltonian_radii,hamiltonian_radius)

        end

        #Checking whether the total Redfield sum has convered at hopping_radius and hamiltonian_radius to a given accuracy.
        if total_hopping_rate_sum[a]/total_hopping_rate_sum[a+1] > accuracy && total_hopping_rate_sum[a+1]/total_hopping_rate_sum[a] > accuracy
            hopping_radius = hopping_radii[a]
            break
        end

    end

    return [hamiltonian_radius, hopping_radius]

end


"""
    collect_exciton_approximating_radii(dimension::Integer,N::Integer,exciton_disorder::Number,
    transition_dipole_moment::Number,epsilon_r::Number,exciton_bath_reorganisation_energy::Number,
    exciton_bath_cutoff_energy::Number,T::Number,kappa::AbstractFloat,site_spacing::Number,accuracy::AbstractFloat,
    landscape_iterations::Integer,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)

Collects the approximating exciton hopping and hamiltonian radii. If the radii have already been calculated, then it 
collects them from their saved location. Otherwise, it will calculate them by running exciton_approximating_radii() on 
many realisations of disorder, splitting the work among available processes. The mean approximating radii are then 
saved for future use.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `exciton_disorder`: Excitonic disorder (in meV).
- `transition_dipole_moment`: Magnitude of the transition dipole moment on every site (in D).
- `exciton_bath_reorganisation_energy`: Reorganisation energy of the exciton's bath (in meV).
- `exciton_exciton_bath_cutoff_energy`: Cutoff energy of the exciton's bath (in meV).
- `T`: Temperature (in K).
- `kappa`: Renormalisation constant for excitonic couplings by the polaron transformation.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `landscape_iterations`: Number of realisations of disorder that approximating radii are calculated and averaged across.
- `K_tot`: Precaclulated K values required for dKMC rate calculations.
- `E_step`: Energy step used for values that K_tot is calculated at.
- `E_limit`: Energy limit used for values that K_tot is calculated at.

# Output:
- `exciton_hamiltonian_radius`: Mean exciton hamiltonian radius cutoff required for outgoing rate sum to converge to a chosen accuracy.
- `exciton_hopping_radius`: Mean exciton_hopping radius cutoff required for outgoing rate sum to converge to a chosen accuracy.

"""
function collect_exciton_approximating_radii(dimension::Integer,N::Integer,exciton_disorder::Number,transition_dipole_moment::Number,epsilon_r::Number,exciton_bath_reorganisation_energy::Number,exciton_bath_cutoff_energy::Number,T::Number,kappa::AbstractFloat,site_spacing::Number,accuracy::AbstractFloat,landscape_iterations::Integer,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)
    
    #Create path for data to be stored.
    DATA_PATH = string(dirname(dirname(dirname(@__DIR__))),"/raw_data")
    if !isdir(DATA_PATH)
        mkdir(DATA_PATH)
    end    
    
    #Checking if the approximating exciton hamiltonian and hopping radii have already been calculated and saved.
    if isfile("$DATA_PATH/hopping_radii/$dimension"*"D_σ$exciton_disorder"*"_μ$transition_dipole_moment"*"_εr$epsilon_r"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$exciton_bath_reorganisation_energy" * "_ωc$exciton_bath_cutoff_energy" *"_T$T"*".dat")
                
        #If they are already calculated, read them from their saved location.
        exciton_hamiltonian_radius = readdlm("$DATA_PATH/hamiltonian_radii/$dimension"*"D_σ$exciton_disorder"*"_μ$transition_dipole_moment"*"_εr$epsilon_r"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$exciton_bath_reorganisation_energy" * "_ωc$exciton_bath_cutoff_energy" *"_T$T"*".dat")[1]
        exciton_hopping_radius = readdlm("$DATA_PATH/hopping_radii/$dimension"*"D_σ$exciton_disorder"*"_μ$transition_dipole_moment"*"_εr$epsilon_r"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$exciton_bath_reorganisation_energy" * "_ωc$exciton_bath_cutoff_energy" *"_T$T"*".dat")[1]

    else
        
        #If they are not already calculated, calculate them.
        results = pmap((exciton_disorders)->exciton_approximating_radii(dimension,N,exciton_disorders,transition_dipole_moment,epsilon_r,exciton_bath_reorganisation_energy,kappa,site_spacing,accuracy,K_tot,E_step,E_limit), exciton_disorder.*ones(landscape_iterations))
        exciton_hamiltonian_radius,exciton_hopping_radius = mean(results[findall(x->x!=[0,0],results)])
   
        #Save the approximating exciton hamiltonian and hopping radii to file for future use.        
        if !isdir("$DATA_PATH/hamiltonian_radii")
            mkdir("$DATA_PATH/hamiltonian_radii")
        end
        writedlm("$DATA_PATH/hamiltonian_radii/$dimension"*"D_σ$exciton_disorder"*"_μ$transition_dipole_moment"*"_εr$epsilon_r"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$exciton_bath_reorganisation_energy" * "_ωc$exciton_bath_cutoff_energy" *"_T$T"*".dat",exciton_hamiltonian_radius)
        
        if !isdir("$DATA_PATH/hopping_radii")
            mkdir("$DATA_PATH/hopping_radii")
        end
        writedlm("$DATA_PATH/hopping_radii/$dimension"*"D_σ$exciton_disorder"*"_μ$transition_dipole_moment"*"_εr$epsilon_r"*"_accuracy$accuracy"*"_iters$landscape_iterations" * "_λ$exciton_bath_reorganisation_energy" * "_ωc$exciton_bath_cutoff_energy" *"_T$T"*".dat",exciton_hopping_radius)
        
    end

    return exciton_hamiltonian_radius, exciton_hopping_radius

end


"""
    exciton_approximating_radii(dimension::Integer,N::Integer,exciton_disorder::Number,
    transition_dipole_moment::Number,epsilon_r::Number,exciton_bath_reorganisation_energy::Number,kappa::AbstractFloat,
    site_spacing::Number,accuracy::AbstractFloat,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)	
    
Calculates the approximating exciton hopping and hamiltonian radii. It does so by periodically increasing them until 
the total outgoing rate sum converges to a chosen accuracy.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `exciton_disorder`: Exciton disorder (in meV).
- `transition_dipole_moment`: Magnitude of the transition dipole moment on every site (in D).
- `epsilon_r`: Dielectric constant.
- `exciton_bath_reorganisation_energy`: Reorganisation energy of the exciton bath (in meV). 
- `kappa`: Renormalisation constant for excitonic couplings by the polaron transformation.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `K_tot`: Precaclulated exciton K values required for dKMC rate calculations.
- `E_step`: Energy step used in K value calculation.
- `E_limit`: Energy limit used in K value calculation.

# Output:
- `exciton_hamiltonian_radius`: Exciton hamiltonian radius cutoff required for outgoing rate sum to converge to a chosen accuracy.
- `exciton_hopping_radius`: Exciton hopping radius cutoff required for outgoing rate sum to converge to a chosen accuracy.

"""
function exciton_approximating_radii(dimension::Integer,N::Integer,exciton_disorder::Number,transition_dipole_moment::Number,epsilon_r::Number,exciton_bath_reorganisation_energy::Number,kappa::AbstractFloat,site_spacing::Number,accuracy::AbstractFloat,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)

    #Assign exciton energies to every site from a Gaussian distribution with standard deviation given by the exciton_disorder.
    exciton_site_energies = exciton_disorder.*randn(N^dimension)

    #Assign a randomly orientated dipole orientation to every site.
    dipole_orientations = setup_hamiltonian.assign_exciton_dipole_orientations(dimension,N)

    #Creating lists to track rate sums and exciton_hopping_radius at each increase of exciton_hopping_radius.
    total_exciton_hopping_rate_sum = [0.0]
    exciton_hopping_radii = [0]

    #Start exciton hopping and hamiltonian radius at 1.
    exciton_hopping_radius = 1
    exciton_hamiltonian_radius = 1

    #Defining the initial position of the exciton, in the middle of the lattice.	
    current_location = ones(dimension).* (N/2)

    #Calculating the original system Hamiltonian, the polaron-trandformed system Hamiltonian, postion vectors, and coupling vector.
    H,Ht,r,transformed_coupling,site_indexes = setup_hamiltonian.current_exciton_transport_hamiltonian(dimension,N,exciton_site_energies,dipole_orientations,transition_dipole_moment,epsilon_r,exciton_bath_reorganisation_energy,kappa,site_spacing,exciton_hamiltonian_radius,current_location)

    #Calculating eigenvectors and eigenvalues in the energy eigenbasis.
    evals,evecs = eigen(Ht)

    #Calculating the expectation value of positions of energy eigenstates.
    centres = setup_hamiltonian.compute_centres(dimension,evecs,r)

    #Choose initial state closest to middle of system.
    current_state = argmin(sum([(centres[:,i] .- current_location[i]).^2 for i in 1:dimension]))
    current_location = centres[current_state,:]

    #Calculating distance from every state to the current state.
    distances_squared = sum([(centres[:,i] .-  current_location[i]).^2 for i in 1:dimension])

    #Create an empty vector to store hamiltonian information at each increase in hamiltonian radius, and then store the current hamiltonian information.
    stored_hamiltonian_information = Vector{Any}(zeros(Int(N/2)))
    stored_hamiltonian_information[exciton_hamiltonian_radius] = [H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared]

    #Increase exciton_hopping_radius until the rate sum converges to a chosen accuracy.
    for a = 1:Int(round(N/2))

        #Increase exciton hopping radius.
        exciton_hopping_radius += 1
        
        #Save the total hopping rate sum from previous increment.
        exciton_hopping_rate_sum = [total_exciton_hopping_rate_sum[end]]

        #Creating vector to track hamiltonian radii at each increase of hopping_radius.
        exciton_hamiltonian_radii = [exciton_hamiltonian_radius]

        #Increase exciton_hamiltonian_radius until the rate sum to states within current exciton_hopping_radius converges to a chosen accuracy.
        for b = 1:Int(round(N/2))

            #Increase exciton hamiltonian radius.
            if b != 1
                exciton_hamiltonian_radius += 1
            elseif b == 1 && exciton_hopping_radius > exciton_hamiltonian_radius
                exciton_hamiltonian_radius += 1
            end

            #Checking whether the Hamiltonian got bigger than the predetermined system size N.
            if exciton_hamiltonian_radius >=  N/2 - 1
                if exciton_hopping_rate_sum[end] != 0
                    @warn "Convergence not reached while calculating approximating radii. Increase the length of the system (N) to avoid this error."
                end 
                return [0, 0]
            end

            #If the hamiltonian information at current hamiltonian_radius hasn't already been calculated and stored, calculate and store it.
            if stored_hamiltonian_information[exciton_hamiltonian_radius] == 0
                previous_eigenstate = evecs[:,current_state]
                previous_site_indexes = copy(site_indexes)
                H,Ht,r,transformed_coupling,site_indexes = setup_hamiltonian.current_exciton_transport_hamiltonian(dimension,N,exciton_site_energies,dipole_orientations,transition_dipole_moment,epsilon_r,exciton_bath_reorganisation_energy,kappa,site_spacing,exciton_hamiltonian_radius,current_location)
                evals,evecs = eigen(Ht)
                centres = setup_hamiltonian.compute_centres(dimension,evecs,r)
                current_state = argmax(abs2.((previous_eigenstate[findall(in(site_indexes),previous_site_indexes)]' * evecs[filter(!iszero,something.(indexin(previous_site_indexes,site_indexes),0)),:])[:]))
                current_location = centres[current_state,:]
                distances_squared = sum([(centres[:,i] .-  current_location[i]).^2 for i in 1:dimension])
                stored_hamiltonian_information[exciton_hamiltonian_radius] = [H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared]   
            else
                H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared = stored_hamiltonian_information[exciton_hamiltonian_radius]
            end

            #Finding states that are within the new exciton_hopping_radius, but not the previous exciton_hopping_radius.
            accessible_states = findall(x-> 0 < x < exciton_hopping_radius^2, distances_squared)

            #If there are no accessible_states within current exciton_hamiltonian_radius, increase the radii and set up and rediagonalise a larger subset of the Hamiltonian.
            if isempty(accessible_states)
                exciton_hopping_radius += 1
                if exciton_hopping_radius > exciton_hamiltonian_radius
                    exciton_hamiltonian_radius += 1
                    if stored_hamiltonian_information[exciton_hamiltonian_radius] == 0
                        previous_eigenstate = evecs[:,current_state]
                        previous_site_indexes = copy(site_indexes)
                        H,Ht,r,transformed_coupling,site_indexes = setup_hamiltonian.current_exciton_transport_hamiltonian(dimension,N,exciton_site_energies,dipole_orientations,transition_dipole_moment,epsilon_r,exciton_bath_reorganisation_energy,kappa,site_spacing,exciton_hamiltonian_radius,current_location)
                        evals,evecs = eigen(Ht)
                        centres = setup_hamiltonian.compute_centres(dimension,evecs,r)
                        current_state = argmax(abs2.((previous_eigenstate[findall(in(site_indexes),previous_site_indexes)]' * evecs[filter(!iszero,something.(indexin(previous_site_indexes,site_indexes),0)),:])[:]))
                        current_location = centres[current_state,:]
                        distances_squared = sum([(centres[:,i] .-  current_location[i]).^2 for i in 1:dimension])
                        stored_hamiltonian_information[exciton_hamiltonian_radius] = [H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared]   
                    else
                        H,Ht,r,transformed_coupling,site_indexes,evals,evecs,centres,current_state,current_location,distances_squared = stored_hamiltonian_information[exciton_hamiltonian_radius]
                    end
                end
                accessible_states = findall(x-> 0 < x < exciton_hopping_radius^2, distances_squared)
            end

            #Calculating hopping rates to all states in accessible_states.
            hopping_rates = zeros(length(accessible_states))
            current_state_relevant_sites = dKMC_hopping_rates.relevant_sites(evecs[:,current_state],accuracy)
            for (f,destination_state) in enumerate(accessible_states)
                destination_state_relevant_sites = dKMC_hopping_rates.relevant_sites(evecs[:,destination_state],accuracy)
                hopping_rate = dKMC_hopping_rates.exciton_transport_dKMC_rate(current_state,destination_state,transformed_coupling,evals,evecs,K_tot,E_step,E_limit,current_state_relevant_sites,destination_state_relevant_sites)
                if hopping_rate > 0
                    hopping_rates[f] = hopping_rate
                end
            end
            push!(exciton_hopping_rate_sum, sum(hopping_rates))

            #Checking whether the exciton_hopping_rate_sum at exciton_hopping_radius has converged at exciton_hamiltonian_radius to a given accuracy.
            if exciton_hopping_rate_sum[b]/exciton_hopping_rate_sum[b+1] > accuracy && exciton_hopping_rate_sum[b+1]/exciton_hopping_rate_sum[b] > accuracy
                push!(total_exciton_hopping_rate_sum,exciton_hopping_rate_sum[b])
                exciton_hamiltonian_radius = exciton_hamiltonian_radii[b]
                push!(exciton_hopping_radii,exciton_hopping_radius)
                break
            end

            #Store the hamiltonian_radius at this iteration.
            push!(exciton_hamiltonian_radii,exciton_hamiltonian_radius)

        end 

        #Checking whether the total Redfield sum has convered at exciton_hopping_radius and exciton_hamiltonian_radius to a given accuracy.
        if total_exciton_hopping_rate_sum[a]/total_exciton_hopping_rate_sum[a+1] > accuracy && total_exciton_hopping_rate_sum[a+1]/total_exciton_hopping_rate_sum[a] > accuracy
            exciton_hopping_radius = exciton_hopping_radii[a]
            break
        end

    end

    return [exciton_hamiltonian_radius, exciton_hopping_radius]

end

end