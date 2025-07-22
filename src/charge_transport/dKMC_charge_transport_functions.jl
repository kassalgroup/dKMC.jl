module dKMC_charge_transport_functions

#Loading the required packages.
include("../shared_functions/package_loading.jl")
include("../shared_functions/approximating_radii_functions.jl")
include("../shared_functions/setup_hamiltonian.jl")
include("../shared_functions/constants.jl")
include("../shared_functions/sPTRE.jl")
include("../shared_functions/dKMC_hopping_rates.jl")


"""
    dKMC_charge_transport_results(dimension::Integer,N::Integer,disorder::Number,electronic_coupling::Number,
    bath_reorganisation_energy::Number,bath_cutoff_energy::Number,T::Number,site_spacing::Number,
    landscape_iterations::Integer,trajectory_iterations::Integer,accuracy::AbstractFloat,end_time::Number,
    number_of_sampling_times::Integer;phi_limit = round(100*constants.hbar/bath_cutoff_energy),
    phi_step = phi_limit/10000,E_limit = round(6*disorder),E_step = E_limit/10000) 

Produces the results for dKMC simulations of charge transport. To do so, it collects (or calculates) kappa 
(renormalisation constant), K_tot (Fourier transform of bath correlation function), and Hamiltonian and hopping radii. 
It then runs iterate_dKMC_charge_transport() on many realisations of disorder, parallelising the work amongst all of 
the available processes. Finally, it collects the mean trajectory information from each realisation of disorder, and 
averages them to return the charge-carrier mobility as well as mean-squared displacements, energies and IPRs at chosen 
sampling times.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorder`: Energetic disorder (in meV).
- `electronic_coupling`: Nearest neighbour electronic coupling (in meV).
- `bath_reorganisation_energy`: Reorganisation energy of the bath (in meV).
- `bath_cutoff_energy`: Cutoff energy of the bath (in meV).
- `T`: Temperature (in K).
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `landscape_iterations`: Number of realisations of disorder (or energetic landscapes) for simulations to be run on.   
- `trajectory_iterations`: Number of trajectories to be simulated on each realisation of disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC). 
- `end_time`: The simulation end time (in s).
- `number_of_sampling_times`: The number of times between 0 and end_time to sample the squared displacement, energies and IPRs.

# Optional arguments:
- `phi_step`: The time step used in calculating phi.
- `phi_limit`: The time limit used in calculating phi.
- `E_step`: Energy step used for values that K_tot is calculated at.
- `E_limit`: Energy limit used for values that K_tot is calculated at.

# Output:
- `mobility`: Charge carrier mobility and its uncertainty (in cm²V⁻¹s⁻¹).
- `sampling_time_squared_displacements`: The mean squared displacements (and standard errors of the mean) across many realisations of disorder at chosen sampling_times (in nm²). 
- `sampling_time_energies`: The mean energies (and standard errors of the mean) across many realisations of disorder at chosen sampling_times (in meV).
- `sampling_time_IPRs`: The mean IPRs (and standard errors of the mean) across many realisations of disorder at chosen sampling_times.
- `total_hops`: The mean number of hops in each trajectory across many realisations of disorder.

"""
function dKMC_charge_transport_results(dimension::Integer,N::Integer,disorder::Number,electronic_coupling::Number,bath_reorganisation_energy::Number,bath_cutoff_energy::Number,T::Number,site_spacing::Number,landscape_iterations::Integer,trajectory_iterations::Integer,accuracy::AbstractFloat,end_time::Number,number_of_sampling_times::Integer;phi_limit = round(100*constants.hbar/bath_cutoff_energy),phi_step = phi_limit/10000,E_limit = round(6*disorder),E_step = E_limit/10000)
        
    #Create the spectral density function, here a super-ohmic spectral density is used.
    J(E) = (bath_reorganisation_energy/2) * (E/bath_cutoff_energy)^3 * exp(-E/bath_cutoff_energy)

    #Calculating the renormalisation constant (kappa) for the electronic couplings following polaron transformation.
    kappa = sPTRE.calculate_kappa(J,T)

    #Reading the precomputed K values, or calculating them if not. K values are stored in a vector and called upon by looking at E and lambda+3 as the indices.
    K_tot = sPTRE.collect_K_results(J,E_step,E_limit,phi_step,phi_limit,bath_reorganisation_energy,bath_cutoff_energy,T)

    #Calculating the approximating hamiltonian and hopping radii for each particle and each material, or reading them from a saved file if previously calculated.
    hamiltonian_radius,hopping_radius = approximating_radii_functions.collect_charge_approximating_radii(dimension,N,disorder,electronic_coupling,bath_reorganisation_energy,bath_cutoff_energy,T,kappa,accuracy,landscape_iterations,K_tot,E_step,E_limit)

    #Creating a list of times which we sample the MSD, E and IPR over.
    sampling_times = collect(0:end_time/number_of_sampling_times:end_time)

    #We parallelise this function over many processes, to repeat it for many (landscape_iterations) realisations of energetic disorder. It calculates the mean squared displacements, energies and IPRs at sampling_times on a single landscape averaged over many (trajectory_iterations) times.
    results = pmap((disorders)->iterate_dKMC_charge_transport(dimension,N,disorders,electronic_coupling,bath_reorganisation_energy,kappa,trajectory_iterations,accuracy,end_time,sampling_times,hopping_radius,hamiltonian_radius,K_tot,E_step,E_limit), disorder.*ones(landscape_iterations))

    #Average the results over all landscape iterations.
    mean_squared_displacements = zeros(landscape_iterations,length(sampling_times))
    mean_energies = zeros(landscape_iterations,length(sampling_times))
    mean_IPRs = zeros(landscape_iterations,length(sampling_times))
    mean_total_hops = zeros(landscape_iterations)
    for i = 1:landscape_iterations
        mean_squared_displacements[i,:] = results[i][1]
        mean_energies[i,:] = results[i][2]
        mean_IPRs[i,:] = results[i][3]
        mean_total_hops[i] = results[i][4]
    end
    sampling_time_squared_displacements = constants.m2_to_nm2(site_spacing^2) .* [mean(mean_squared_displacements,dims=1); std(mean_squared_displacements,dims=1)./sqrt(landscape_iterations)]
    sampling_time_energies = [mean(mean_energies,dims=1); std(mean_energies,dims=1)./sqrt(landscape_iterations)]
    sampling_time_IPRs = [mean(mean_IPRs,dims=1); std(mean_IPRs,dims=1)./sqrt(landscape_iterations)]
    total_hops = [mean(mean_total_hops); std(mean_total_hops)./sqrt(landscape_iterations)]

    #Calculate the mobility.
    mobility = calculate_mobility(sampling_time_squared_displacements[1,:],sampling_time_squared_displacements[2,:],sampling_times,dimension,T,site_spacing)

    return mobility, sampling_time_squared_displacements, sampling_time_energies, sampling_time_IPRs, total_hops

end


"""
    iterate_dKMC_charge_transport(dimension::Integer,N::Integer,disorder::Number,electronic_coupling::Number,
    bath_reorganisation_energy::Number,kappa::AbstractFloat,trajectory_iterations::Integer,end_time::Number,
    sampling_times::Vector{<:AbstractFloat},accuracy::AbstractFloat,hopping_radius::AbstractFloat,
    hamiltonian_radius::AbstractFloat,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)

Repeats dKMC_charge_transport() for many trajectories on a single realisation of energetic disorder, before averaging 
mean squared displacements, energies and IPRs at chosen sampling times.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorder`: Energetic disorder (in meV).
- `electronic_coupling`: Nearest neighbour electronic coupling (in meV).
- `bath_reorganisation_energy`: Reorganisation energy of the bath (in meV).  
- `kappa`: Renormalisation constant for electronic couplings by the polaron transformation.
- `trajectory_iterations`: Number of trajectories simulated on each realsation of energetic disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `end_time`: The simulation end time (in s).
- `sampling_times`: The times at which KMC trajectories are sampled for averaging.
- `hopping_radius`: Precalculated distance that charges can hop. 
- `hamiltonian_radius`: Precalculated distance for sites that are included in Hamiltonian subset.
- `K_tot`: Precaclulated K values required for dKMC rate calculations.
- `E_step`: Energy step used for values that K_tot is calculated at.
- `E_limit`: Energy limit used for values that K_tot is calculated at.

# Output:
- `mean_squared_displacements`: The mean squared displacements (in site spacing) across many trajectory iterations at chosen sampling_times.
- `mean_energies`: The mean energies (in meV) across many trajectory iterations at chosen sampling_times.
- `mean_IPRs`: The mean IPRs across many trajectory iterations at chosen sampling_times.
- `mean_total_hops`: The mean total number of hops across many trajectory iterations.

"""
function iterate_dKMC_charge_transport(dimension::Integer,N::Integer,disorder::Number,electronic_coupling::Number,bath_reorganisation_energy::Number,kappa::AbstractFloat,trajectory_iterations::Integer,accuracy::AbstractFloat,end_time::Number,sampling_times::Vector{<:AbstractFloat},hopping_radius::AbstractFloat,hamiltonian_radius::AbstractFloat,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)

    #Assign energies to every site from a Gaussian distribution with standard deviation given by the disorder.
    site_energies = disorder.*randn(N^dimension)

    #Create lists to track the squared displacements, energies, and IPRs at chosen sampling times, as well as total number of hops.
    sample_time_squared_displacements = zeros(trajectory_iterations,length(sampling_times))
    sample_time_energies = zeros(trajectory_iterations,length(sampling_times))
    sample_time_IPRs = zeros(trajectory_iterations,length(sampling_times))
    total_hops = Vector{Integer}(zeros(trajectory_iterations))

    #Looping over the same landscape a total of trajectory_iterations amount of times.
    for i = 1:trajectory_iterations

        #Performing the dKMC procedure for a single trajectory on a single realisation of energetic disorder.
        squared_displacements,energies,IPRs,times = dKMC_charge_transport(dimension,N,site_energies,electronic_coupling,bath_reorganisation_energy,kappa,accuracy,end_time,hopping_radius,hamiltonian_radius,K_tot,E_step,E_limit)

        #Converting results from irregular time intervals times to those at regular time intervals sampling_times.
        close_times = Vector{Integer}(zeros(length(sampling_times)))
        for j = eachindex(sampling_times)
            close_times[j] = Int(findlast(x-> x.<=sampling_times[j],times))
        end
        sample_time_squared_displacements[i,:] = squared_displacements[close_times]
        sample_time_energies[i,:] = energies[close_times]
        sample_time_IPRs[i,:] = IPRs[close_times]
        total_hops[i] = length(times) - 1

    end

    #Averaging over all of the trajectories.
    mean_squared_displacements = mean(sample_time_squared_displacements,dims=1)
    mean_energies = mean(sample_time_energies,dims=1)
    mean_IPRs = mean(sample_time_IPRs,dims=1)
    mean_total_hops = mean(total_hops)

    return mean_squared_displacements, mean_energies, mean_IPRs, mean_total_hops

end


"""
    dKMC_charge_transport(dimension::Integer,N::Integer,site_energies::Vector{<:AbstractFloat},
    electronic_coupling::Number,bath_reorganisation_energy::Number,kappa::AbstractFloat,end_time::Number,
    accuracy::AbstractFloat,hopping_radius::AbstractFloat,hamiltonian_radius::AbstractFloat,
    K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)

Propagates the dKMC procedure for a single charge transport trajectory. It starts with the charge in the state closest 
to the middle of the landscape, and terminates once the time is updated past the chosen end time. It returns discrete 
values of squared displacement, energy, IPR, and time at each hop.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `site_energies`: List of energies of every site in (meV), of length N^dimension.
- `electronic_coupling`: Nearest neighbour electronic coupling (in meV).
- `bath_reorganisation_energy`: Reorganisation energy of the bath (in meV).  
- `kappa`: Renormalisation constant for electronic couplings by the polaron transformation.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `end_time`: The simulation end time (in s).
- `hopping_radius`: Precalculated distance that charges can hop. 
- `hamiltonian_radius`: Precalculated distance for sites that are included in Hamiltonian subset.
- `K_tot`: Precaclulated K values required for dKMC rate calculations.
- `E_step`: Energy step used for values that K_tot is calculated at.
- `E_limit`: Energy limit used for values that K_tot is calculated at.

# Output:
- `squared_displacements`: Squared displacements (in site spacing) after each hop.
- `energies`: Energies (in meV) after each hop.
- `IPRs`: IPRs after each hop.
- `times`: Times (in s) after each hop.

"""
function dKMC_charge_transport(dimension::Integer,N::Integer,site_energies::Vector{<:AbstractFloat},electronic_coupling::Number,bath_reorganisation_energy::Number,kappa::AbstractFloat,accuracy::AbstractFloat,end_time::Number,hopping_radius::AbstractFloat,hamiltonian_radius::AbstractFloat,K_tot::Matrix{ComplexF64},E_step::Number,E_limit::Number)

    #Defining the initial position of the charge, in the middle of the lattice.	
    current_location = ones(dimension).* (N/2)

    #Calculate the original system Hamiltonian, the polaron-trandformed system Hamiltonian, postion vectors, and coupling vector.
    H,Ht,r,transformed_coupling,site_indexes = setup_hamiltonian.current_charge_transport_hamiltonian(dimension,N,site_energies,electronic_coupling,bath_reorganisation_energy,kappa,hamiltonian_radius,current_location)

    #Calculating eigenvectors and eigenvalues in the energy eigenbasis.
    evals,evecs = eigen(Ht)

    #Calculating the expectation value of positions of energy eigenstates.
    centres = setup_hamiltonian.compute_centres(dimension,evecs,r)

    #Choose initial state closest to middle of system.
    current_state = argmin(sum([(centres[:,i] .-  current_location[i]).^2 for i in 1:dimension]))

    #Setting the initial conditions.
    t = 0.0

    #Creating lists to track time, positions, squared displacement, energy, and IPR.
    times = [t]
    current_location = centres[current_state,:]
    initial_location = copy(current_location)
    squared_displacements = [sum((initial_location .- current_location).^2)]
    energies = [evals[current_state]]
    IPRs = [1/sum(evecs[:,current_state].^4)]

    #Iterate the dKMC procedure until the time is updated past a predetermined end_time.
    while t < end_time

        #Check whether we are too close to the edges of the system.
        if iszero(current_location .- hamiltonian_radius .< 1) == false || iszero(current_location .+ hamiltonian_radius .> N) == false 
            @warn "Simulation terminated as the charge carrier got too close to the edge of the system. Increase the length of the system (N) to avoid this error."
            break
        end

        #Rediagonalise a new subset of the Hamiltonian.
        previous_eigenstate = evecs[:,current_state]
        previous_site_indexes = copy(site_indexes)
        H,Ht,r,transformed_coupling,site_indexes = setup_hamiltonian.current_charge_transport_hamiltonian(dimension,N,site_energies,electronic_coupling,bath_reorganisation_energy,kappa,hamiltonian_radius,current_location)
        evals,evecs = eigen(Ht)
        centres = setup_hamiltonian.compute_centres(dimension,evecs,r)
        current_state = argmax(((previous_eigenstate[findall(in(site_indexes),previous_site_indexes)]' * evecs[findall(in(previous_site_indexes),site_indexes),:])[:]).^2)
        current_location = centres[current_state,:]

        #Finding which states are accessible from the current state.
        accessible_states = findall(x-> 0 .< x .< hopping_radius^2, sum([(centres[:,i] .-  current_location[i]).^2 for i in 1:dimension]))

        #Calculating hopping rates to all states in accessible_states.
        hopping_rates = zeros(length(accessible_states))
        current_state_relevant_sites = sortperm(abs.(evecs[:,current_state]),rev=true)[1:findfirst(x->x>accuracy,cumsum(sort(abs.(evecs[:,current_state])./sum(abs.(evecs[:,current_state])),rev=true)))]
        for f = eachindex(accessible_states)
            destination_state_relevant_sites = sortperm(abs.(evecs[:,accessible_states[f]]),rev=true)[1:findfirst(x->x>accuracy,cumsum(sort(abs.(evecs[:,accessible_states[f]])./sum(abs.(evecs[:,accessible_states[f]])),rev=true)))]
            hopping_rate = dKMC_hopping_rates.charge_transport_dKMC_rate(current_state,accessible_states[f],transformed_coupling,evals,evecs,K_tot,E_step,E_limit,current_state_relevant_sites,destination_state_relevant_sites)
            if hopping_rate > 0
                hopping_rates[f] = hopping_rate
            end
        end
        
        #Calculating sum of all rates.
        hopping_rate_sum = sum(hopping_rates)

        #Prevents an error if the rate sum comes out to zero.
        if hopping_rate_sum == 0
            break
        end

        #Selected hop is chosen propabalistically in proportion to the hopping rate. 
        ran = rand()
        current_state = accessible_states[findfirst(x->x>ran*hopping_rate_sum,cumsum(hopping_rates))]
        current_location = centres[current_state,:]

        #Time updated by calculating elapsed time.
        ran2 = rand()
        t += constants.from_femto(hopping_rate_sum^(-1)*log(1/ran2))
        
        #Recording the new time, positions, squared displacement, energy, and IPR.
        push!(times,t)
        push!(squared_displacements,sum((initial_location .- current_location).^2))
        push!(energies,evals[current_state])
        push!(IPRs,1/sum(evecs[:,current_state].^4))

    end

    return squared_displacements, energies, IPRs, times

end


"""
    calculate_mobility(mean_squared_displacement::Vector{<:AbstractFloat},
    standard_error_of_mean::Vector{<:AbstractFloat},sampling_times::Vector{<:AbstractFloat},dimension::Integer,
    T::Number,site_spacing::Number)

Calculates the charge carrier mobillity from the provided mean squared displacement data, by traking the 
gradient over the second half of the sampling times.

# Arguments:
- `mean_squared_displacement`: Mean squared displacements at sampling times (in site spacings).
- `standard_error_of_mean`: Standard errors of the mean squared displacement provided above. 
- `sampling_times`: The times at which KMC trajectories are sampled for averaging.
- `dimension`: Dimension of the system (1, 2, or 3).
- `T`: Temperature (in K).
- `site_spacing`: Distance between sites in the cubic lattice (in m).

# Output:
- `mobility`: Charge carrier mobility and its uncertainty (in cm²V⁻¹s⁻¹).

"""
function calculate_mobility(mean_squared_displacement::Vector{<:AbstractFloat},standard_error_of_mean::Vector{<:AbstractFloat},sampling_times::Vector{<:AbstractFloat},dimension::Integer,T::Number,site_spacing::Number)

    #Use only mean squared displacements over the second half of the sampling times.
    starting_index = Int(round(length(mean_squared_displacement)/2))
    x_data = collect(sampling_times)[starting_index:length(mean_squared_displacement)]
    y_data = mean_squared_displacement[starting_index:length(mean_squared_displacement)] .± standard_error_of_mean[starting_index:length(mean_squared_displacement)]
    
    #Find the line of best fit.
    result = Polynomials.fit(x_data,y_data,1)
    
    #Calculate the mobility using the gradient.
    mobility = result[1].*((constants.m2_to_cm2(site_spacing^2)*constants.q_si)/(2*dimension*constants.kb_si*T))

    return [Measurements.value(mobility), Measurements.uncertainty(mobility)]

end

end
