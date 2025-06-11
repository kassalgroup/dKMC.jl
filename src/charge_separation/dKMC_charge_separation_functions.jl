module dKMC_charge_separation_functions

#Loading the required packages.
include("../shared_functions/package_loading.jl")
include("../shared_functions/setup_hamiltonian.jl")
include("../shared_functions/constants.jl")
include("../shared_functions/sPTRE.jl")
include("../shared_functions/dKMC_hopping_rates.jl")
include("../shared_functions/approximating_radii_functions.jl")


"""
    dKMC_charge_separation_results(dimension::Integer,N::Integer,disorders::Vector{<:Number},
    donor_HOMO_acceptor_LUMO_gap::Number,electronic_couplings::Vector{<:Number},epsilon_r::Number,
    site_spacing::Number,CT_lifetime::Number,bath_reorganisation_energies::Vector{<:Number},
    bath_cutoff_energies::Vector{<:Number},T::Number,landscape_iterations::Integer,trajectory_iterations::Integer,
    accuracy::AbstractFloat,maximum_hops_cutoff::Integer,separation_cutoff::Number;
    phi_limits=round.(100.*constants.hbar./bath_cutoff_energies),phi_steps=phi_limits./10000,
    E_limits=round.(6.*disorders),E_steps=E_limits./10000)

Produces the results for dKMC simulations of charge separation. To do so, it collects (or calculates) kappas 
(renormalisation constant), K_tots (Fourier transform of bath correlation function), and Hamiltonian and hopping 
radii. It then runs iterate_dKMC_charge_separation() on many realisations of disorder, parallelising the work amongst 
all of the available processes. Finally, from each realisation of disorder it collects and averages the mean outcomes, 
separation times and initial and final state characteristics.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorders`: Vector of energetic disorders (in meV) of electrons in the acceptor and holes in the donor.
- `donor_HOMO_acceptor_LUMO_gap`: Energy gap between donor HOMO and acceptor LUMO levels (in meV).
- `electronic_couplings`: Vector of nearest neighbour electronic couplings (in meV) for electrons in the acceptor and holes in the donor.
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `CT_lifetime`: CT state lifetime (in s).
- `bath_reorganisation_energies`: Vector of reorganisation energies of the bath for electrons and holes (in meV).  
- `bath_cutoff_energies`: Vector of cutoff energies of the bath for electrons and holes (in meV).
- `T`: Temperature (in K).
- `landscape_iterations`: Number of realisations of disorder (or energetic landscapes) for simulations to be run on.
- `trajectory_iterations`: Number of trajectories to be simulated on each realisation of disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `maximum_hops_cutoff`: Maximum number of hops the simulation will run before terminating.
- `separation_cutoff`: Distance (in m) at which charges are considered separated, and the simulation terminates. 

# Optional arguments:
- `phi_steps`: Vector of the time steps used in calculating phi for electrons and holes.
- `phi_limits`: Vector of the time limits used in calculating phi for electrons and holes.
- `E_steps`: Vector of the energy step used for values that K_tot is calculated at for electrons and holes.
- `E_limits`: Vector of the energy limits used for values that K_tot is calculated at for electrons and holes.

# Output:
- `mean_outcomes`: Vector of percentage of trajectories that end in each of the four possible outcomes.
- `mean_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in m) of the intial state.
- `mean_separation_time`: The mean separation time (in s).
- `mean_separated_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in m) of the intial state for trajectories where charges separate.
- `mean_separated_final_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in m) of the final state for trajectories where charges separate.
    
"""
function dKMC_charge_separation_results(dimension::Integer,N::Integer,disorders::Vector{<:Number},donor_HOMO_acceptor_LUMO_gap::Number,electronic_couplings::Vector{<:Number},epsilon_r::Number,site_spacing::Number,CT_lifetime::Number,bath_reorganisation_energies::Vector{<:Number},bath_cutoff_energies::Vector{<:Number},T::Number,landscape_iterations::Integer,trajectory_iterations::Integer,accuracy::AbstractFloat,maximum_hops_cutoff::Integer,separation_cutoff::Number;phi_limits=round.(100 .* constants.hbar./bath_cutoff_energies),phi_steps=phi_limits./10000,E_limits=round.(6 .* disorders),E_steps=E_limits./10000)

    #Create the spectral density function, here a super-ohmic spectral density is used.
    electron_J(E) = (bath_reorganisation_energies[1]/2) * (E/bath_cutoff_energies[1])^3 * exp(-E/bath_cutoff_energies[1])
    hole_J(E) = (bath_reorganisation_energies[2]/2) * (E/bath_cutoff_energies[2])^3 * exp(-E/bath_cutoff_energies[2])
            
    #Calculating the renormalisation constant (kappa) for the electronic couplings following polaron transformation.
    electron_kappa = sPTRE.calculate_kappa(electron_J,T)
    hole_kappa = sPTRE.calculate_kappa(hole_J,T)
    kappas = [electron_kappa, hole_kappa]

    #Reading the precomputed K values, or calculating them if not. K values are stored in a vector and called upon by looking at E and lambda+3 as the indices.
    electron_K_tot = sPTRE.collect_K_results(electron_J,E_steps[1],E_limits[1],phi_steps[1],phi_limits[1],bath_reorganisation_energies[1],bath_cutoff_energies[1],T)
    hole_K_tot = sPTRE.collect_K_results(hole_J,E_steps[2],E_limits[2],phi_steps[2],phi_limits[2],bath_reorganisation_energies[2],bath_cutoff_energies[2],T)
    K_tots = [electron_K_tot, hole_K_tot]

    #Calculating the approximating hamiltonian and hopping radii for each particle, or reading them from a saved file if previously calculated.
    acceptor_electron_hamiltonian_radius,acceptor_electron_hopping_radius = approximating_radii_functions.collect_charge_approximating_radii(dimension,N,disorders[1],electronic_couplings[1],bath_reorganisation_energies[1],bath_cutoff_energies[1],T,kappas[1],accuracy,landscape_iterations,K_tots[1],E_steps[1],E_limits[1])
    donor_hole_hamiltonian_radius,donor_hole_hopping_radius = approximating_radii_functions.collect_charge_approximating_radii(dimension,N,disorders[2],electronic_couplings[2],bath_reorganisation_energies[2],bath_cutoff_energies[2],T,kappas[2],accuracy,landscape_iterations,K_tots[2],E_steps[2],E_limits[2])
    hamiltonian_radii = [acceptor_electron_hamiltonian_radius, donor_hole_hamiltonian_radius]
    hopping_radii = [acceptor_electron_hopping_radius, donor_hole_hopping_radius]

    #We parallelise this function over many processes, to repeat it for many (landscape_iterations) realisations of energetic disorder. 
    results = pmap((realisations_of_disorder)->dKMC_charge_separation_functions.iterate_dKMC_charge_separation(realisations_of_disorder,N,disorders,donor_HOMO_acceptor_LUMO_gap,electronic_couplings,epsilon_r,site_spacing,CT_lifetime,bath_reorganisation_energies,kappas,trajectory_iterations,accuracy,maximum_hops_cutoff,separation_cutoff,hopping_radii,hamiltonian_radii,K_tots,E_steps,E_limits), Int.(dimension.*ones(landscape_iterations)))

    #Average the results over all landscape iterations.
    outcomes = []
    initial_state_characteristics = []
    separation_times = []
    separated_initial_state_characteristics = []
    separated_final_state_characteristics = []
    for i=1:landscape_iterations
        push!(outcomes,results[i][1])
        push!(initial_state_characteristics,results[i][2])
        if results[i][3] != 0 
            push!(separation_times,results[i][3])
            push!(separated_initial_state_characteristics,results[i][4])
            push!(separated_final_state_characteristics,results[i][5])
        end
    end
    mean_outcomes = [mean(outcomes); std(outcomes)/sqrt(landscape_iterations)]
    mean_initial_state_characteristics = [mean(initial_state_characteristics); std(initial_state_characteristics)/sqrt(landscape_iterations)]
    if isempty(separation_times)
        mean_separation_time = 0
        mean_separated_initial_state_characteristics = 0
        mean_separated_final_state_characteristics = 0
    else
        mean_separation_time = [mean(separation_times); std(separation_times)/sqrt(length(separation_times))]
        mean_separated_initial_state_characteristics = [mean(separated_initial_state_characteristics); std(separated_initial_state_characteristics)./sqrt(length(separated_initial_state_characteristics))]
        mean_separated_final_state_characteristics = [mean(separated_final_state_characteristics); std(separated_final_state_characteristics)./sqrt(length(separated_final_state_characteristics))]
    end

    #Print warning message if too many trajectories exceeded the maximum hops cutoff.
    if mean_outcomes[1,3] > 0.05
        @warn "A large proportion (>5%) of trajectories exceeded the maximum hops cutoff. Consider increasing the cutoff and repeating the calculation."
    end    

    return 	mean_outcomes, mean_initial_state_characteristics, mean_separation_time, mean_separated_initial_state_characteristics, mean_separated_final_state_characteristics

end


"""
    iterate_dKMC_charge_separation(dimension::Integer,N::Integer,disorders::Vector{<:Number},
    donor_HOMO_acceptor_LUMO_gap::Number,electronic_couplings::Vector{<:Number},epsilon_r::Number,
    site_spacing::Number,CT_lifetime::Number,bath_reorganisation_energies::Vector{<:Number},
    kappas::Vector{<:AbstractFloat},trajectory_iterations::Integer,accuracy::AbstractFloat,
    maximum_hops_cutoff::Integer,separation_cutoff::Number,hopping_radii::Vector{<:AbstractFloat},
    hamiltonian_radii::Vector{<:AbstractFloat},K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:Number},
    E_limits::Vector{<:Number})	
    
Repeats dKMC_charge_separation() for many trajectories on a single realisation of energetic disorder, before averaging 
outcomes, initial state characteristics (energy, IPR and separation), separation time, and characteristics of initial 
and final states of trajectories that successfully separated.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorders`: Vector of energetic disorders (in meV) of electrons in the acceptor and holes in the donor.
- `donor_HOMO_acceptor_LUMO_gap`: Energy gap between donor HOMO and acceptor LUMO levels (in meV).
- `electronic_couplings`: Vector of nearest neighbour electronic couplings (in meV) for electrons in the acceptor and holes in the donor.
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `CT_lifetime`: CT state lifetime (in s).
- `bath_reorganisation_energies`: Vector of reorganisation energies of the bath for electrons and holes (in meV).  
- `kappas`: Vector containing renormalisation constants for electron and hole couplings.
- `trajectory_iterations`: Number of trajectories to be simulated on each realisation of disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `maximum_hops_cutoff`: Maximum number of hops the simulation will run before terminating.
- `separation_cutoff`: Distance (in m) at which charges are considered separated, and the simulation terminates. 
- `hopping_radii`: Vector of precalculated distances that electrons and holes can hop. 
- `hamiltonian_radii`: Vector of precalculated distance for sites that are included in Hamiltonian subset.
- `K_tots`: Vector of precaclulated K_tot values for electrons and holes required for dKMC rate calculations.
- `E_steps`: Vector of the energy step used for values that K_tot is calculated at for electrons and holes.
- `E_limits`: Vector of the energy limits used for values that K_tot is calculated at for electrons and holes.

# Output:
- `mean_outcomes`: Vector of percentage of trajectories that end in each of the four possible outcomes.
- `mean_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in m) of the intial state.
- `mean_separation_time`: The mean separation time (in s).
- `mean_separated_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in m) of the intial state for trajectories where charges separate.
- `mean_separated_final_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in m) of the final state for trajectories where charges separate.

"""
function iterate_dKMC_charge_separation(dimension::Integer,N::Integer,disorders::Vector{<:Number},donor_HOMO_acceptor_LUMO_gap::Number,electronic_couplings::Vector{<:Number},epsilon_r::Number,site_spacing::Number,CT_lifetime::Number,bath_reorganisation_energies::Vector{<:Number},kappas::Vector{<:AbstractFloat},trajectory_iterations::Integer,accuracy::AbstractFloat,maximum_hops_cutoff::Integer,separation_cutoff::Number,hopping_radii::Vector{<:AbstractFloat},hamiltonian_radii::Vector{<:AbstractFloat},K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:Number},E_limits::Vector{<:Number})

    #Calculate a list of energies of the electron and hole sites.
    acceptor_electron_energies = randn(Int(N^dimension/2))*disorders[1] .+ donor_HOMO_acceptor_LUMO_gap
    donor_hole_energies = randn(Int(N^dimension/2))*disorders[2]
    energies = vcat(donor_hole_energies,acceptor_electron_energies)

    #We want to store how many trajectories conclude with each possible outcome.
    outcomes = zeros(1,4)
    
    #Create empty vectors for storing information about the initial state, separation time, and initial and final states for trajectories that separate.
    initial_states = []
    seperation_times = []
    separated_initial_states = []
    separated_final_states = []

    #Looping over the same landscape a total of trajectory_iterations amount of times.
    for iters = 1:trajectory_iterations

        #Performing the dKMC procedure for a single trajectory on a single realisation of energetic disorder.
        outcome,t,initial_state_characteristics,final_state_characteristics = dKMC_charge_separation(dimension,N,energies,electronic_couplings,epsilon_r,site_spacing,CT_lifetime,bath_reorganisation_energies,kappas,accuracy,maximum_hops_cutoff,separation_cutoff,hopping_radii,hamiltonian_radii,K_tots,E_steps,E_limits)
        
        #Storing outcomes and trajectory information.
        if outcome == "separation"
            outcomes[1] += 1
            push!(seperation_times,t)
            push!(separated_initial_states,initial_state_characteristics)
            push!(separated_final_states,final_state_characteristics)
        elseif outcome == "recombination"
            outcomes[2] += 1
        elseif outcome == "exceeded maximum hops cutoff"
            outcomes[3] +=1
        elseif outcome == "reached the edge of the system"
            outcomes[4] += 1
        end
        push!(initial_states,initial_state_characteristics)

    end

    #Averaging outcomes over all of the trajectories.
    mean_outcomes = outcomes./trajectory_iterations
    
    #Average initial state characteristics over all trajctories.
    mean_initial_state_characteristics = mean(initial_states)

    #Average separation characteristics over trajectories that underwent separation.
    if isempty(seperation_times)
        mean_separation_time = 0
        mean_separated_initial_state_characteristics = 0
        mean_separated_final_state_characteristics = 0
    else
        mean_separation_time = mean(seperation_times)
        mean_separated_initial_state_characteristics = mean(separated_initial_states)
        mean_separated_final_state_characteristics = mean(separated_final_states)
    end

    return mean_outcomes, mean_initial_state_characteristics, mean_separation_time, mean_separated_initial_state_characteristics, mean_separated_final_state_characteristics

end


"""
    dKMC_charge_separation(dimension::Integer,N::Integer,energies::Vector{<:AbstractFloat},
    electronic_couplings::Vector{<:Number},epsilon_r::Number,site_spacing::Number,CT_lifetime::Number,
    bath_reorganisation_energies::Vector{<:Number},kappas::Vector{<:AbstractFloat},accuracy::AbstractFloat,
    maximum_hops_cutoff::Integer,separation_cutoff::Number,hopping_radii::Vector{<:AbstractFloat},
    hamiltonian_radii::Vector{<:AbstractFloat},K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:Number},
    E_limits::Vector{<:Number})

Propagates the dKMC procedure for a single charge separation trajectory. It starts with a CT state closest to the 
centre of the landscape, where electrons are restricted to the acceptor and holes to the donor. Using dKMC rates it 
propagates the separation dynamics until a) the charges are separated by greater than the concluding separation, b) 
they recombine from interfacial CT states or c) the number of hops exceeds maximum_hops_cutoff. It returns the outcome 
(separated, recombined, or over hopped), the elapsed time, and the characteritics (energy, IPR, and separation) of the 
initial and final state.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `energies`: List of site energies of every site, of length N^dimension and containing donor HOMO energies followed by acceptor LUMO energies.
- `electronic_couplings`: Vector of nearest neighbour electronic couplings (in meV) for electrons in the acceptor and holes in the donor.
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `CT_lifetime`: CT state lifetime (in s).
- `bath_reorganisation_energies`: Vector of reorganisation energies of the bath for electrons and holes (in meV).  
- `kappas`: Vector containing renormalisation constants for electron and hole couplings.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `maximum_hops_cutoff`: Maximum number of hops the simulation will run before terminating.
- `separation_cutoff`: Distance (in m) at which charges are considered separated, and the simulation terminates. 
- `hopping_radii`: Vector of precalculated distances that electrons and holes can hop. 
- `hamiltonian_radii`: Vector of precalculated distance for sites that are included in Hamiltonian subset.
- `K_tots`: Vector of precaclulated K_tot values for electrons and holes required for dKMC rate calculations.
- `E_steps`: Vector of the energy step used for values that K_tot is calculated at for electrons and holes.
- `E_limits`: Vector of the energy limits used for values that K_tot is calculated at for electrons and holes.
    
# Output:
- `outcome`: String describing the trajectories outcome. 
- `t`: Elapsed time (in s).
- `initial_state_characteristics`: Vector containing the energy (in meV), IPR, and separation (in nm) of the initial state.
- `final_state_characteristics`: Vector containing the energy (in meV), IPR, and separation (in nm) of the final state.

"""
function dKMC_charge_separation(dimension::Integer,N::Integer,energies::Vector{<:AbstractFloat},electronic_couplings::Vector{<:Number},epsilon_r::Number,site_spacing::Number,CT_lifetime::Number,bath_reorganisation_energies::Vector{<:Number},kappas::Vector{<:AbstractFloat},accuracy::AbstractFloat,maximum_hops_cutoff::Integer,separation_cutoff::Number,hopping_radii::Vector{<:AbstractFloat},hamiltonian_radii::Vector{<:AbstractFloat},K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:Number},E_limits::Vector{<:Number})

    #Defining the initial position of the electron and hole, an interfacial CT state in the middle of the lattice.	
    current_hole_location = ones(dimension).* (N/2)
    current_electron_location = copy(current_hole_location)
    current_electron_location[1] += 1

    #Calculating the current Hamiltonian, the polaron transformed Hamiltonian, the reference matrices.
    H,Ht,electron_r,hole_r,transformed_coupling,electron_index,hole_index,bath_index,site_pair_indexes = setup_hamiltonian.current_charge_separation_hamiltonian(dimension,N,energies,electronic_couplings,epsilon_r,site_spacing,bath_reorganisation_energies,kappas,hamiltonian_radii,current_electron_location,current_hole_location)

    #Finding the indexs corresponding to interfacial CT site-pairs.
    interfacial_sites = findall(x->x==1,sum([(electron_r[:,i] .- hole_r[:,i]).^2 for i in 1:dimension]))

    #Calculating eigenvectors and eigenvalues in the energy eigenbasis.
    evals,evecs = eigen(Ht)

    #Calculating the expectation value of positions of energy eigenstates.
    electron_centres = [[evecs[:,i]' * Diagonal(electron_r[:,j]) * evecs[:,i] for i=eachindex(evals)] for j in 1:dimension]
    hole_centres = [[evecs[:,i]' * Diagonal(hole_r[:,j]) * evecs[:,i] for i=eachindex(evals)] for j in 1:dimension]

    #Choose initial CT state closest to middle of system.
    current_state = argmin(sqrt.(sum([(electron_centres[i] .- current_electron_location[i]).^2 for i in 1:dimension])) .+ sqrt.(sum([(hole_centres[i] .- current_hole_location[i]).^2 for i in 1:dimension])))

    #Recording characteristics of the current state.
    current_electron_location = [electron_centres[i][current_state] for i in 1:dimension]
    current_hole_location = [hole_centres[i][current_state] for i in 1:dimension]
    current_energy = evals[current_state]
    current_IPR = 1/sum(evecs[:,current_state].^4)
    current_separation = sqrt(sum((current_electron_location .- current_hole_location).^2)) * site_spacing
    
    #Recording characteristics of the initial state.
    initial_state_characteristics = [current_energy current_IPR constants.to_nano(current_separation)] 

    #Setting the initial conditions.
    t = 0
    hops = 0
    outcome = ""

    #Iterate the dKMC procedure until one of the terminating conditions.
    while hops < maximum_hops_cutoff

        #Checking if either charge is too close to the edge of the landscape to create a new subset of the Hamiltonian.
        if iszero(current_electron_location .+ hamiltonian_radii[2] .> N) == false || iszero(current_hole_location .- hamiltonian_radii[1] .< 1) == false 
            @warn "Simulation terminated as a charge carrier got too close to the edge of the system. Increase the length of the system (N) to avoid this error." 
            outcome = "reached the edge of the system"
            break
        end

        #Rediagonalise a new subset of the Hamiltonian.
        previous_eigenstate = evecs[:,current_state]
        previous_site_pair_indexes = copy(site_pair_indexes)
        H,Ht,electron_r,hole_r,transformed_coupling,electron_index,hole_index,bath_index,site_pair_indexes = setup_hamiltonian.current_charge_separation_hamiltonian(dimension,N,energies,electronic_couplings,epsilon_r,site_spacing,bath_reorganisation_energies,kappas,hamiltonian_radii,current_electron_location,current_hole_location)
        interfacial_sites = findall(x->x==1,sum([(electron_r[:,i] .- hole_r[:,i]).^2 for i in 1:dimension]))
        evals,evecs = eigen(Ht)
        electron_centres = [[evecs[:,i]' * Diagonal(electron_r[:,j]) * evecs[:,i] for i=eachindex(evals)] for j in 1:dimension]
        hole_centres = [[evecs[:,i]' * Diagonal(hole_r[:,j]) * evecs[:,i] for i=eachindex(evals)] for j in 1:dimension]
        current_state = argmax(((previous_eigenstate[findall(x->x in site_pair_indexes,previous_site_pair_indexes)]' * evecs[findall(x->x in previous_site_pair_indexes,site_pair_indexes),:])[:]).^2)
        current_electron_location = [electron_centres[i][current_state] for i in 1:dimension]
        current_hole_location = [hole_centres[i][current_state] for i in 1:dimension]

        #Finding which states are accessible from the current state.
        accessible_states = findall(x-> 0 < x < 1,(sum([(electron_centres[i] .- current_electron_location[i]).^2 for i in 1:dimension]))./(hopping_radii[1]^2) .+ (sum([(hole_centres[i] .- current_hole_location[i]).^2 for i in 1:dimension]))./(hopping_radii[2]^2))

        #Calculating hopping rates to all states in accessible_states.
        hopping_rates = zeros(length(accessible_states))
        current_state_relevant_sites = sortperm(abs.(evecs[:,current_state]),rev=true)[1:findfirst(x->x>accuracy,cumsum(sort(abs.(evecs[:,current_state])./sum(abs.(evecs[:,current_state])),rev=true)))]
        for f = eachindex(accessible_states)
            destination_state_relevant_sites = sortperm(abs.(evecs[:,accessible_states[f]]),rev=true)[1:findfirst(x->x>accuracy,cumsum(sort(abs.(evecs[:,accessible_states[f]])./sum(abs.(evecs[:,accessible_states[f]])),rev=true)))]
            hopping_rate = dKMC_hopping_rates.charge_separation_dKMC_rate(current_state,accessible_states[f],transformed_coupling,evals,evecs,K_tots,E_steps,E_limits,current_state_relevant_sites,destination_state_relevant_sites,electron_index,hole_index,bath_index)
            if hopping_rate > 0
                hopping_rates[f] = hopping_rate
            end
        end

        #Calculating the recombinate rate.
        push!(hopping_rates,(1/constants.to_femto(CT_lifetime))*(sum(evecs[:,current_state][interfacial_sites]).^2))
        
        #Calculating sum of all rates.
        hopping_rate_sum = sum(hopping_rates)

        #Prevents an error if the rate sum comes out to zero.
        if hopping_rate_sum == 0
            break
        end

        #Selected hop is chosen propabalistically in proportion to the hopping rate. 
        ran = rand()
        index = findfirst(x->x>ran*hopping_rate_sum,cumsum(hopping_rates))

        #Time updated by calculating elapsed time.
        ran2 = rand()
        t += constants.from_femto(hopping_rate_sum^(-1)*log(1/ran2))
        hops += 1

        #Checking if recombined, if so we break.
        if index == length(hopping_rates) 
            outcome = "recombination"
            break
        end

        #Updating the current state.
        current_state = accessible_states[index]

        #Recording the characteristics of the new state.
        current_electron_location = [electron_centres[i][current_state] for i in 1:dimension]
        current_hole_location = [hole_centres[i][current_state] for i in 1:dimension]
        current_energy = evals[current_state]
        current_IPR = 1/sum(evecs[:,current_state].^4)
        current_separation = sqrt(sum((current_electron_location .- current_hole_location).^2)) * site_spacing

        #Checking if charges are separated by greater than separation_cutoff, if so we break.
        if current_separation > separation_cutoff
            outcome = "separation"
            break
        end
    end

    #Checking whether we have had more hops than maximum_hops_cutoff.
    if hops == maximum_hops_cutoff
        outcome = "exceeded maximum hops cutoff"
    end

    #Saving the characteristics of the final state.
    final_state_characteristics = [current_energy current_IPR constants.to_nano(current_separation)]
    
    return outcome, t, initial_state_characteristics, final_state_characteristics

end

end
