module dKMC_charge_generation_functions

#Loading the required packages.
include("../shared_functions/package_loading.jl")
include("../shared_functions/constants.jl")
include("../shared_functions/sPTRE.jl")
include("../shared_functions/setup_hamiltonian.jl")
include("../shared_functions/dKMC_hopping_rates.jl")
include("../shared_functions/approximating_radii_functions.jl")


"""
    dKMC_charge_generation_results(dimension::Integer,N::Integer,disorders::Matrix{<:Number},
    exciton_disorders::Vector{<:Number},donor_HOMO_LUMO_gap::Number,LUMO_offset::Number,HOMO_offset::Number,
    exciton_binding_energies::Vector{<:Number},electronic_couplings::Matrix{<:Number},
    transition_dipole_moments::Vector{<:Number},epsilon_r::Number,site_spacing::Number,
    exciton_lifetimes::Vector{<:Number},CT_lifetimes::Vector{<:Number},
    bath_reorganisation_energies::Vector{<:Number},bath_cutoff_energies::Vector{<:Number},T::Number,
    landscape_iterations::Integer,trajectory_iterations::Integer,accuracy::AbstractFloat,
    maximum_hops_cutoff::Integer,separation_cutoff::Number,exciton_population_cutoff::AbstractFloat,
    maximum_excitation_distance::Number;
    phi_limits=round.(100 .* constants.hbar./bath_cutoff_energies),phi_steps=phi_limits./10000,
    E_limits=round.(6 .* vcat(maximum(disorders,dims=2)[:],maximum(exciton_disorders)) .+ [LUMO_offset,HOMO_offset,HOMO_offset-LUMO_offset]),
    E_steps = round.(6 .* vcat(minimum(disorders,dims=2)[:],minimum(exciton_disorders)))./10000)

Produces the results for dKMC simulations of charge generation. To do so, it collects (or calculates) kappas 
(renormalisation constant), K_tots (Fourier transform of bath correlation function), and Hamiltonian and hopping 
radii. It then runs iterate_dKMC_charge_generation() on many realisations of disorder, parallelising the work amongst 
all of the available processes. Finally, from each realisation of disorder it collects and averages the mean outcomes, 
separation times and initial and final state characteristics.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorders`: Matrix of energetic disorders (in meV) of the electrons (row 1) and holes (row 2) in the donor (column 1) and acceptor (column 2).
- `exciton_disorders`: Vector of excitonic disorders of the donor and acceptor (in meV).
- `donor_HOMO_LUMO_gap`: Energy gap between donor HOMO and LUMO levels (in meV).
- `LUMO_offset`: Energy offset between donor and acceptor LUMO levels (in meV).
- `HOMO_offset`: Energy offset between donor and acceptor HOMO levels (in meV).
- `exciton_binding_energies`: Vector of exciton binding energies for the donor and acceptor (in meV).
- `electronic_couplings`: Matrix of nearest neighbour electronic couplings (in meV) for electrons (row 1) and holes (row 2) in the donor (column 1), at the interface (column 2), and in the acceptor (column 3).
- `transition_dipole_moments`: Vector of magnitude of the transition dipole moments on donor and acceptor sites (in D).
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `exciton_lifetimes`: Vector of exciton lifetimes (in s) in the donor and acceptor.
- `CT_lifetimes`: Vector of CT state lifetimes (in s) in the donor, at the interface, and in the acceptor.
- `bath_reorganisation_energies`: Vector of reorganisation energies of the bath for electrons, holes, and excitons (in meV).  
- `bath_cutoff_energies`: Vector of cutoff energies of the bath for electrons, holes, and excitons (in meV).
- `T`: Temperature (in K).
- `landscape_iterations`: Number of realisations of disorder (or energetic landscapes) for simulations to be run on.
- `trajectory_iterations`: Number of trajectories to be simulated on each realisation of disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `maximum_hops_cutoff`: Maximum number of hops the simulation will run before terminating.
- `separation_cutoff`: Distance (in m) at which charges are considered separated, and the simulation terminates. 
- `exciton_population_cutoff`: Minimum population on exciton site-pairs to treat the current state as an exciton.
- `maximum_excitation_distance`: Maximum distance (in m) that an exciton can be excited from the interface.

# Optional arguments:
- `phi_steps`: Vector of the time steps used in calculating phi for electrons, holes, and excitons.
- `phi_limits`: Vector of the time limits used in calculating phi for electrons, holes, and excitons.
- `E_steps`: Vector of the energy steps used for values that K_tot is calculated at for electrons, holes, and excitons.
- `E_limits`: Vector of the energy limits used for values that K_tot is calculated at for electrons, holes, and excitons.

# Output:
- `mean_outcomes`: Vector of percentage of trajectories that end in each of the seven possible outcomes.
- `mean_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the intial state.
- `mean_separation_time`: The mean separation time (in s).
- `mean_separated_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the intial state for trajectories where charges separate.
- `mean_separated_final_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the final state for trajectories where charges separate.
- `mean_interfacial_separation_time`: The mean separation time (in s) for trajectories underdoing interfacial separation.
- `mean_interfacial_separated_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the intial state for trajectories underdoing interfacial separation.
- `mean_interfacial_separated_final_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the final state for trajectories underdoing interfacial separation.
- `mean_bulk_separation_time`: The mean separation time (in s) for trajectories underdoing bulk separation.
- `mean_bulk_separated_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the intial state for trajectories underdoing bulk separation.
- `mean_bulk_separated_final_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the final state for trajectories underdoing bulk separation.

"""
function dKMC_charge_generation_results(dimension::Integer,N::Integer,disorders::Matrix{<:Number},exciton_disorders::Vector{<:Number},donor_HOMO_LUMO_gap::Number,LUMO_offset::Number,HOMO_offset::Number,exciton_binding_energies::Vector{<:Number},electronic_couplings::Matrix{<:Number},transition_dipole_moments::Vector{<:Number},epsilon_r::Number,site_spacing::Number,exciton_lifetimes::Vector{<:Number},CT_lifetimes::Vector{<:Number},bath_reorganisation_energies::Vector{<:Number},bath_cutoff_energies::Vector{<:Number},T::Number,landscape_iterations::Integer,trajectory_iterations::Integer,accuracy::AbstractFloat,maximum_hops_cutoff::Integer,separation_cutoff::Number,exciton_population_cutoff::AbstractFloat,maximum_excitation_distance::Number;phi_limits=round.(100 .* constants.hbar./bath_cutoff_energies),phi_steps=phi_limits./10000,E_limits=round.(6 .* vcat(maximum(disorders,dims=2)[:],maximum(exciton_disorders)) .+ [LUMO_offset,HOMO_offset,HOMO_offset-LUMO_offset]),E_steps=round.(6 .* vcat(minimum(disorders,dims=2)[:],minimum(exciton_disorders)))./10000)

    #Create the spectral density function for each particle type, here a super-ohmic spectral density is used.
    electron_J(E) = (bath_reorganisation_energies[1]/2) * (E/bath_cutoff_energies[1])^3 * exp(-E/bath_cutoff_energies[1])
    hole_J(E) = (bath_reorganisation_energies[2]/2) * (E/bath_cutoff_energies[2])^3 * exp(-E/bath_cutoff_energies[2])
    exciton_J(E) = (bath_reorganisation_energies[3]/2) * (E/bath_cutoff_energies[3])^3 * exp(-E/bath_cutoff_energies[3])

    #Calculating expressions for the mixed spectral densities for cross terms in the damping rates.
    electron_and_hole_J(E) = sqrt(electron_J(E)*hole_J(E))
    electron_and_exciton_J(E) = sqrt(electron_J(E)*exciton_J(E))
    hole_and_exciton_J(E) = sqrt(hole_J(E)*exciton_J(E))

    #Calculating for each particle type a renormalisation constant (kappa) for the couplings following polaron transformation.
    electron_kappa = sPTRE.calculate_kappa(electron_J,T)
    hole_kappa = sPTRE.calculate_kappa(hole_J,T)
    exciton_kappa = sPTRE.calculate_kappa(exciton_J,T)
    kappas = [electron_kappa, hole_kappa, exciton_kappa]

    #Reading the precomputed K values for each particle, or calculating them if not already computed. K values are stored in a vector and called upon by looking at E and lambda+3 as the indices.
    electron_K_tot = sPTRE.collect_K_results(electron_J,E_steps[1],E_limits[1],phi_steps[1],phi_limits[1],bath_reorganisation_energies[1],bath_cutoff_energies[1],T)
    hole_K_tot = sPTRE.collect_K_results(hole_J,E_steps[2],E_limits[2],phi_steps[2],phi_limits[2],bath_reorganisation_energies[2],bath_cutoff_energies[2],T)
    exciton_K_tot = sPTRE.collect_K_results(hole_J,E_steps[3],E_limits[3],phi_steps[3],phi_limits[3],bath_reorganisation_energies[3],bath_cutoff_energies[3],T)
    electron_and_hole_K_tot = sPTRE.collect_K_results(electron_and_hole_J,min(E_steps[1],E_steps[2]),max(E_limits[1],E_limits[2]),min(phi_steps[1],phi_steps[2]),max(phi_limits[1],phi_limits[2]),[bath_reorganisation_energies[1] bath_reorganisation_energies[2]] ,[bath_cutoff_energies[1] bath_cutoff_energies[2]],T)
    electron_and_exciton_K_tot = sPTRE.collect_K_results(electron_and_exciton_J,min(E_steps[1],E_steps[3]),max(E_limits[1],E_limits[3]),min(phi_steps[1],phi_steps[3]),max(phi_limits[1],phi_limits[3]),[bath_reorganisation_energies[1] bath_reorganisation_energies[3]] ,[bath_cutoff_energies[1] bath_cutoff_energies[3]],T)
    hole_and_exciton_K_tot = sPTRE.collect_K_results(hole_and_exciton_J,min(E_steps[2],E_steps[3]),max(E_limits[2],E_limits[3]),min(phi_steps[2],phi_steps[3]),max(phi_limits[2],phi_limits[3]),[bath_reorganisation_energies[2] bath_reorganisation_energies[3]] ,[bath_cutoff_energies[2] bath_cutoff_energies[3]],T)
    K_tots = [electron_K_tot, hole_K_tot, exciton_K_tot, electron_and_hole_K_tot, electron_and_exciton_K_tot, hole_and_exciton_K_tot]

    #Calculating the approximating hamiltonian and hopping radii for each particle and each material, or reading them from a saved file if previously calculated.
    donor_electron_hamiltonian_radius,donor_electron_hopping_radius = approximating_radii_functions.collect_charge_approximating_radii(dimension,N,disorders[1,1],electronic_couplings[1,1],bath_reorganisation_energies[1],bath_cutoff_energies[1],T,kappas[1],accuracy,landscape_iterations,K_tots[1],E_steps[1],E_limits[1])
    acceptor_electron_hamiltonian_radius,acceptor_electron_hopping_radius = approximating_radii_functions.collect_charge_approximating_radii(dimension,N,disorders[1,2],electronic_couplings[1,2],bath_reorganisation_energies[1],bath_cutoff_energies[1],T,kappas[1],accuracy,landscape_iterations,K_tots[1],E_steps[1],E_limits[1])
    donor_hole_hamiltonian_radius,donor_hole_hopping_radius = approximating_radii_functions.collect_charge_approximating_radii(dimension,N,disorders[2,1],electronic_couplings[2,1],bath_reorganisation_energies[2],bath_cutoff_energies[2],T,kappas[2],accuracy,landscape_iterations,K_tots[2],E_steps[2],E_limits[2])
    acceptor_hole_hamiltonian_radius,acceptor_hole_hopping_radius = approximating_radii_functions.collect_charge_approximating_radii(dimension,N,disorders[2,2],electronic_couplings[2,2],bath_reorganisation_energies[2],bath_cutoff_energies[2],T,kappas[2],accuracy,landscape_iterations,K_tots[2],E_steps[2],E_limits[2])
    donor_exciton_hamiltonian_radius,donor_exciton_hopping_radius = approximating_radii_functions.collect_exciton_approximating_radii(dimension,N,exciton_disorders[1],transition_dipole_moments[1],epsilon_r,bath_reorganisation_energies[3],bath_cutoff_energies[3],T,kappas[3],site_spacing,accuracy,landscape_iterations,K_tots[3],E_steps[3],E_limits[3])
    acceptor_exciton_hamiltonian_radius,acceptor_exciton_hopping_radius = approximating_radii_functions.collect_exciton_approximating_radii(dimension,N,exciton_disorders[2],transition_dipole_moments[2],epsilon_r,bath_reorganisation_energies[3],bath_cutoff_energies[3],T,kappas[3],site_spacing,accuracy,landscape_iterations,K_tots[3],E_steps[3],E_limits[3])
    hamiltonian_radii = [donor_electron_hamiltonian_radius acceptor_electron_hamiltonian_radius; donor_hole_hamiltonian_radius acceptor_hole_hamiltonian_radius; donor_exciton_hamiltonian_radius acceptor_exciton_hamiltonian_radius]
    hopping_radii = [donor_electron_hopping_radius acceptor_electron_hopping_radius; donor_hole_hopping_radius acceptor_hole_hopping_radius; donor_exciton_hopping_radius acceptor_exciton_hopping_radius]

    #We parallelise this function over many processes, to repeat it for many (landscape_iterations) realisations of energetic disorder.
    results = pmap((realisations_of_disorder)->iterate_dKMC_charge_generation(realisations_of_disorder,N,disorders,exciton_disorders,donor_HOMO_LUMO_gap,LUMO_offset,HOMO_offset,exciton_binding_energies,electronic_couplings,kappas,transition_dipole_moments,epsilon_r,site_spacing,exciton_lifetimes,CT_lifetimes,bath_reorganisation_energies,trajectory_iterations,accuracy,maximum_hops_cutoff,separation_cutoff,exciton_population_cutoff,maximum_excitation_distance,hopping_radii,hamiltonian_radii,K_tots,E_steps,E_limits), Int.(dimension.*ones(landscape_iterations)))

    #Average the results over all landscape iterations.
    outcomes = []
    initial_state_characteristics = []
    separation_times = []
    separated_initial_state_characteristics = []
    separated_final_state_characteristics = []
    interfacial_separation_times = []
    interfacial_separated_initial_state_characteristics = []
    interfacial_separated_final_state_characteristics = []
    bulk_separation_times = []
    bulk_separated_initial_state_characteristics = []
    bulk_separated_final_state_characteristics = []
    for i=1:landscape_iterations
        push!(outcomes,results[i][1])
        push!(initial_state_characteristics,results[i][2])
        if results[i][3] != 0 
            push!(separation_times,results[i][3])
            push!(separated_initial_state_characteristics,results[i][4])
            push!(separated_final_state_characteristics,results[i][5])
        end
        if results[i][6] != 0 
            push!(interfacial_separation_times,results[i][6])
            push!(interfacial_separated_initial_state_characteristics,results[i][7])
            push!(interfacial_separated_final_state_characteristics,results[i][8])
        end
        if results[i][9] != 0 
            push!(bulk_separation_times,results[i][9])
            push!(bulk_separated_initial_state_characteristics,results[i][10])
            push!(bulk_separated_final_state_characteristics,results[i][11])
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
    if isempty(interfacial_separation_times)
        mean_interfacial_separation_time = 0
        mean_interfacial_separated_initial_state_characteristics = 0
        mean_interfacial_separated_final_state_characteristics = 0
    else
        mean_interfacial_separation_time = [mean(interfacial_separation_times); std(interfacial_separation_times)/sqrt(length(interfacial_separation_times))]
        mean_interfacial_separated_initial_state_characteristics = [mean(interfacial_separated_initial_state_characteristics); std(interfacial_separated_initial_state_characteristics)./sqrt(length(interfacial_separated_initial_state_characteristics))]
        mean_interfacial_separated_final_state_characteristics = [mean(interfacial_separated_final_state_characteristics); std(interfacial_separated_final_state_characteristics)./sqrt(length(interfacial_separated_final_state_characteristics))]
    end
    if isempty(bulk_separation_times)
        mean_bulk_separation_time = 0
        mean_bulk_separated_initial_state_characteristics = 0
        mean_bulk_separated_final_state_characteristics = 0
    else
        mean_bulk_separation_time = [mean(bulk_separation_times); std(bulk_separation_times)/sqrt(length(bulk_separation_times))]
        mean_bulk_separated_initial_state_characteristics = [mean(bulk_separated_initial_state_characteristics); std(bulk_separated_initial_state_characteristics)./sqrt(length(bulk_separated_initial_state_characteristics))]
        mean_bulk_separated_final_state_characteristics = [mean(bulk_separated_final_state_characteristics); std(bulk_separated_final_state_characteristics)./sqrt(length(bulk_separated_final_state_characteristics))]
    end

    #Print warning message if too many trajectories exceeded the maximum hops cutoff.
    if mean_outcomes[1,6] > 0.05
        @warn "A large proportion (>5%) of trajectories exceeded the maximum hops cutoff. Consider increasing the cutoff and repeating the calculation."
    end

    return 	mean_outcomes, mean_initial_state_characteristics, mean_separation_time, mean_separated_initial_state_characteristics, mean_separated_final_state_characteristics, mean_interfacial_separation_time, mean_interfacial_separated_initial_state_characteristics, mean_interfacial_separated_final_state_characteristics, mean_bulk_separation_time, mean_bulk_separated_initial_state_characteristics, mean_bulk_separated_final_state_characteristics

end


"""
    iterate_dKMC_charge_generation(dimension::Integer,N::Integer,disorders::Matrix{<:Number},
    exciton_disorders::Vector{<:Number},donor_HOMO_LUMO_gap::Number,LUMO_offset::Number,HOMO_offset::Number,
    exciton_binding_energies::Vector{<:Number},electronic_couplings::Matrix{<:Number},
    kappas::Vector{<:AbstractFloat},transition_dipole_moments::Vector{<:Number},epsilon_r::Number,
    site_spacing::Number,exciton_lifetimes::Vector{<:Number},CT_lifetimes::Vector{<:Number},
    bath_reorganisation_energies::Vector{<:Number},trajectory_iterations::Integer,accuracy::AbstractFloat,
    maximum_hops_cutoff::Integer,separation_cutoff::Number,exciton_population_cutoff::AbstractFloat,
    maximum_excitation_distance::Number,hopping_radii::Matrix{<:AbstractFloat},
    hamiltonian_radii::Matrix{<:AbstractFloat},K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:Number},
    E_limits::Vector{<:Number})

Repeats dKMC_charge_generation() for many trajectories on a single realisation of energetic disorder, before averaging 
outcomes, initial state characteristics (energy, IPR and separation), separation time, and characteristics of initial 
and final states of trajectories that successfully separated.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorders`: Matrix of energetic disorders (in meV) of the electrons (row 1) and holes (row 2) in the donor (column 1) and acceptor (column 2).
- `exciton_disorders`: Vector of excitonic disorders of the donor and acceptor (in meV).
- `donor_HOMO_LUMO_gap`: Energy gap between donor HOMO and LUMO levels (in meV).
- `LUMO_offset`: Energy offset between donor and acceptor LUMO levels (in meV).
- `HOMO_offset`: Energy offset between donor and acceptor HOMO levels (in meV).
- `exciton_binding_energies`: Vector of exciton binding energies for the donor and acceptor (in meV).
- `electronic_couplings`: Matrix of nearest neighbour electronic couplings (in meV) for electrons (row 1) and holes (row 2) in the donor (column 1), at the interface (column 2), and in the acceptor (column 3).
- `kappas`: Vector containing renormalisation constants for electron, hole, and exciton couplings.
- `transition_dipole_moments`: Vector of magnitude of the transition dipole moments on donor and acceptor sites (in D).
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `exciton_lifetimes`: Vector of exciton lifetimes (in s) in the donor and acceptor.
- `CT_lifetimes`: Vector of CT state lifetimes (in s) in the donor, at the interface, and in the acceptor.
- `bath_reorganisation_energies`: Vector of reorganisation energies of the bath for electrons, holes, and excitons (in meV).  
- `trajectory_iterations`: Number of trajectories to be simulated on each realisation of disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `maximum_hops_cutoff`: Maximum number of hops the simulation will run before terminating.
- `separation_cutoff`: Distance (in m) at which charges are considered separated, and the simulation terminates. 
- `exciton_population_cutoff`: Minimum population on exciton site-pairs to treat the current state as an exciton.
- `maximum_excitation_distance`: Maximum distance (in m) that an exciton can be excited from the interface.
- `hopping_radii`: Matrix of precalculated distances that electrons, holes, and excitons can hop in the donor (row 1) and acceptor (row 2). 
- `hamiltonian_radii`: Matrix of precalculated distance for electron, hole, and exciton sites in the donor (row 1) and acceptor (row 2) that are included as site-pairs in Hamiltonian subset.
- `K_tots`: Vector of precaclulated K_tot values for electrons, holes, and excitons required for dKMC rate calculations.
- `E_steps`: Vector of the energy step used for values that K_tot is calculated at for electrons, holes, and excitons.
- `E_limits`: Vector of the energy limits used for values that K_tot is calculated at for electrons, holes, and excitons.

# Output:
- `mean_outcomes`: Vector of percentage of trajectories that end in each of the seven possible outcomes.
- `mean_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the intial state.
- `mean_separation_time`: The mean separation time (in s).
- `mean_separated_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the intial state for trajectories where charges separate.
- `mean_separated_final_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the final state for trajectories where charges separate.
- `mean_interfacial_separation_time`: The mean separation time (in s) for trajectories underdoing interfacial separation.
- `mean_interfacial_separated_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the intial state for trajectories underdoing interfacial separation.
- `mean_interfacial_separated_final_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the final state for trajectories underdoing interfacial separation.
- `mean_bulk_separation_time`: The mean separation time (in s) for trajectories underdoing bulk separation.
- `mean_bulk_separated_initial_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the intial state for trajectories underdoing bulk separation.
- `mean_bulk_separated_final_state_characteristics`: Vector of the mean energy (in meV), IPR, and electron-hole separation (in nm) of the final state for trajectories underdoing bulk separation.

"""
function iterate_dKMC_charge_generation(dimension::Integer,N::Integer,disorders::Matrix{<:Number},exciton_disorders::Vector{<:Number},donor_HOMO_LUMO_gap::Number,LUMO_offset::Number,HOMO_offset::Number,exciton_binding_energies::Vector{<:Number},electronic_couplings::Matrix{<:Number},kappas::Vector{<:AbstractFloat},transition_dipole_moments::Vector{<:Number},epsilon_r::Number,site_spacing::Number,exciton_lifetimes::Vector{<:Number},CT_lifetimes::Vector{<:Number},bath_reorganisation_energies::Vector{<:Number},trajectory_iterations::Integer,accuracy::AbstractFloat,maximum_hops_cutoff::Integer,separation_cutoff::Number,exciton_population_cutoff::AbstractFloat,maximum_excitation_distance::Number,hopping_radii::Matrix{<:AbstractFloat},hamiltonian_radii::Matrix{<:AbstractFloat},K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:Number},E_limits::Vector{<:Number})

    #Assigning HOMO and LUMO energies to every site.
    LUMO_HOMO_energies = setup_hamiltonian.LUMO_HOMO_energies(dimension,N,disorders,exciton_disorders,donor_HOMO_LUMO_gap,HOMO_offset,LUMO_offset)

    #Assigning a random orientation for a dipole of exciton on each site.
    dipole_orientations = setup_hamiltonian.assign_exciton_dipole_orientations(dimension,N)

    #We want to store how many trajectories conclude with each of the 7 possible outcomes.
    outcomes = zeros(1,7)

    #Create empty vectors for storing information about the initial state, separation time, and initial and final states for trajectories that separate.
    initial_states = []
    separation_times = []
    separated_initial_states = []
    separated_final_states = []
    interfacial_separation_times = []
    interfacial_separated_initial_states = []
    interfacial_separated_final_states = []
    bulk_separation_times = []
    bulk_separated_initial_states = []
    bulk_separated_final_states = []

    #Looping over the same landscape a total of trajectory_iterations amount of times.
    for iterations = 1:trajectory_iterations

        #Performing the dKMC procedure for a single trajectory on a single realisation of energetic disorder.
        outcome, t, initial_state_characteristics, final_state_characteristics = dKMC_charge_generation(dimension,N,LUMO_HOMO_energies,exciton_binding_energies,electronic_couplings,kappas,transition_dipole_moments,dipole_orientations,epsilon_r,site_spacing,exciton_lifetimes,CT_lifetimes,bath_reorganisation_energies,accuracy,maximum_hops_cutoff,separation_cutoff,exciton_population_cutoff,maximum_excitation_distance,hopping_radii,hamiltonian_radii,K_tots,E_steps,E_limits)

        #Storing outcomes and trajectory information.
        if outcome == "interfacial separation"
            outcomes[1] += 1
            push!(separation_times,t)
            push!(separated_initial_states,initial_state_characteristics)
            push!(separated_final_states,final_state_characteristics)
            push!(interfacial_separation_times,t)
            push!(interfacial_separated_initial_states,initial_state_characteristics)
            push!(interfacial_separated_final_states,final_state_characteristics)
        elseif outcome == "bulk separation"
            outcomes[2] += 1
            push!(separation_times,t)
            push!(separated_initial_states,initial_state_characteristics)
            push!(separated_final_states,final_state_characteristics)
            push!(bulk_separation_times,t)
            push!(bulk_separated_initial_states,initial_state_characteristics)
            push!(bulk_separated_final_states,final_state_characteristics)
        elseif outcome == "exciton recombination"
            outcomes[3] += 1
        elseif outcome == "interfacial CT recombination"
            outcomes[4] += 1
        elseif outcome == "bulk CT recombination"
            outcomes[5] += 1
        elseif outcome == "exceeded maximum hops cutoff"
            outcomes[6] += 1
        elseif outcome == "reached the edge of the system"
            outcomes[7] += 1
        end
        push!(initial_states,initial_state_characteristics)

    end

    #Averaging outcomes over all of the trajectories.
    mean_outcomes = outcomes./trajectory_iterations

    #Average initial state characteristics over all trajctories.
    mean_initial_state_characteristics = mean(initial_states)

    #Average separation characteristics over trajectories that underwent separation.
    if isempty(separation_times)
        mean_separation_time = 0
        mean_separated_initial_state_characteristics = 0
        mean_separated_final_state_characteristics = 0
    else
        mean_separation_time = mean(separation_times)
        mean_separated_initial_state_characteristics = mean(separated_initial_states)
        mean_separated_final_state_characteristics = mean(separated_final_states)
    end

    #Average separation characteristics over trajectories that underwent interfacial separation.
    if isempty(interfacial_separation_times)
        mean_interfacial_separation_time = 0
        mean_interfacial_separated_initial_state_characteristics = 0
        mean_interfacial_separated_final_state_characteristics = 0
    else
        mean_interfacial_separation_time = mean(interfacial_separation_times)
        mean_interfacial_separated_initial_state_characteristics = mean(interfacial_separated_initial_states)
        mean_interfacial_separated_final_state_characteristics = mean(interfacial_separated_final_states)
    end
    
    #Average separation characteristics over trajectories that underwent bulk separation.
    if isempty(bulk_separation_times)
        mean_bulk_separation_time = 0
        mean_bulk_separated_initial_state_characteristics = 0
        mean_bulk_separated_final_state_characteristics = 0
    else
        mean_bulk_separation_time = mean(bulk_separation_times)
        mean_bulk_separated_initial_state_characteristics = mean(bulk_separated_initial_states)
        mean_bulk_separated_final_state_characteristics = mean(bulk_separated_final_states)
    end

    return mean_outcomes, mean_initial_state_characteristics,mean_separation_time, mean_separated_initial_state_characteristics, mean_separated_final_state_characteristics, mean_interfacial_separation_time, mean_interfacial_separated_initial_state_characteristics, mean_interfacial_separated_final_state_characteristics, mean_bulk_separation_time, mean_bulk_separated_initial_state_characteristics, mean_bulk_separated_final_state_characteristics

end


"""
    dKMC_charge_generation(dimension::Integer,N::Integer,LUMO_HOMO_energies::Matrix{<:AbstractFloat},
    exciton_binding_energies::Vector{<:Number},electronic_couplings::Matrix{<:Number},
    kappas::Vector{<:AbstractFloat},transition_dipole_moments::Vector{<:Number},
    dipole_orientations::Matrix{<:AbstractFloat},epsilon_r::Number,site_spacing::Number,
    exciton_lifetimes::Vector{<:Number},CT_lifetimes::Vector{<:Number},
    bath_reorganisation_energies::Vector{<:Number},accuracy::AbstractFloat,maximum_hops_cutoff::Integer,
    separation_cutoff::Number,exciton_population_cutoff::AbstractFloat,maximum_excitation_distance::Number,
    hopping_radii::Matrix{<:Number},hamiltonian_radii::Matrix{<:Number},K_tots::Vector{Matrix{ComplexF64}},
    E_steps::Vector{<:Number},E_limits::Vector{<:Number})

Propagates the dKMC procedure for a single charge generation trajectory. It starts with an exciton state in the donor 
material at a distance from the interface with the acceptor material of between 1 and maximum_excitation_distance 
sites. It then propagates the separation dynamics until the charges either recombine or separate. It returns the 
outcome of the trajectory, the time, and the intial and final state characterstics.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `LUMO_HOMO_energies`: Matrix of LUMO (row 1) and HOMO energies (row 2) on every site (in meV).
- `exciton_binding_energies`: Vector of exciton binding energies for the donor and acceptor (in meV).
- `electronic_couplings`: Matrix of nearest neighbour electronic couplings (in meV) for electrons (row 1) and holes (row 2) in the donor (column 1), at the interface (column 2), and in the acceptor (column 3).
- `kappas`: Vector containing renormalisation constants for electron, hole, and exciton couplings.
- `transition_dipole_moments`: Vector of magnitude of the transition dipole moments on donor and acceptor sites (in D).
- `dipole_orientations`: List of unit vectors describing the orientation of the transition dipole moment of an exciton on every site in the system. 
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `exciton_lifetimes`: Vector of exciton lifetimes (in s) in the donor and acceptor.
- `CT_lifetimes`: Vector of CT state lifetimes (in s) in the donor, at the interface, and in the acceptor.
- `bath_reorganisation_energies`: Vector of reorganisation energies of the bath for electrons, holes, and excitons (in meV).  
- `accuracy`: The accuracy of dKMC calculations (a_dKMC).
- `maximum_hops_cutoff`: Maximum number of hops the simulation will run before terminating.
- `separation_cutoff`: Distance (in m) at which charges are considered separated, and the simulation terminates. 
- `exciton_population_cutoff`: Minimum population on exciton site-pairs to treat the current state as an exciton.
- `maximum_excitation_distance`: Maximum distance (in m) that an exciton can be excited from the interface.
- `hopping_radii`: Matrix of precalculated distances that electrons, holes, and excitons can hop in the donor (row 1) and acceptor (row 2). 
- `hamiltonian_radii`: Matrix of precalculated distance for electron, hole, and exciton sites in the donor (row 1) and acceptor (row 2) that are included as site-pairs in Hamiltonian subset.
- `K_tots`: Vector of precaclulated K_tot values for electrons, holes, and excitons required for dKMC rate calculations.
- `E_steps`: Vector of the energy step used for values that K_tot is calculated at for electrons, holes, and excitons.
- `E_limits`: Vector of the energy limits used for values that K_tot is calculated at for electrons, holes, and excitons.        

# Output:
- `outcome`: String describing the trajectories outcome. 
- `t`: Elapsed time (in s).
- `initial_state_characteristics`: Vector containing the energy (in meV), IPR, and separation (in nm) of the initial state.
- `final_state_characteristics`: Vector containing the energy (in meV), IPR, and separation (in nm) of the final state.

"""
function dKMC_charge_generation(dimension::Integer,N::Integer,LUMO_HOMO_energies::Matrix{<:AbstractFloat},exciton_binding_energies::Vector{<:Number},electronic_couplings::Matrix{<:Number},kappas::Vector{<:AbstractFloat},transition_dipole_moments::Vector{<:Number},dipole_orientations::Matrix{<:AbstractFloat},epsilon_r::Number,site_spacing::Number,exciton_lifetimes::Vector{<:Number},CT_lifetimes::Vector{<:Number},bath_reorganisation_energies::Vector{<:Number},accuracy::AbstractFloat,maximum_hops_cutoff::Integer,separation_cutoff::Number,exciton_population_cutoff::AbstractFloat,maximum_excitation_distance::Number,hopping_radii::Matrix{<:Number},hamiltonian_radii::Matrix{<:Number},K_tots::Vector{Matrix{ComplexF64}},E_steps::Vector{<:Number},E_limits::Vector{<:Number})

    #Randomly selecting how far an exciton is excited from the interface, between 1 and maximum_excitation_distance sites.
    excitation_distance = Int(ceil(maximum_excitation_distance*rand()/site_spacing))

    #Defining the initial position of the electron, hole and exciton.
    initial_location = ones(dimension).* (N/2)
    initial_location[1] -= excitation_distance
    current_electron_location = current_hole_location = current_exciton_location = initial_location
    exciton = true 
    
    #Calculating the current Hamiltonian, the polaron transformed Hamiltonian, the reference matrices.
    H,Ht,electron_r,hole_r,dipoles,transformed_coupling,electron_index,hole_index,exciton_index,bath_index,site_pair_indexes = setup_hamiltonian.current_charge_generation_hamiltonian(dimension,N,exciton,LUMO_HOMO_energies,exciton_binding_energies,electronic_couplings,transition_dipole_moments,dipole_orientations,epsilon_r,site_spacing,bath_reorganisation_energies,kappas,hamiltonian_radii,current_electron_location,current_hole_location,current_exciton_location)

    #Finding which sites correspond to exciton and CT site pairs.
    donor_exciton_sites, acceptor_exciton_sites, donor_CT_sites, acceptor_CT_sites, interfacial_CT_sites = find_exciton_and_CT_site_pairs(N,electron_r,hole_r)
    exciton_sites = vcat(donor_exciton_sites,acceptor_exciton_sites)

    #Calculating eigenvectors and eigenvalues in the energy eigenbasis.
    evals,evecs = eigen(Ht)

    #Calculating list of exciton states, as those with populations on exciton sites greater than exciton_population_cutoff.
    excitons = findall(x->x>exciton_population_cutoff,sum(evecs[exciton_sites,:].^2,dims=1)[:])

    #Calculating the average positions of energy eigenstates.
    electron_centres = [[evecs[:,i]' * Diagonal(electron_r[:,j]) * evecs[:,i] for i=eachindex(evals)] for j in 1:dimension]
    hole_centres = [[evecs[:,i]' * Diagonal(hole_r[:,j]) * evecs[:,i] for i=eachindex(evals)] for j in 1:dimension]
    exciton_centres = [[(evecs[exciton_sites,i]' * Diagonal(electron_r[:,j][exciton_sites]) * evecs[exciton_sites,i])/(evecs[exciton_sites,i]'*evecs[exciton_sites,i]) for i=eachindex(evals)] for j in 1:dimension]
    
    #Finding donor exciton states within a radius of 1 of the chosen initial position of the electron.
    excitation_radius = 1
    accessible_excitons = excitons[findall(x->x<=excitation_radius,sqrt.(sum([(exciton_centres[i][excitons] .- current_exciton_location[i]).^2 for i in 1:dimension])))]
    while isempty(accessible_excitons)
        excitation_radius += 1
        accessible_excitons = excitons[findall(x->x<=excitation_radius,sqrt.(sum([(exciton_centres[i][excitons] .- current_exciton_location[i]).^2 for i in 1:dimension])))]
    end

    #Choosing the initial state from accessible_excitons probabilistically based on the expectation value of each state's dipole moment.
    current_state = choosing_initial_exciton_state(accessible_excitons,evecs,dipoles)

    #Determining whether the current state is an exciton.
    if current_state in excitons
        exciton = true
    else
        exciton = false
    end

    #Recording characteristics of the current state.
    current_electron_location = [electron_centres[i][current_state] for i in 1:dimension]
    current_hole_location = [hole_centres[i][current_state] for i in 1:dimension]
    current_exciton_location = [exciton_centres[i][current_state] for i in 1:dimension]
    current_energy = evals[current_state]
    current_IPR = 1/sum(evecs[:,current_state].^4)
    current_separation = sqrt(sum((current_electron_location .- current_hole_location).^2)) * site_spacing
    
    #Saving the characteristics of the initial state.
    initial_state_characteristics = [current_energy current_IPR constants.to_nano(current_separation)] 

    #Setting the initial conditions.
    t = 0
    hops = 0
    outcome = ""

    #Iterate the dKMC procedure until one of the terminating conditions.
    while hops < maximum_hops_cutoff
        
        #Checking if either charge is too close to the edge of the landscape to create a new subset of the Hamiltonian.
        if iszero(current_electron_location .- hamiltonian_radii[1,1] .< 1) == false || iszero(current_electron_location .+ hamiltonian_radii[1,2] .> N) == false || iszero(current_hole_location .- hamiltonian_radii[2,1] .< 1) == false || iszero(current_hole_location .+ hamiltonian_radii[2,2] .> N) == false
            @warn "Simulation terminated as a charge carrier got too close to the edge of the system. Increase the length of the system (N) to avoid this error."
            outcome = "reached the edge of the system"
            break
        end

        #If currently in an exciton state, checking if the exciton is too close to the edge of the landscape to create a new subset of the Hamiltonian.
        if exciton == true  
            if iszero(current_exciton_location .- hamiltonian_radii[3,1] .< 1) == false || iszero(current_exciton_location .+ hamiltonian_radii[3,2] .> N) == false
                @warn "Simulation terminated as a charge carrier got too close to the edge of the system. Increase the length of the system (N) to avoid this error."
                outcome = "reached the edge of the system"
                break
            end 
        end 

        #Rediagonalise a new subset of the Hamiltonian.
        previous_eigenstate = evecs[:,current_state]
        previous_site_pair_indexes = copy(site_pair_indexes)
        H,Ht,electron_r,hole_r,dipoles,transformed_coupling,electron_index,hole_index,exciton_index,bath_index,site_pair_indexes = setup_hamiltonian.current_charge_generation_hamiltonian(dimension,N,exciton,LUMO_HOMO_energies,exciton_binding_energies,electronic_couplings,transition_dipole_moments,dipole_orientations,epsilon_r,site_spacing,bath_reorganisation_energies,kappas,hamiltonian_radii,current_electron_location,current_hole_location,current_exciton_location)
        donor_exciton_sites, acceptor_exciton_sites, donor_CT_sites, acceptor_CT_sites, interfacial_CT_sites = find_exciton_and_CT_site_pairs(N,electron_r,hole_r)
        exciton_sites = vcat(donor_exciton_sites,acceptor_exciton_sites)
        evals,evecs = eigen(Ht)
        excitons = findall(x->x>exciton_population_cutoff,sum(evecs[exciton_sites,:].^2,dims=1)[:])
        electron_centres = [[evecs[:,i]' * Diagonal(electron_r[:,j]) * evecs[:,i] for i=eachindex(evals)] for j in 1:dimension]
        hole_centres = [[evecs[:,i]' * Diagonal(hole_r[:,j]) * evecs[:,i] for i=eachindex(evals)] for j in 1:dimension]
        current_state = argmax(((previous_eigenstate[findall(x->x in site_pair_indexes,previous_site_pair_indexes)]' * evecs[findall(x->x in previous_site_pair_indexes,site_pair_indexes),:])[:]).^2)
        current_electron_location = [electron_centres[i][current_state] for i in 1:dimension]
        current_hole_location = [hole_centres[i][current_state] for i in 1:dimension]

        #Assigning which approximating_radii to use for each charge.
        if current_electron_location[1] <= N/2
            electron_hopping_radius = hopping_radii[1,1]
        else
            electron_hopping_radius = hopping_radii[1,2]
        end
        if current_hole_location[1] <= N/2
            hole_hopping_radius = hopping_radii[2,1]
        else
            hole_hopping_radius = hopping_radii[2,2]
        end
        if current_exciton_location[1] <= N/2
            exciton_hopping_radius = hopping_radii[3,1]
        else
            exciton_hopping_radius = hopping_radii[3,2]
        end

        #Finding which states are accessible from the current state.
        if exciton == true
            exciton_centres = [[(evecs[exciton_sites,i]' * Diagonal(electron_r[exciton_sites,j]) * evecs[exciton_sites,i])/(evecs[exciton_sites,i]'*evecs[exciton_sites,i]) for i=eachindex(evals)] for j in 1:dimension]
            current_exciton_location = [exciton_centres[i][current_state] for i in 1:dimension]
            accessible_exciton_states = excitons[findall(x-> 0 < x <1,(sum([(exciton_centres[i][excitons] .- current_exciton_location[i]).^2 for i in 1:dimension]))./(exciton_hopping_radius^2))]
            accessible_charges_states = setdiff(findall(x-> 0 < x < 1,(sum([(electron_centres[i] .- current_electron_location[i]).^2 for i in 1:dimension]))./(electron_hopping_radius^2) .+ (sum([(hole_centres[i] .- current_hole_location[i]).^2 for i in 1:dimension]))./(hole_hopping_radius^2)),accessible_exciton_states)
            accessible_states = vcat(accessible_exciton_states,accessible_charges_states)
        elseif exciton == false
            accessible_states = findall(x-> 0 < x < 1,(sum([(electron_centres[i] .- current_electron_location[i]).^2 for i in 1:dimension]))./(electron_hopping_radius^2) .+ (sum([(hole_centres[i] .- current_hole_location[i]).^2 for i in 1:dimension]))./(hole_hopping_radius^2))
        end

        #Calculating hopping rates to all states in accessible_states.
        hopping_rates = zeros(length(accessible_states))
        current_state_relevant_sites = sortperm(abs.(evecs[:,current_state]),rev=true)[1:findfirst(x->x>accuracy,cumsum(sort(abs.(evecs[:,current_state])./sum(abs.(evecs[:,current_state])),rev=true)))]
        for f = eachindex(accessible_states)
            destination_state_relevant_sites = sortperm(abs.(evecs[:,accessible_states[f]]),rev=true)[1:findfirst(x->x>accuracy,cumsum(sort(abs.(evecs[:,accessible_states[f]])./sum(abs.(evecs[:,accessible_states[f]])),rev=true)))]
            hopping_rate = dKMC_hopping_rates.charge_generation_dKMC_rate(current_state,accessible_states[f],transformed_coupling,evals,evecs,K_tots,E_steps,E_limits,current_state_relevant_sites,destination_state_relevant_sites,electron_index,hole_index,exciton_index,bath_index)
            if hopping_rate > 0
                hopping_rates[f] = hopping_rate
            end 
        end 

        #Calculating the exciton recombination rates, first donor and then acceptor.
        push!(hopping_rates,(1/constants.to_femto(exciton_lifetimes[1]))*(sum(evecs[donor_exciton_sites,current_state])^2))
        push!(hopping_rates,(1/constants.to_femto(exciton_lifetimes[2]))*(sum(evecs[acceptor_exciton_sites,current_state])^2))

        #Calculating the CT recombination rates, first donor, and then interfacial, and then acceptor.
        push!(hopping_rates,(1/constants.to_femto(CT_lifetimes[1]))*(sum(evecs[donor_CT_sites,current_state])^2))
        push!(hopping_rates,(1/constants.to_femto(CT_lifetimes[2]))*(sum(evecs[interfacial_CT_sites,current_state])^2))
        push!(hopping_rates,(1/constants.to_femto(CT_lifetimes[3]))*(sum(evecs[acceptor_CT_sites,current_state])^2))

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

        #Checking if state has recombined, either through exciton or CT recombination pathways, if so we break.
        if index == length(hopping_rates)-4
            outcome = "exciton recombination"
            break
        elseif index == length(hopping_rates)-3
            outcome = "exciton recombination"
            break	
        elseif index == length(hopping_rates)-2
            outcome = "bulk CT recombination"
            break
        elseif index == length(hopping_rates)-1
            outcome = "interfacial CT recombination"
            break
        elseif index == length(hopping_rates)
            outcome = "bulk CT recombination"
            break
        end

        #Updating the current state.
        current_state = accessible_states[index]

        #Recording the characteristics of the new state.
        if current_state in excitons
            exciton = true
        else
            exciton = false
        end	
        current_electron_location = [electron_centres[i][current_state] for i in 1:dimension]
        current_hole_location = [hole_centres[i][current_state] for i in 1:dimension]
        current_energy = evals[current_state]
        current_IPR = 1/sum(evecs[:,current_state].^4)
        current_separation = sqrt(sum((current_electron_location .- current_hole_location).^2)) * site_spacing

        #Checking if charges are separated by greater than separation_cutoff, if so we break.
        if current_separation > separation_cutoff
            if current_electron_location[1] < N/2 && current_hole_location[1] < N/2
                outcome = "bulk separation"
            elseif current_electron_location[1] > N/2 + 1 && current_hole_location[1] > N/2 + 1
                outcome = "bulk separation"
            else
                outcome = "interfacial separation"
            end
            break
        end
    end

    #Checking whether we have had more hops than maximum_hops_cutoff.
    if hops >= maximum_hops_cutoff
        outcome = "exceeded maximum hops cutoff"
    end

    #Saving the characteristics of the final state.
    final_state_characteristics = [current_energy current_IPR constants.to_nano(current_separation)]

    return outcome, t, initial_state_characteristics, final_state_characteristics

end


"""
    choosing_initial_exciton_state(excitons::Vector{<:Integer}},evecs::Matrix{<:AbstractFloat},
    dipoles::Matrix{<:AbstractFloat})

Chooses which of the provided exciton states is excited. The excited exciton state is chosen with a 
probability proportional to the square of the exectation value of the states transition dipole moment.

# Arguments:
- `excitons`: List of indexes for states defined as excitons
- `evecs`: Energy eigenvectors of Hamiltonian.
- `dipoles`: List of dipole orientations of all sites in the Hamiltonian.

# Output:
- `initial_exciton_state`: Index for the initially excited exciton.

"""
function choosing_initial_exciton_state(excitons::Vector{<:Integer},evecs::Matrix{<:AbstractFloat},dipoles::Matrix{<:AbstractFloat})
    
    #Calculating the square of the expectation value of the dipole moment of the energy_accessible_excitons
    dipole_expectations = [(evecs[:,i]' * Diagonal(dipoles[:,1]) * evecs[:,i]).^2 .+ (evecs[:,i]' * Diagonal(dipoles[:,2]) * evecs[:,i]).^2 .+ (evecs[:,i]' * Diagonal(dipoles[:,3]) * evecs[:,i]).^2 for i in excitons]

    #Total sum of these dipole expectations squared
    dipole_sum = sum(dipole_expectations)

    #Choose a random number and choose the initial state probabalistically based on this
    random_number = rand()
    index = findfirst(x->x>random_number*dipole_sum,cumsum(dipole_expectations))
    initial_exciton_state = excitons[index]

    return initial_exciton_state

end


"""
    find_exciton_and_CT_site_pairs(N::Integer,electron_r::Matrix{<:Integer},hole_r::Matrix{<:Integer})

Uses the lists of positions of electrons and holes on every site-pair to find which site-pairs correspond to different 
types of exciton and CT site-pairs, including donor and acceptor exciton site-pairs, and donor, acceptor and 
interfacial CT site-pairs.

# Arguments:
- `N`: Length of system, i.e., number of sites in each direction.
- `electron_r`: Electron positions of every site; r_e_x in 1D, [r_e_x r_e_y] in 2D, [r_e_x r_e_y r_e_z] in 3D.
- `hole_r`: Hole positions of every site; r_h_x in 1D, [r_h_x r_h_y] in 2D, [r_h_x r_h_y r_h_z] in 3D.

# Output:
- `donor_exciton_sites:	Vector of indices to site-pairs where the electron and hole are on the same donor site.
- `acceptor_exciton_sites`: Vector of indices to site-pairs where the electron and hole are on the same acceptor site.
- `donor_CT_sites`: Vector of indices to site-pairs where the electron and hole are on adjacent donor sites.
- `acceptor_CT_sites`: Vector of indices to site-pairs where the electron and hole are on adjacent acceptor sites.
- `interfacial_CT_sites`: Vector of indices to site-pairs where the electron and hole are on adjacent sites across the interface between the donor and acceptor.

"""
function find_exciton_and_CT_site_pairs(N::Integer,electron_r::Matrix{<:Integer},hole_r::Matrix{<:Integer})

    #Calculate the electron-hole separation for every site-pair.
    site_separations = sum([(electron_r[:,i] .- hole_r[:,i]).^2 for i in 1:size(electron_r)[2]]) 

    #Find all exciton site-pairs, those with an electron-hole separation of zero.
    exciton_sites = findall(x->x == 0,site_separations)

    #Split exciton site-pairs into those in the donor and those in the acceptor.
    donor_exciton_sites = intersect(findall(x->x<=N/2,electron_r[:,1]),exciton_sites)
    acceptor_exciton_sites = intersect(findall(x->x>N/2,electron_r[:,1]),exciton_sites)

    #Find all adjacent site-pairs, those with a electron-hole separation of one.
    adjacent_pair_sites = findall(x->x == 1,site_separations)

    #Split adjacent site pairs into those in the donor, those in the acceptor and those across the interface.
    electron_in_donor_interfacial_CT_sites = intersect(findall(x->x==Int(N/2),electron_r[:,1]),findall(x->x==Int(N/2+1),hole_r[:,1]),adjacent_pair_sites)
    electron_in_acceptor_interfacial_CT_sites = intersect(findall(x->x==Int(N/2+1),electron_r[:,1]),findall(x->x==Int(N/2),hole_r[:,1]),adjacent_pair_sites)
    interfacial_CT_sites = vcat(electron_in_donor_interfacial_CT_sites,electron_in_acceptor_interfacial_CT_sites)
    donor_adjacent_pair_sites = intersect(findall(x->x<=N/2,electron_r[:,1]),adjacent_pair_sites)
    acceptor_adjacent_pair_sites = intersect(findall(x->x>N/2,electron_r[:,1]),adjacent_pair_sites)
    donor_CT_sites = setdiff(donor_adjacent_pair_sites,electron_in_donor_interfacial_CT_sites)
    acceptor_CT_sites = setdiff(acceptor_adjacent_pair_sites,electron_in_acceptor_interfacial_CT_sites)

    return donor_exciton_sites, acceptor_exciton_sites, donor_CT_sites, acceptor_CT_sites, interfacial_CT_sites

end

end