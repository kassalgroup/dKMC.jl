#Print a message to indicate the calculation is running.
println("Welcome to dKMC, the calculation is now running.")

#Activate the dKMC package.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."), io=devnull)

#Loading the required packages.
include("../shared_functions/package_loading.jl")

#Change directory to location of this file.
cd(@__DIR__)

#Define the input and output files.
if length(ARGS) == 2
    input_file, output_file = ARGS
else
    error("Input and output files were not supplied.")
end

#Read input files.
raw_input = open(input_file) do file
    read(file, String)
end
inputs = YAML.load(raw_input)

#Group together like parameters into vector inputs.
disorders = [inputs["donor_electron_disorder"] inputs["acceptor_electron_disorder"]; inputs["donor_hole_disorder"] inputs["acceptor_hole_disorder"]]
exciton_disorders = [inputs["donor_exciton_disorder"], inputs["acceptor_exciton_disorder"]]
exciton_binding_energies = [inputs["donor_exciton_binding_energy"], inputs["acceptor_exciton_binding_energy"]]
electronic_couplings = [inputs["donor_electron_coupling"] inputs["interface_electron_coupling"] inputs["acceptor_electron_coupling"]; inputs["donor_hole_coupling"] inputs["interface_hole_coupling"] inputs["acceptor_hole_coupling"]]
transition_dipole_moments = [inputs["donor_transition_dipole_moment"], inputs["acceptor_transition_dipole_moment"]]
exciton_lifetimes = [inputs["donor_exciton_lifetime"], inputs["acceptor_exciton_lifetime"]]
CT_lifetimes = [inputs["donor_CT_lifetime"], inputs["interfacial_CT_lifetime"], inputs["acceptor_CT_lifetime"]]
bath_reorganisation_energies = [inputs["electron_bath_reorganisation_energy"], inputs["hole_bath_reorganisation_energy"], inputs["exciton_bath_reorganisation_energy"]]
bath_cutoff_energies = [inputs["electron_bath_cutoff_energy"], inputs["hole_bath_cutoff_energy"], inputs["exciton_bath_cutoff_energy"]] 

#Record the starting time.
start_time = DateTime(now())

#Record function introduction and input file to output file.
output = open(output_file, "w");
print(output,"""
────────────────────────────────────────────────────────────────────────────────
dKMC.jl v1.0
Charge Generation Module
Github: https://github.com/kassalgroup/dKMC.jl

Daniel Balzer
University of Sydney
Kassal Group: https://www.kassal.group

This dKMC module simulates the complete process of the generation of partially 
delocalised free charges from excitons. The simulations start from an exciton 
in the donor material and near the interface with the acceptor material. The 
simulation uses dKMC to propagate the full two-paricle dynamics, where the 
electron and hole can occupy any donor or acceptor site, until the charges 
either recombine or become separated.

If this module is used, please cite the following:

1. Balzer, D.; Kassal, I. Delocalisation enables efficient charge generation in 
organic photovoltaics, even with little to no energetic offset. Chemical Science 
2024, 15, 4779.

────────────────────────────────────────────────────────────────────────────────

Julia version:  $VERSION
Start time:     $start_time

────────────────────────────────────────────────────────────────────────────────
User input:
────────────────────────────────────────────────────────────────────────────────
$raw_input

""")

#If there are multiple processes available, add these proccesses for parallel computing.
if inputs["number_of_processes"] > 1
    addprocs(inputs["number_of_processes"])
end

#Make sure the appropriate functions are available to all processes.
@everywhere include("dKMC_charge_generation_functions.jl")

#Run the dKMC charge generation calculations.
mean_outcomes, mean_initial_state_characteristics,mean_separation_time, mean_separated_initial_state_characteristics, mean_separated_final_state_characteristics, mean_interfacial_separation_time, mean_interfacial_separated_initial_state_characteristics, mean_interfacial_separated_final_state_characteristics, mean_bulk_separation_time, mean_bulk_separated_initial_state_characteristics, mean_bulk_separated_final_state_characteristics = dKMC_charge_generation_functions.dKMC_charge_generation_results(inputs["dimension"],inputs["N"],disorders,exciton_disorders,inputs["donor_HOMO_LUMO_gap"],inputs["LUMO_offset"],inputs["HOMO_offset"],exciton_binding_energies,electronic_couplings,transition_dipole_moments,inputs["epsilon_r"],inputs["site_spacing"],exciton_lifetimes,CT_lifetimes,bath_reorganisation_energies,bath_cutoff_energies,inputs["T"],inputs["landscape_iterations"],inputs["trajectory_iterations"],inputs["accuracy"],inputs["maximum_hops_cutoff"],inputs["separation_cutoff"],inputs["exciton_population_cutoff"],inputs["maximum_excitation_distance"])

#Record results and timer output to the output file.
print(output,"""
────────────────────────────────────────────────────────────────────────────────
Results:
────────────────────────────────────────────────────────────────────────────────
All uncertainties included below are standard errors of the mean.

Internal quantum efficiency (IQE): $(mean_outcomes[1,1]+mean_outcomes[1,2]) ± $(sqrt(mean_outcomes[2,1]^2 + mean_outcomes[2,2]^2))

Proportion of all possible outcomes:
    1. Separation: $(mean_outcomes[1,1]+mean_outcomes[1,2]) ± $(sqrt(mean_outcomes[2,1]^2 + mean_outcomes[2,2]^2))
        a. Interfacial separation: $(mean_outcomes[1,1]) ± $(mean_outcomes[2,1])
        b. Bulk separation: $(mean_outcomes[1,2]) ± $(mean_outcomes[2,2])
    2. Recombination: $(mean_outcomes[1,3]+mean_outcomes[1,4]+mean_outcomes[1,5]) ± $(sqrt(mean_outcomes[2,3]^2 + mean_outcomes[2,4]^2 + mean_outcomes[2,5]^2))
        a. Exciton recombination: $(mean_outcomes[1,3]) ± $(mean_outcomes[2,3])
        b. CT recombination: $(mean_outcomes[1,4]+mean_outcomes[1,5]) ± $(sqrt(mean_outcomes[2,4]^2+mean_outcomes[2,5]^2))
            i. Interfacial CT recombination: $(mean_outcomes[1,4]) ± $(mean_outcomes[2,6])
            ii. Bulk CT recombination: $(mean_outcomes[1,5]) ± $(mean_outcomes[2,5])
    3. Exceeded maximum hops cutoff: $(mean_outcomes[1,6]) ± $(mean_outcomes[2,6])
    4. Reached the edge of the system: $(mean_outcomes[1,7]) ± $(mean_outcomes[2,7])

Mean characteristics of the initial state:
    Energy (meV): $(mean_initial_state_characteristics[1,1]) ± $(mean_initial_state_characteristics[2,1])
    Inverse participation ratio: $(mean_initial_state_characteristics[1,2]) ± $(mean_initial_state_characteristics[2,2])
    Electron-hole separation (nm): $(mean_initial_state_characteristics[1,3]) ± $(mean_initial_state_characteristics[2,3])

""")

if mean_separation_time != 0
    print(output,"""
    For all trajectories where charges separate:
        Mean separation time (s): $(mean_separation_time[1]) ± $(mean_separation_time[2])
        Mean characteristics of the initial state:
            Energy (meV): $(mean_separated_initial_state_characteristics[1,1]) ± $(mean_separated_initial_state_characteristics[2,1])
            Inverse participation ratio: $(mean_separated_initial_state_characteristics[1,2]) ± $(mean_separated_initial_state_characteristics[2,2])
            Electron-hole separation (nm): $(mean_separated_initial_state_characteristics[1,3]) ± $(mean_separated_initial_state_characteristics[2,3])
        Mean characteristics of the final state:
            Energy (meV): $(mean_separated_final_state_characteristics[1,1]) ± $(mean_separated_final_state_characteristics[2,1])
            Inverse participation ratio: $(mean_separated_final_state_characteristics[1,2]) ± $(mean_separated_final_state_characteristics[2,2])
            Electron-hole separation (nm): $(mean_separated_final_state_characteristics[1,3]) ± $(mean_separated_final_state_characteristics[2,3])
    
    """)
end

if mean_interfacial_separation_time != 0.0
    print(output,"""
    For trajectories undergoing interfacial separation:
        Mean interfacial separation time (s): $(mean_interfacial_separation_time[1]) ± $(mean_interfacial_separation_time[2])
        Mean characteristics of the initial state:
            Energy (meV): $(mean_interfacial_separated_initial_state_characteristics[1,1]) ± $(mean_interfacial_separated_initial_state_characteristics[2,1])
            Inverse participation ratio: $(mean_interfacial_separated_initial_state_characteristics[1,2]) ± $(mean_interfacial_separated_initial_state_characteristics[2,2])
            Electron-hole separation (nm): $(mean_interfacial_separated_initial_state_characteristics[1,3]) ± $(mean_interfacial_separated_initial_state_characteristics[2,3])
        Mean characteristics of the final state:
            Energy (meV): $(mean_interfacial_separated_final_state_characteristics[1,1]) ± $(mean_interfacial_separated_final_state_characteristics[2,1])
            Inverse participation ratio: $(mean_interfacial_separated_final_state_characteristics[1,2]) ± $(mean_interfacial_separated_final_state_characteristics[2,2])
            Electron-hole separation (nm): $(mean_interfacial_separated_final_state_characteristics[1,3]) ± $(mean_interfacial_separated_final_state_characteristics[2,3])
    
    """)
end

if mean_bulk_separation_time != 0
    print(output,"""
    For trajectories undergoing bulk separation:
        Mean bulk separation time (s): $(mean_bulk_separation_time[1]) ± $(mean_bulk_separation_time[2])
        Mean characteristics of the initial state:
            Energy (meV): $(mean_bulk_separated_initial_state_characteristics[1,1]) ± $(mean_bulk_separated_initial_state_characteristics[2,1])
            Inverse participation ratio: $(mean_bulk_separated_initial_state_characteristics[1,2]) ± $(mean_bulk_separated_initial_state_characteristics[2,2])
            Electron-hole separation (nm): $(mean_bulk_separated_initial_state_characteristics[1,3]) ± $(mean_bulk_separated_initial_state_characteristics[2,3])
        Mean characteristics of the final state:
            Energy (meV): $(mean_bulk_separated_final_state_characteristics[1,1]) ± $(mean_bulk_separated_final_state_characteristics[2,1])
            Inverse participation ratio: $(mean_bulk_separated_final_state_characteristics[1,2]) ± $(mean_bulk_separated_final_state_characteristics[2,2])
            Electron-hole separation (nm): $(mean_bulk_separated_final_state_characteristics[1,3]) ± $(mean_bulk_separated_final_state_characteristics[2,3])
    
    """)
end

if mean_outcomes[1,6] > 0.05
    print(output,"""
    ────────────────────────────────────────────────────────────────────────────────

    WARNING: "A large proportion ($(mean_outcomes[1,6])) of trajectories exceeded the maximum hops cutoff. Consider increasing maximum_hops_cutoff and repeating the calculation."

    """)
end


if mean_outcomes[1,7] > 0.01
    print(output,"""
    ────────────────────────────────────────────────────────────────────────────────

    WARNING: "A large proportion ($(mean_outcomes[1,7])) of trajectories terminated early as a charge carrier got too close to the edge of the system. Consider increasing N and repeating the calculation."

    """)
end

#Record the end time and run time and print to output.
end_time = DateTime(now())
run_time = canonicalize(end_time - start_time)
print(output,"""
────────────────────────────────────────────────────────────────────────────────

Evaluation completed successfully.
End time:   $end_time
Run time:   $run_time 

────────────────────────────────────────────────────────────────────────────────
""")
close(output)

#Print a message to indicate the calculation is finished.
println("The dKMC calculation is finished.")