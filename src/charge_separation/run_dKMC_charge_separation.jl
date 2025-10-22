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
disorders = [inputs["acceptor_electron_disorder"], inputs["donor_hole_disorder"]]
electronic_couplings = [inputs["acceptor_electron_coupling"], inputs["donor_hole_coupling"]]
bath_reorganisation_energies = [inputs["electron_bath_reorganisation_energy"], inputs["hole_bath_reorganisation_energy"]]
bath_cutoff_energies = [inputs["electron_bath_cutoff_energy"], inputs["hole_bath_cutoff_energy"]] 

#Record the starting time.
start_time = DateTime(now())

#Record function introduction and input file to output file.
output = open(output_file, "w");
print(output,"""
────────────────────────────────────────────────────────────────────────────────
dKMC.jl v1.0
Charge Separation Module
Github: https://github.com/kassalgroup/dKMC.jl

Daniel Balzer
University of Sydney
Kassal Group: https://www.kassal.group

This dKMC module simulates the separation of a partially delocalised 
electron-hole pair from an interfacial charge-transfer (CT) state to free 
charges. The simulations start from an interfacial CT state, i.e., an electron 
and a hole on neighbouring molecules on the opposite sides of an interface 
between the electron-donor and electron-acceptor materials. The simulation then 
uses dKMC to propagate the separation dynamics, with the electron constrained 
to the acceptor and the hole to the donor, until the charges either recombine 
or separate.

If this module is used, please cite the following:

1. Balzer, D.; Kassal, I. Even a Little Delocalisation Produces Large Kinetic 
Enhancements of Charge-Separation Efficiency in Organic Photovoltaics. Science 
Advances 2022, 8, eabl9692.

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
@everywhere include("dKMC_charge_separation_functions.jl")

#Run the dKMC charge separation calculations.
mean_outcomes,mean_initial_state_characteristics,mean_separation_time,mean_separated_initial_state_characteristics,mean_separated_final_state_characteristics = dKMC_charge_separation_functions.dKMC_charge_separation_results(inputs["dimension"],inputs["N"],disorders,inputs["donor_HOMO_acceptor_LUMO_gap"],electronic_couplings,inputs["epsilon_r"],inputs["site_spacing"],inputs["CT_lifetime"],bath_reorganisation_energies,bath_cutoff_energies,inputs["T"],inputs["landscape_iterations"],inputs["trajectory_iterations"],inputs["accuracy"],inputs["maximum_hops_cutoff"],inputs["separation_cutoff"])

#Record results and timer output to the output file.
print(output,"""
────────────────────────────────────────────────────────────────────────────────
Results:
────────────────────────────────────────────────────────────────────────────────
All uncertainties included below are standard errors of the mean.

Internal quantum efficiency: $(mean_outcomes[1,1]) ± $(mean_outcomes[2,1])

Proportion of all possible outcomes:
    1. Separation: $(mean_outcomes[1,1]) ± $(mean_outcomes[2,1])
    2. Recombination: $(mean_outcomes[1,2]) ± $(mean_outcomes[2,2])
    3. Exceeded maximum hops cutoff: $(mean_outcomes[1,3]) ± $(mean_outcomes[2,3])
    4. Reached the edge of the system: $(mean_outcomes[1,4]) ± $(mean_outcomes[2,4])

Mean characteristics of the initial state:
    Energy (meV): $(mean_initial_state_characteristics[1,1]) ± $(mean_initial_state_characteristics[2,1])
    Inverse participation ratio: $(mean_initial_state_characteristics[1,2]) ± $(mean_initial_state_characteristics[2,2])
    Electron-hole separation (nm): $(mean_initial_state_characteristics[1,3]) ± $(mean_initial_state_characteristics[2,3])

""")

if mean_separation_time != 0
    print(output,"""
    For trajectories where charges separate:
        Mean separation time (s): $(mean_separation_time[1]) ± $(mean_separation_time[2])
        Mean characteristics of the initial state:
            Energy (meV): $(mean_separated_initial_state_characteristics[1,1]) ± $(mean_separated_initial_state_characteristics[2,1])
            Inverse participation ratio: $(mean_separated_initial_state_characteristics[1,2]) ± $(mean_separated_initial_state_characteristics[2,2])
            Electron-hole separation: $(mean_separated_initial_state_characteristics[1,3]) ± $(mean_separated_initial_state_characteristics[2,3])
        Mean characteristics of the final state:
            Energy (meV): $(mean_separated_final_state_characteristics[1,1]) ± $(mean_separated_final_state_characteristics[2,1])
            Inverse participation ratio: $(mean_separated_final_state_characteristics[1,2]) ± $(mean_separated_final_state_characteristics[2,2])
            Electron-hole separation (nm): $(mean_separated_final_state_characteristics[1,3]) ± $(mean_separated_final_state_characteristics[2,3])

    """)
end

if mean_outcomes[1,3] > 0.05
    print(output,"""
    ────────────────────────────────────────────────────────────────────────────────

    WARNING: "A large proportion ($(mean_outcomes[1,3])) of trajectories exceeded the maximum hops cutoff. Consider increasing maximum_hops_cutoff and repeating the calculation."

    """)
end

if mean_outcomes[1,4] > 0.01
    print(output,"""
    ────────────────────────────────────────────────────────────────────────────────

    WARNING: "A large proportion ($(mean_outcomes[1,4])) of trajectories terminated early as a charge carrier got too close to the edge of the system. Consider increasing N and repeating the calculation."

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