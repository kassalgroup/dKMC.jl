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

#Record function introduction and input file to output file.
output = open(output_file, "w");                                                                                                                    
print(output,"""
────────────────────────────────────────────────────────────────────────────────
dKMC.jl v1.0 
Charge Transport Module
Github: https://github.com/kassalgroup/dKMC.jl

Daniel Balzer
University of Sydney
Kassal Group: https://www.kassal.group

This dKMC module simulates the movement of a partially delocalised charge 
carrier, either an electron or a hole, in a disordered material. The simulations 
start with a charge in the middle of the lattice and propagate the transport 
dynamics using dKMC.

If this module is used, please cite the following:

1. Balzer, D.; Smolders, T. J. A. M.; Blyth, D.; Hood, S. N.; Kassal, I. 
Delocalised Kinetic Monte Carlo for Simulating Delocalisation-Enhanced Charge 
and Exciton Transport in Disordered Materials. Chemical Science 2021, 12(6), 
2276–2285. 

────────────────────────────────────────────────────────────────────────────────

Julia version:  $VERSION
Start time:     $(DateTime(now()))

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
@everywhere include("dKMC_charge_transport_functions.jl")

#Run the dKMC charge transport calculations.
mobility,sampling_time_squared_displacements,sampling_time_energies,sampling_time_IPRs,total_hops,reached_boundary_proportion = dKMC_charge_transport_functions.dKMC_charge_transport_results(inputs["dimension"],inputs["N"],inputs["disorder"],inputs["electronic_coupling"],inputs["bath_reorganisation_energy"],inputs["bath_cutoff_energy"],inputs["T"],inputs["site_spacing"],inputs["landscape_iterations"],inputs["trajectory_iterations"],inputs["accuracy"],inputs["end_time"],inputs["number_of_sampling_times"])

#Record results and timer output to the output file.
print(output,"""
────────────────────────────────────────────────────────────────────────────────
Results:
────────────────────────────────────────────────────────────────────────────────
Mobility: $(mobility[1]) cm²V⁻¹s⁻¹
Mobility error: $(mobility[2]) cm²V⁻¹s⁻¹

Sampling times (s): $(collect(0:inputs["end_time"]/inputs["number_of_sampling_times"]:inputs["end_time"]))

Mean squared displacements (in nm²) at sampling times: $(sampling_time_squared_displacements[1,:])
Standard errors of the squared displacements (in nm²) at sampling times: $(sampling_time_squared_displacements[2,:])

Mean energies (in meV) at sampling times: $(sampling_time_energies[1,:])
Standard errors of the energies (in meV) at sampling times: $(sampling_time_energies[2,:])

Mean IPRs at sampling times: $(sampling_time_IPRs[1,:])
Standard errors of the IPRs at sampling times: $(sampling_time_IPRs[2,:])

Mean number of hops: $(total_hops[1]) ± $(total_hops[2])
""")

if reached_boundary_proportion > 0.01
    print(output,"""
    ────────────────────────────────────────────────────────────────────────────────

    WARNING: "A large proportion ($(reached_boundary_proportion)) of trajectories terminated early as a charge carrier got too close to the edge of the system. Consider increasing N and repeating the calculation."

    """)
end

print(output,"""
────────────────────────────────────────────────────────────────────────────────

Evaluation completed successfully.
End time:       $(DateTime(now()))

────────────────────────────────────────────────────────────────────────────────
""")
close(output)