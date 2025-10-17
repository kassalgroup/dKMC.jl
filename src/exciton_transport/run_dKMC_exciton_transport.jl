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
Exciton Transport Module
Github: https://github.com/kassalgroup/dKMC.jl

Daniel Balzer
University of Sydney
Kassal Group: https://www.kassal.group

This dKMC module simulates the movement of a partially delocalised exciton in a 
disordered material. The simulations start with an exciton in the middle of the 
lattice and propagate the transport dynamics using dKMC. 

If this module is used, please cite the following:

1. Balzer, D.; Kassal, I. Mechanism of Delocalization-Enhanced Exciton Transport 
in Disordered Organic Semiconductors. Journal of Physical Chemisry Letters 2023, 
14(8), 2155-2162.

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
@everywhere include("dKMC_exciton_transport_functions.jl")

#Run the dKMC exciton transport calculations.
diffusion_coefficient,sampling_time_displacements,sampling_time_energies,sampling_time_IPRs,total_hops,reached_boundary_proportion = dKMC_exciton_transport_functions.dKMC_exciton_transport_results(inputs["dimension"],inputs["N"],inputs["exciton_disorder"],inputs["transition_dipole_moment"],inputs["epsilon_r"],inputs["exciton_bath_reorganisation_energy"],inputs["exciton_bath_cutoff_energy"],inputs["T"],inputs["site_spacing"],inputs["landscape_iterations"],inputs["trajectory_iterations"],inputs["accuracy"],inputs["end_time"],inputs["number_of_sampling_times"])

#Record results and timer output to the output file.
print(output,"""
────────────────────────────────────────────────────────────────────────────────
Results:
────────────────────────────────────────────────────────────────────────────────
Diffusion coefficient: $(diffusion_coefficient[1]) cm²s⁻¹
Diffusion coefficient error: $(diffusion_coefficient[2]) cm²s⁻¹

Sampling times (s): $(collect(0:inputs["end_time"]/inputs["number_of_sampling_times"]:inputs["end_time"]))

Mean squared displacements (in nm²) at sampling times: $(sampling_time_displacements[1,:])
Standard errors of the mean squared displacements (in nm²) at sampling times: $(sampling_time_displacements[2,:])

Mean energies (in meV) at sampling times: $(sampling_time_energies[1,:])
Standard errors of the mean energies (in meV) at sampling times: $(sampling_time_energies[2,:])

Mean IPRs at sampling times: $(sampling_time_IPRs[1,:])
Standard errors of the mean IPRs at sampling times: $(sampling_time_IPRs[2,:])

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