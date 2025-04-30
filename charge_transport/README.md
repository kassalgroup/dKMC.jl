# dKMC Charge Transport:
This dKMC module simulates the movement of a partially delocalised charge carrier, either an electron or a hole, in a disordered material. The simulations start with a charge in the delocalised state closest to the middle of the lattice and propagate the transport dynamics using dKMC. The simulation takes in parameters describing the material, and returns the charge carrier's mobility, as well as its mean-squared displacment, energy, and delocalisation (IPR) as a function of time.

For details, see:
[Balzer, D.; Smolders, T. J. A. M.; Blyth, D.; Hood, S. N.; Kassal, I. Delocalised Kinetic Monte Carlo for Simulating Delocalisation-Enhanced Charge and Exciton Transport in Disordered Materials. *Chemical Science* **2021**, *12*(6), 2276–2285.](https://pubs.rsc.org/en/content/articlelanding/2021/sc/d0sc04116e)

## Input File:
The input parameters are to be supplied in YAML format; for an example, see `charge_transport_input.in`.

Parameters required in the input file:
- `number_of_processes`: The number of computer processes available for parallelisation.
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorder`: Energetic disorder (in meV).
- `electronic_coupling`: Nearest neighbour electronic coupling (in meV).
- `bath_reorganisation_energy`: Reorganisation energy of the bath (in meV). This is used in the bath spectral density, which is assumed to be super-Ohmic.
- `bath_cutoff_energy`: Cutoff energy of the bath (in meV). This is used in the bath spectral density, which is assumed to be super-Ohmic.
- `T`: Temperature (in K).
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `landscape_iterations`: Number of realisations of disorder (or energetic landscapes) for simulations to be run on.
- `trajectory_iterations`: Number of trajectories to be simulated on each realisation of disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC). A number between 0 and 1, typically chosen to be 0.99. Increasing the accuracy will reduce the error introduced by dKMC's approximations but will increase the computational cost.
- `end_time`: The simulation end time (in s).
- `number_of_sampling_times`: The number of times between 0 and end_time to sample the squared displacement, energies and IPRs.

## How to use:
1. Fill in the input file with the desired simulation parameters and save the file.
2. In terminal, run 'julia run_dKMC_charge_transport.jl charge_transport_input.in charge_transport_out.out'.
A sample output file *charge_transport_out.out* is also provided.