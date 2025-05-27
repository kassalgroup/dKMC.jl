# dKMC Exciton Transport:
This dKMC module simulates the movement of a partially delocalised exciton in a disordered material. These simulations are similar to those for charge transport, except that transport is mediated by long-range dipole-dipole excitonic couplings instead of nearest-neighbour electronic couplings. The simulations start with an exciton in the delocalised state closest to the middle of the lattice and propagate the transport dynamics using dKMC. The simulation takes in parameters describing the material, and returns the exciton's diffusion coefficient, as well as its mean-squared displacment, energy, and delocalisation (IPR) as a function of time.

For details, see:
[Balzer, D.; Kassal, I. Mechanism of Delocalization-Enhanced Exciton Transport in Disordered Organic Semiconductors. *Journal of Physical Chemisry Letters* **2023**, *14*(8), 2155-2162.](https://pubs.acs.org/doi/10.1021/acs.jpclett.2c03886)

## Input File:
The input parameters are to be supplied in YAML format; for an example, see `exciton_transport_input.in`.

Parameters required in the input file:
- `number_of_processes`: The number of computer processes available for parallelisation.
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `exciton_disorder`: Excitonic disorder (in meV).
- `transition_dipole_moment`: Magnitude of the transition dipole moment on every site (in D).
- `exciton_bath_reorganisation_energy`: Reorganisation energy of the exciton's bath (in meV). This is used in the bath spectral density, which is assumed to be super-Ohmic.
- `exciton_bath_cutoff_energy`: Cutoff energy of the exciton's bath (in meV). This is used in the bath spectral density, which is assumed to be super-Ohmic.
- `T`: Temperature (in K).
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `landscape_iterations`: Number of realisations of disorder (or energetic landscapes) for simulations to be run on.
- `trajectory_iterations`: Number of trajectories to be simulated on each realisation of disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC). A number between 0 and 1, typically chosen to be 0.99. Increasing the accuracy will reduce the error introduced by dKMC's approximations but will increase the computational cost.
- `end_time`: The simulation end time (in s).
- `number_of_sampling_times`: The number of times between 0 and end_time to sample the squared displacement, energies and IPRs.

## How to use:
1. Fill in the input file with the desired simulation parameters and save the file.
2. In terminal, run 'julia run_dKMC_exciton_transport.jl exciton_transport_input.in exciton_transport_out.out'.
A sample output file *exciton_transport_out.out* is also provided.