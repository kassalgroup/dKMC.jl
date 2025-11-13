# dKMC charge generation:
This dKMC module simulates the complete process of the generation of partially delocalised free charges from excitons. The simulations start from a delocalised exciton in the donor material and near the interface with the acceptor material. The simulation uses dKMC to propagate the full two-particle dynamics, where the electron and hole can occupy any donor or acceptor site, until the charges either recombine or become separated. The simulation takes in parameters describing both materials and returns the internal quantum efficiency (IQE) of charge generation, as well as the mean characteristics (energy, separation, and IPR) of the initial and final states.

For details, see:
[Balzer, D.; Kassal, I. Delocalisation enables efficient charge generation in organic photovoltaics, even with little to no energetic offset. *Chemical Science* **2024**, *15*, 4779.](https://pubs.rsc.org/en/content/articlelanding/2024/sc/d3sc05409h)

The algorithm is the same as is published above, except that following rediagonalisation the new state is chosen as the state with the biggest overlap with the previous state rather than as the state closest in position.

## Input File:
The input parameters are to be supplied in YAML format; for an example, see `charge_transport_input.in`.

Parameters required in the input file:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `donor_electron_disorder`: Energetic disorder (in meV) of the electrons in the LUMO of the donor.
- `acceptor_electron_disorder`: Energetic disorder (in meV) of the electrons in the LUMO of the acceptor.
- `donor_hole_disorder`: Energetic disorder (in meV) of the holes in the HOMO of the donor.
- `acceptor_hole_disorder`: Energetic disorder (in meV) of the holes in the HOMO of the acceptor.
- `donor_exciton_disorder`: Excitonic disorder (in meV) in the donor.
- `acceptor_exciton_disorder`: Excitonic disorder (in meV) in the acceptor.
- `donor_HOMO_LUMO_GAP`: Energy gap between donor HOMO and LUMO levels (in meV).
- `LUMO_offset`: Energy offset between donor and acceptor LUMO levels (in meV).
- `HOMO_offset`: Energy offset between donor and acceptor HOMO levels (in meV).
- `donor_exciton_binding_energy`: Exciton binding energies for the donor (in meV).
- `acceptor_exciton_binding_energy`: Exciton binding energies for the acceptor (in meV).
- `donor_electron_coupling`: Nearest neighbour electronic couplings (in meV) for electron transfer in the donor.
- `interface_electron_coupling`: Nearest neighbour electronic couplings (in meV) for electron transfer across the interface.
- `acceptor_electron_coupling`: Nearest neighbour electronic couplings (in meV) for electron transfer in the acceptor.
- `donor_hole_coupling`: Nearest neighbour electronic couplings (in meV) for hole transfer in the donor.
- `interface_hole_coupling`: Nearest neighbour electronic couplings (in meV) for hole transfer across the interface.
- `acceptor_hole_coupling`: Nearest neighbour electronic couplings (in meV) for hole transfer in the acceptor.
- `donor_transition_dipole_moment`: Magnitude of the transition dipole moments on donor sites (in D).
- `acceptor_transition_dipole_moment`: Magnitude of the transition dipole moments on acceptor sites (in D).
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `donor_exciton_lifetime`: Exciton lifetime in the donor (in s).
- `acceptor_exciton_lifetime`: Exciton lifetime in the acceptor (in s).
- `donor_CT_lifetime`: CT state lifetime in the donor (in s).
- `interfacial_CT_lifetime`: CT state lifetime at the interface(in s).
- `acceptor_CT_lifetime`: CT state lifeime in the acceptor (in s).
- `electron_bath_reorganisation_energy`: Reorganisation energies of the bath for electrons (in meV). This is used in the electron's bath spectral density, which is assumed to be super-Ohmic.
- `hole_bath_reorganisation_energy`: Reorganisation energies of the bath for holes (in meV). This is used in the hole's bath spectral density, which is assumed to be super-Ohmic.
- `exciton_bath_reorganisation_energy`: Reorganisation energies of the bath for excitons (in meV). This is used in the exciton's bath spectral density, which is assumed to be super-Ohmic.
- `electron_bath_cutoff_energy`: Cutoff energies of the bath for electrons (in meV). This is used in the electron's bath spectral density, which is assumed to be super-Ohmic.
- `hole_bath_cutoff_energy`: Cutoff energies of the bath for holes (in meV). This is used in the hole's bath spectral density, which is assumed to be super-Ohmic.
- `exciton_bath_cutoff_energy`: Cutoff energies of the bath for excitons (in meV). This is used in the exciton's bath spectral density, which is assumed to be super-Ohmic.
- `T`: Temperature (in K).
- `landscape_iterations`: Number of realisations of disorder (or energetic landscapes) for simulations to be run on.
- `trajectory_iterations`: Number of trajectories to be simulated on each realisation of disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC). A number between 0 and 1, typically chosen to be 0.99. Increasing the accuracy will reduce the error introduced by dKMC's approximations but will increase the computational cost.
- `maximum_hops_cutoff`: Maximum number of hops the simulation will run before terminating. This cutoff is imposed to prevent infinite loops between, for example, two low-lying traps. A sensible choice for this cutoff will mean that only a small proportion (<5%) of trajectories terminate this way.
- `separation_cutoff`: Distance (in m) at which charges are considered separated, and the simulation terminates. This is typically chosen to be 5 nm in organic photovoltaics, representing the separation distance at which charges are unlikely to recombine.
- `exciton_population_cutoff`: Minimum population on exciton site-pairs to treat the current state as an exciton. A number between 0 and 1, typically chosen to be 0.1, which can be decreased for improved accuracy at an increased computational cost.
- `maximum_excitation_distance`: Maximum distance (in m) that an exciton can be excited from the interface. This should be chosen to be half the length of typical domain sizes.

## How to use:
1. Fill in the input file with the desired simulation parameters and save the file.
2. In terminal, run 'julia run_dKMC_charge_generation.jl charge_generation_input.in charge_generation_out.out'.
A sample output file *charge_generation_out.out* is also provided.