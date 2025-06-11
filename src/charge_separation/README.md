# dKMC Charge Separation:
This dKMC module simulates the separation of a partially delocalised electron-hole pair from an interfacial charge-transfer (CT) state to free charges. The simulations start from a delocalised interfacial CT state, i.e., an electron and a hole on the opposite sides of an interface between the electron-donor and electron-acceptor materials. The simulation then uses dKMC to propagate the two-particle separation dynamics, with the electron constrained to the acceptor and the hole to the donor, until the charges either recombine or separate. The simulation takes in parameters describing both materials and returns the internal quantum efficiency (IQE) of charge separation, as well as the mean characteristics (energy, separation, and IPR) of the initial and final states.

For details, see:
[Balzer, D.; Kassal, I. Even a Little Delocalisation Produces Large Kinetic Enhancements of Charge-Separation Efficiency in Organic Photovoltaics. *Science Advances* **2022**, *8*, eabl9692.](https://www.science.org/doi/10.1126/sciadv.abl9692)

## Input File:
The input parameters are to be supplied in YAML format; for an example, see `charge_separation_input.in` included in this folder

Parameters required in the input file:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `acceptor_electron_disorder`: Energetic disorder (in meV) of electrons in the acceptor material.
- `donor_hole_disorder`: Energetic disorder (in meV) of holes in the donor material.
- `donor_HOMO_acceptor_LUMO_GAP`: Energy gap between donor HOMO and acceptor LUMO levels (in meV).
- `acceptor_electron_coupling`: Nearest neighbour electronic couplings (in meV) for electrons in the acceptor.
- `donor_hole_coupling`: Nearest neighbour electronic couplings (in meV) for holes in the donor.
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `CT_lifetime`: CT state lifetime (in s).
- `electron_bath_reorganisation_energy`: Reorganisation energy of the bath for electrons (in meV). This is used in the electron's bath spectral density, which is assumed to be super-Ohmic.
- `hole_bath_reorganisation_energy`: Reorganisation energy of the bath for holes (in meV). This is used in the hole's bath spectral density, which is assumed to be super-Ohmic.
- `electron_bath_cutoff_energy`: Cutoff energy of the bath for electrons (in meV). This is used in the electron's bath spectral density, which is assumed to be super-Ohmic.
- `hole_bath_cutoff_energy`: Cutoff energy of the bath for holes (in meV). This is used in the hole's bath spectral density, which is assumed to be super-Ohmic.
- `T`: Temperature (in K).
- `landscape_iterations`: Number of realisations of disorder (or energetic landscapes) for simulations to be run on.
- `trajectory_iterations`: Number of trajectories to be simulated on each realisation of disorder.
- `accuracy`: The accuracy of dKMC calculations (a_dKMC). A number between 0 and 1, typically chosen to be 0.99. Increasing the accuracy will reduce the error introduced by dKMC's approximations but will increase the computational cost.
- `maximum_hops_cutoff`: Maximum number of hops the simulation will run before terminating. This cutoff is imposed to prevent infinite loops between, for example, two low-lying traps. A sensible choice for this cutoff will mean that only a small proportion (<5%) of trajectories terminate this way.
- `separation_cutoff`: Distance (in m) at which charges are considered separated, and the simulation terminates. This is typically chosen to be 5 nm in organic photovoltaics, representing the separation distance at which charges are unlikely to recombine.

## How to use:
1. Fill in the input file with the desired simulation parameters and save the file.
2. In terminal, run 'julia run_dKMC_charge_separation.jl charge_separation_input.in charge_separation_out.out'.
A sample output file *charge_separation_out.out* is also provided.