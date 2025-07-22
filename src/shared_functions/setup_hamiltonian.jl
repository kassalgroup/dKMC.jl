module setup_hamiltonian

#Loading the required packages.
include("../shared_functions/package_loading.jl")
include("constants.jl")


"""
    current_charge_transport_hamiltonian(dimension::Integer,N::Integer,site_energies::Vector{<:AbstractFloat},
    electronic_coupling::Number,bath_reorganisation_energy::Number,kappa::AbstractFloat,hamiltonian_radius::Number,
    current_location::Vector{<:AbstractFloat})

Creates a hamiltonian for a subset of sites within a hamiltonian radius of the current location.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `site_energies`: List of energies of every site, of length N^dimension.
- `electronic_coupling`: Nearest neighbour electronic coupling (in meV).
- `bath_reorganisation_energy`: Reorganisation energy of the bath (in meV).
- `kappa`: Renormalisation constant for electronic couplings by the polaron transformation.
- `hamiltonian_radius`: Precalculated distance for sites that are included in Hamiltonian subset.
- `current_location`: Vector of length dimension containing the current position of the charge.

# Output:
- `H`: The system hamiltonian for the subset of sites.  
- `Ht`: The polaron-transformed system hamiltonian for the subset of sites. 
- `r`: Matrix containing the position of every site in the hamiltonian.
- `transformed_coupling`: Matrix containing only polaron-transformed coupling terms, used for dKMC rate calculations.
- `site_indexes`: Vector containing the indexes of each site included in the Hamiltonian subset.
    
"""
function current_charge_transport_hamiltonian(dimension::Integer,N::Integer,site_energies::Vector{<:AbstractFloat},electronic_coupling::Number,bath_reorganisation_energy::Number,kappa::AbstractFloat,hamiltonian_radius::Number,current_location::Vector{<:AbstractFloat})
    
    #Find independent co-ordinates that are within hamiltonian radius.
    accessible_coordinates = [Int(ceil(current_location[i]-hamiltonian_radius)):Int(floor(current_location[i]+hamiltonian_radius)) for i in 1:dimension]
    
    #Create a list of the coordinates of sites that are within the hamiltonian radius.
    site_locations = Vector{Integer}[]
    if dimension == 1
        for x in accessible_coordinates[1]
            push!(site_locations,[x])
        end
    elseif dimension == 2
        for x in accessible_coordinates[1], y in accessible_coordinates[2] 
            if (x-current_location[1])^2 + (y-current_location[2])^2 <= hamiltonian_radius^2
                push!(site_locations,[x,y])
            end
        end
    elseif dimension == 3
        for x in accessible_coordinates[1], y in accessible_coordinates[2], z in accessible_coordinates[3] 
            if (x-current_location[1])^2 + (y-current_location[2])^2 + (z-current_location[3])^2 <= hamiltonian_radius^2
                push!(site_locations,[x,y,z])
            end
        end
    end

    #Create a list of indexes to all sites within the hamiltonian subset.
    site_indexes = location_to_index.(dimension,N,site_locations)

    #Predefine the size of hamiltonian matrix and position vectors.
    hamiltonian_size = length(site_locations)
    H = zeros(hamiltonian_size,hamiltonian_size)
    r = Matrix{Integer}(zeros(hamiltonian_size,dimension))

    #Loop over all sites to fill the diagonal elements of the Hamiltonian with energies and the position matrix with site coordinates.
    for (i,site_location) in enumerate(site_locations)

        #Fill the diagonal elements of H with the appropriate site energy.
        H[i,i] = site_energies[site_indexes[i]]
        
        #Fill the position matrix with coordinates of the site.
        r[i,:] = site_location

        #Loop over all sites again to assign electronic couplings between sites to off diagonal elements of the hamiltionian.
        for (i_2,site_location_2) in enumerate(site_locations)

            #Only pairs of sites that are nearest neighbours are coupled, with strength electronic_coupling.
            if i_2 > i 
                if separation_squared(dimension,site_location,site_location_2) == 1
                    H[i,i_2] = H[i_2,i] = electronic_coupling
                end
            end
        
        end

    end

    #Create a polaron trasformed Hamiltonian with renormalised electronic couplings, and a matrix containing only polaron-transformed electronic couplings.
    Ht = H .* kappa
    transformed_coupling = copy(Ht)
    for i = eachindex(H[1,:])
        transformed_coupling[i,i] = 0
        Ht[i,i] = H[i,i] - bath_reorganisation_energy
    end

    return H,Ht,r,transformed_coupling,site_indexes

end


"""
    location_to_index(dimension::Integer,N::Integer,site_location::Vector{<:Integer})

Converts a vector containing the location of a site to the index of that site in the full lattice.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `site_location`: Vector containing the coordinates of the location of the site.

# Output:
- `index`: Index for the site in the full lattice.			

"""
function location_to_index(dimension::Integer,N::Integer,site_location::Vector{<:Integer})

    if dimension == 1
        index = site_location[1]
    elseif dimension == 2
        index = N*(site_location[1]-1) + site_location[2]
    elseif dimension == 3
        index = N^2*(site_location[1]-1) + N*(site_location[2]-1) + site_location[3]
    end

    return  index

end


"""
    separation_squared(dimension::Integer,location::Vector{<:Number},location_2::Vector{<:Number})

Calculates the separation squared between two locations.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `location`: Vector containing coordinates of location 1.
- `location_2`: Vector containing coordinates of location 2.

# Output:
- `sep_squared`: The distance squared that separates the two locations.			

"""
function separation_squared(dimension::Integer,location::Vector{<:Number},location_2::Vector{<:Number})

    if dimension == 1
        sep_squared = (location[1] - location_2[1])^2 
    elseif dimension == 2
        sep_squared = (location[1] - location_2[1])^2 + (location[2] - location_2[2])^2 
    elseif dimension == 3
        sep_squared = (location[1] - location_2[1])^2 + (location[2] - location_2[2])^2 + (location[3] - location_2[3])^2
    end

    return sep_squared

end


"""
    assign_exciton_dipole_orientations(dimension::Integer,N::Integer)

Assigns a random orientation (unit vector) to the dipole on each site in the lattice.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.

# Output:
- `dipole_orientations`: Matrix containing the orientation of the transition dipole moment of an exciton on every site. 

"""
function assign_exciton_dipole_orientations(dimension::Integer,N::Integer)

    #Create an empty matrix of the same legth as number of sites.
    dipole_orientations = zeros(N^dimension,3)

    #Assign a randomly orientated unit vector to each site.
    for i=1:N^dimension
        rand_1=1
        rand_2=1
        while rand_1^2 + rand_2^2 >= 1
            rand_1=2*rand()-1
            rand_2=2*rand()-1
        end
        dipole_orientations[i,:]=[2*rand_1*sqrt(1-rand_1^2-rand_2^2) 2*rand_2*sqrt(1-rand_1^2-rand_2^2) 1-2*(rand_1^2+rand_2^2)]
    end

    return dipole_orientations

end


"""
    dipole_coupling(transition_dipole_moment_1::Number,transition_dipole_moment_2::Number,
    dipole_orientation_1::Vector{<:AbstractFloat},dipole_orientation_2::Vector{<:AbstractFloat},
    dipole_location_1::Vector{<:Integer},dipole_location_2::Vector{<:Integer},epsilon_r::Number,site_spacing::Number)

Calculates the coupling between two dipoles under the point-dipole approximation.

# Arguments:
- `transition_dipole_moment_1`: Magnitude of the transition dipole moment on site 1 (in D).
- `transition_dipole_moment_2`: Magnitude of the transition dipole moment on site 2 (in D).
- `dipole_orientation_1`: Unit vector pointing in direction of transition dipole moment 1.
- `dipole_orientation_2`: Unit vector pointing in direction of transition dipole moment 2.
- `dipole_location_1`: Location of transition dipole moment 1 in 3D space.
- `dipole_location_2`: Location of transition dipole moment 2 in 3D space.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `epsilon_r`: Dielectric constant.

# Output:
- `coupling`: Coupling between the two dipoles under the point-dipole approximation (in meV).

"""
function dipole_coupling(transition_dipole_moment_1::Number,transition_dipole_moment_2::Number,dipole_orientation_1::Vector{<:AbstractFloat},dipole_orientation_2::Vector{<:AbstractFloat},dipole_location_1::Vector{<:Integer},dipole_location_2::Vector{<:Integer},epsilon_r::Number,site_spacing::Number)

    if length(dipole_location_1) == 1
        full_dipole_location_1 = vcat(vcat(dipole_location_1,0),0)
        full_dipole_location_2 = vcat(vcat(dipole_location_2,0),0)
    elseif length(dipole_location_1) == 2
        full_dipole_location_1 = vcat(dipole_location_1,0)
        full_dipole_location_2 = vcat(dipole_location_2,0)
    elseif length(dipole_location_1) == 3
        full_dipole_location_1 = dipole_location_1
        full_dipole_location_2 = dipole_location_2
    end

    coupling = constants.J_to_meV(((sum(dipole_orientation_1 .* dipole_orientation_2) - 3*(sum(normalize(full_dipole_location_2.-full_dipole_location_1).*dipole_orientation_1)*(sum(normalize(full_dipole_location_2.-full_dipole_location_1).*dipole_orientation_2)))) * transition_dipole_moment_1*transition_dipole_moment_2*constants.debye^2)/(4*pi*constants.epsilon_0*epsilon_r*(sqrt(sum((full_dipole_location_1 .- full_dipole_location_2).^2))*site_spacing)^3))

    return coupling

end


"""
    current_exciton_transport_hamiltonian(dimension::Integer,N::Integer,exciton_site_energies::Vector{<:AbstractFloat},
    dipole_orientations::Matrix{<:AbstractFloat},transition_dipole_moment::Number,epsilon_r::Number,
    bath_reorganisation_energy::Number,kappa::AbstractFloat,site_spacing::Number,exciton_hamiltonian_radius::Number,
    current_location::Vector{<:AbstractFloat})

Creates a hamiltonian for a subset of exciton sites within an exciton hamiltonian radius of the current location.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `exciton_site_energies`: List of energies of every site, of length N^dimension.
- `dipole_orientations`: List of unit vectors describing the orientation of the transition dipole moment of an exciton on every site.
- `transition_dipole_moment`: Magnitude of the transition dipole moment on every site (in D).
- `epsilon_r`: Dielectric constant.
- `bath_reorganisation_energy`: Reorganisation energy of the exciton bath (in meV).
- `kappa`: Renormalisation constant for excitonic couplings by the polaron transformation.
- `exciton_hamiltonian_radius`: Precalculated distance for sites that are included in Hamiltonian subset.
- `current_location`: Vector of length dimension containing the current position of the exciton.

# Output:
- `H`: The system hamiltonian for the subset of sites.  
- `Ht`: The polaron-transformed system hamiltonian for the subset of sites. 
- `r`: Matrix containing the position of every site in the hamiltonian.
- `transformed_coupling`: Matrix containing only polaron-transformed coupling terms, used for dKMC rate calculations.
- `site_indexes`: Vector containing the indexes of each site included in the Hamiltonian subset.

"""
function current_exciton_transport_hamiltonian(dimension::Integer,N::Integer,exciton_site_energies::Vector{<:AbstractFloat},dipole_orientations::Matrix{<:AbstractFloat},transition_dipole_moment::Number,epsilon_r::Number,bath_reorganisation_energy::Number,kappa::AbstractFloat,site_spacing::Number,exciton_hamiltonian_radius::Number,current_location::Vector{<:AbstractFloat})

    #Find independent co-ordinates that are within hamiltonian radius.
    accessible_coordinates = [Int.(collect(ceil(current_location[i]-exciton_hamiltonian_radius):floor(current_location[i]+exciton_hamiltonian_radius))) for i in 1:dimension]

    #Create a list of the coordinates of sites that are within the hamiltonian radius.
    site_locations = Vector{Integer}[]
    if dimension == 1
        for x in accessible_coordinates[1]
            push!(site_locations,[x])
        end
    elseif dimension == 2
        for x in accessible_coordinates[1], y in accessible_coordinates[2] 
            if (x-current_location[1])^2 + (y-current_location[2])^2 <= exciton_hamiltonian_radius^2
                push!(site_locations,[x,y])
            end
        end
    elseif dimension == 3
        for x in accessible_coordinates[1], y in accessible_coordinates[2], z in accessible_coordinates[3] 
            if (x-current_location[1])^2 + (y-current_location[2])^2 + (z-current_location[3])^2 <= exciton_hamiltonian_radius^2
                push!(site_locations,[x,y,z])
            end
        end
    end

    #Create a list of indexes to all sites within the hamiltonian subset.
    site_indexes = location_to_index.(dimension,N,site_locations)

    #Predefine the size of hamiltonian matrix and position vectors.
    hamiltonian_size = length(site_locations)
    H = zeros(hamiltonian_size,hamiltonian_size)
    r = Matrix{Integer}(zeros(hamiltonian_size,dimension))

    #Loop over all sites to fill the diagonal elements of the Hamiltonian with energies and the position matrix with site coordinates.
    for (i,site_location) in enumerate(site_locations)

        #Fill the diagonal elements of H with the appropriate site energy.
        H[i,i] = exciton_site_energies[site_indexes[i]]
        
        #Fill the position matrix with coordinates of the site.
        r[i,:] = site_location

        #Loop over all sites again to assign electronic couplings between sites to off diagonal elements of the hamiltionian.
        for (i_2,site_location_2) in enumerate(site_locations)
            if i_2 > i
                #The coupling strength is given by the dipole-dipole interaction.
                H[i,i_2] = H[i_2,i] = dipole_coupling(transition_dipole_moment,transition_dipole_moment,dipole_orientations[site_indexes[i],:],dipole_orientations[site_indexes[i_2],:],site_location,site_location_2,epsilon_r,site_spacing)
            end
        end
    end

    #Create a polaron trasformed Hamiltonian with renormalised electronic couplings, and a matrix containing only polaron-transformed electronic couplings.
    Ht = H .* kappa
    transformed_coupling = copy(Ht)
    for i = eachindex(H[1,:])
        transformed_coupling[i,i] = 0
        Ht[i,i] = H[i,i] - bath_reorganisation_energy
    end

    return H,Ht,r,transformed_coupling,site_indexes

end


"""
    charge_separation_couloumb_interaction(dimension,electron_site_location::Vector{<:Integer},
    hole_site_location::Vector{<:Integer},epsilon_r::Number,site_spacing::Number)

Calculates the coulomb interaction between the electron and hole in a site-pair.

# Arguments:
- `electron_site_location`: Coordinates for the electron's location in the site-pair.
- `hole_site_location`: Coordinates for the hole's location in the site-pair.
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).

# Output:
- `U`: The coulomb interaction strength (in meV).

"""
function charge_separation_couloumb_interaction(dimension,electron_site_location::Vector{<:Integer},hole_site_location::Vector{<:Integer},epsilon_r::Number,site_spacing::Number)

    #Calculate the seperation between the electron and hole in the site pair.
    r = sqrt(separation_squared(dimension,electron_site_location,hole_site_location))

    #Calculate the coulomb interaction strength between the charges.
    U = constants.J_to_meV(-constants.q_si^2/(4*pi*constants.epsilon_0*epsilon_r*(r*site_spacing)))

    return U

end


"""
    current_charge_separation_hamiltonian(dimension::Integer,N::Integer,energies::Vector{<:AbstractFloat},
    electronic_couplings::Vector{<:Number},epsilon_r::Number,site_spacing::Number,
    bath_reorganisation_energies::Vector{<:Number},kappas::Vector{<:AbstractFloat},
    hamiltonian_radii::Vector{<:AbstractFloat},current_electron_location::Vector{<:AbstractFloat},
    current_hole_location::Vector{<:AbstractFloat})

Creates a hamiltonian for a subset of site-pairs within hamiltonian radii of the current locations of the 
electron and hole.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `energies`: List of site energies of every site, of length N^dimension and containing donor HOMO energies followed by acceptor LUMO energies.
- `electronic_couplings`: Vector of nearest neighbour electronic couplings (in meV) for electrons in the acceptor and holes in the donor.
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `bath_reorganisation_energies`: Vector of reorganisation energies of the bath for electrons and holes (in meV).  
- `kappas`: Vector containing renormalisation constants for electron and hole couplings.
- `hamiltonian_radii`: Vector of precalculated distance for sites that are included in Hamiltonian subset.
- `current_electron_location`: Vector of length dimension containing the current position of the electron.
- `current_hole_location`: Vector of length dimension containing the current position of the hole.

# Output:
- `H`: The system hamiltonian for the subset of sites.  
- `Ht`: The polaron-transformed system hamiltonian for the subset of sites. 
- `electron_r`: Matrix containing the electron's position for every site-pair in the hamiltonian.
- `hole_r`: Matrix containing the hole's position for every site-pair in the hamiltonian.
- `transformed_coupling`: Matrix containing only polaron-transformed coupling terms, used for dKMC rate calculations.
- `electron_index`: Vector of indices to only the electron position of every site-pair.
- `hole_index`: Vector of indices to only the hole position of every site-pair.
- `bath_index`: Matrix of references to which particle type, and therefore which bath, corresponds to the coupling elements in Hamiltonian (1 for electrons, 2 for holes).
- `site_pair_indexes`: Vector containing vectors of the indexes of the electron and hole sites for the site-pairs included in the Hamiltonian subset.

"""
function current_charge_separation_hamiltonian(dimension::Integer,N::Integer,energies::Vector{<:AbstractFloat},electronic_couplings::Vector{<:Number},epsilon_r::Number,site_spacing::Number,bath_reorganisation_energies::Vector{<:Number},kappas::Vector{<:AbstractFloat},hamiltonian_radii::Vector{<:AbstractFloat},current_electron_location::Vector{<:AbstractFloat},current_hole_location::Vector{<:AbstractFloat})
    
    #Find independent co-ordinates that are within hamiltonian radius of each charge.
    accessible_electron_coordinates = [Int.(collect(ceil(current_electron_location[i]-hamiltonian_radii[1]):floor(current_electron_location[i]+hamiltonian_radii[1]))) for i in 1:dimension]
    accessible_hole_coordinates = [Int.(collect(ceil(current_hole_location[i]-hamiltonian_radii[2]):floor(current_hole_location[i]+hamiltonian_radii[2]))) for i in 1:dimension]
    
    #Remove co-ordinates where the electron is in the donor or the hole is in the acceptor.
    accessible_electron_coordinates[1] = accessible_electron_coordinates[1][findall(x->x>N/2,accessible_electron_coordinates[1])]
    accessible_hole_coordinates[1] = accessible_hole_coordinates[1][findall(x->x<=N/2,accessible_hole_coordinates[1])]

    #Create a list of the coordinates the electron and hole for site-pairs that are within the hamiltonian radius.
    site_locations = Vector{Vector{Integer}}[] 
    if dimension == 1
        for e_x in accessible_electron_coordinates[1], h_x in accessible_hole_coordinates[1]  
            if sqrt((e_x - current_electron_location[1])^2)/hamiltonian_radii[1] + sqrt((h_x - current_hole_location[1])^2)/hamiltonian_radii[2] <= 1
                push!(site_locations,[[e_x],[h_x]])
            end
        end
    elseif dimension == 2
        for e_x in accessible_electron_coordinates[1], e_y in accessible_electron_coordinates[2], h_x in accessible_hole_coordinates[1], h_y in accessible_hole_coordinates[2] 
            if sqrt((e_x - current_electron_location[1])^2 + (e_y - current_electron_location[2])^2)/hamiltonian_radii[1] + sqrt((h_x - current_hole_location[1])^2 + (h_y - current_hole_location[2])^2)/hamiltonian_radii[2] <= 1
                push!(site_locations,[[e_x,e_y],[h_x,h_y]])
            end
        end
    elseif dimension == 3
        for e_x in accessible_electron_coordinates[1], e_y in accessible_electron_coordinates[2], e_z in accessible_electron_coordinates[3], h_x in accessible_hole_coordinates[1], h_y in accessible_hole_coordinates[2], h_z in accessible_hole_coordinates[3] 
            if sqrt((e_x - current_electron_location[1])^2 + (e_y - current_electron_location[2])^2 + (e_z - current_electron_location[3])^2)/hamiltonian_radii[1] + sqrt((h_x - current_hole_location[1])^2 + (h_y - current_hole_location[2])^2 + (h_z - current_hole_location[3])^2)/hamiltonian_radii[2] <= 1
                push!(site_locations,[[e_x,e_y,e_z],[h_x,h_y,h_z]])
            end
        end
    end

    #Create a list of vectors of indexes to electron and hole sites for site-pairs within the hamiltonian subset.
    site_pair_indexes = [[location_to_index(dimension,N,site_locations[i][1]), location_to_index(dimension,N,site_locations[i][2])] for i in eachindex(site_locations)]

    #Predefine the size of hamiltonian matrix, charge position reference vectors, charge index vectors, and bath index vector.
    hamiltonian_size = length(site_locations)
    H = zeros(hamiltonian_size,hamiltonian_size)
    Ht = zeros(hamiltonian_size,hamiltonian_size)
    electron_r = Matrix{Integer}(zeros(hamiltonian_size,dimension))
    hole_r = Matrix{Integer}(zeros(hamiltonian_size,dimension))
    electron_index = Vector{Integer}(zeros(hamiltonian_size))
    hole_index = Vector{Integer}(zeros(hamiltonian_size))
    bath_index = Matrix{Integer}(zeros(hamiltonian_size,hamiltonian_size))

    #Loop over all site-pairs to fill the diagonal elements of the Hamiltonian with energies, and to record the electron and hole coordinates and electron and hole indexes for each site-pair.
    for (i,(electron_site_location,hole_site_location)) in enumerate(site_locations)

        #Fill the charge position reference vectors.
        electron_r[i,:] = electron_site_location
        hole_r[i,:] = hole_site_location

        #Fill the charge index vectors.
        electron_index[i] = site_pair_indexes[i][1]
        hole_index[i] = site_pair_indexes[i][2]

        #Assign energy of pair site to diagonal element in Hamiltonian.
        H[i,i] = energies[electron_index[i]] - energies[hole_index[i]] + charge_separation_couloumb_interaction(dimension,electron_site_location,hole_site_location,epsilon_r,site_spacing)
        
        #Calculate the energy of the site-pair following polaron transformation and assign to the diagonal element of Ht.
        Ht[i,i] = H[i,i] - bath_reorganisation_energies[1] - bath_reorganisation_energies[2]

        #Loop over all site-pairs again to assign electronic couplings between sites-pairs to off diagonal elements of the hamiltionian.
        for (i_2,(electron_site_location_2,hole_site_location_2)) in enumerate(site_locations)
            if i_2 > i
                if separation_squared(dimension,electron_site_location,electron_site_location_2) == 1 && separation_squared(dimension,hole_site_location,hole_site_location_2) == 0 #electron is adjacent only
                    bath_index[i,i_2] = bath_index[i_2,i] = 1
                    H[i,i_2] = H[i_2,i] = electronic_couplings[1]
                    Ht[i,i_2] = Ht[i_2,i] = electronic_couplings[1]*kappas[1]
                elseif separation_squared(dimension,electron_site_location,electron_site_location_2) == 0 && separation_squared(dimension,hole_site_location,hole_site_location_2) == 1 #hole is adjacent only
                    bath_index[i,i_2] = bath_index[i_2,i] = 2
                    H[i,i_2] = H[i_2,i] = electronic_couplings[2]
                    Ht[i,i_2] = Ht[i_2,i] = electronic_couplings[2]*kappas[2]
                end
            end
        end
    end

    #Create a matrix containing only polaron-transformed electronic couplings.
    transformed_coupling = copy(Ht)
    for i = eachindex(H[1,:])
        transformed_coupling[i,i] = 0
    end

    return H,Ht,electron_r,hole_r,transformed_coupling,electron_index,hole_index,bath_index,site_pair_indexes

end


"""
    LUMO_HOMO_energies(dimension::Integer,N::Integer,disorders::Matrix{<:Number},exciton_disorders::Vector{<:Number},
    donor_HOMO_LUMO_GAP::Number,HOMO_offset::Number,LUMO_offset::Number)

Assigns a HOMO and LUMO energy to every site in your system. The HOMO and LUMO energies are drawn from gaussian 
distributions that are correlated such that their difference (exciton energy) forms a gaussian distribution with a 
smaller excitonic disorder.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `disorders`: Matrix of energetic disorders (in meV) of the electrons (row 1) and holes (row 2) in the donor (column 1) and acceptor (column 2).
- `exciton_disorders`: Vector of excitonic disorders of the donor and acceptor (in meV).
- `donor_HOMO_LUMO_GAP`: Energy gap between donor HOMO and LUMO levels (in meV).
- `LUMO_offset`: Energy offset between donor and acceptor LUMO levels (in meV).
- `HOMO_offset`: Energy offset between donor and acceptor HOMO levels (in meV).

# Output:
- `LUMO_HOMO_energies`: Matrix containing LUMO (row 1) and HOMO energies (row 2) for every site in the system. These are referenced in the following way: in 1D [:,x], in 2D [:,N*(x-1)+y], or in 3D [:,N^2*(x-1)+N*(y-1)+z].

"""
function LUMO_HOMO_energies(dimension::Integer,N::Integer,disorders::Matrix{<:Number},exciton_disorders::Vector{<:Number},donor_HOMO_LUMO_GAP::Number,HOMO_offset::Number,LUMO_offset::Number)

    #Mean of the bivariate normal distribution.
    mu_donor = [donor_HOMO_LUMO_GAP; 0]
    mu_acceptor = [donor_HOMO_LUMO_GAP-LUMO_offset; -HOMO_offset]

    #Calculating the HOMO_LUMO_correlations from the charge and exciton disorders
    HOMO_LUMO_correlations = [(disorders[1,i]^2 + disorders[2,i]^2 - exciton_disorders[i]^2)/(2*disorders[1,i]*disorders[2,i]) for i in [1 2]]

    #Covariance of the bivariate normal distribution.
    sigma_donor = [disorders[1,1]^2 HOMO_LUMO_correlations[1]*disorders[1,1]*disorders[2,1]; HOMO_LUMO_correlations[1]*disorders[1,1]*disorders[2,1] disorders[2,1]^2]
    sigma_acceptor = [disorders[1,2]^2 HOMO_LUMO_correlations[2]*disorders[1,2]*disorders[2,2]; HOMO_LUMO_correlations[2]*disorders[1,2]*disorders[2,2] disorders[2,2]^2]

    #Bivariate normal distribution.
    distribution_donor = Distributions.MvNormal(mu_donor,sigma_donor)
    distribution_acceptor = Distributions.MvNormal(mu_acceptor,sigma_donor)

    #Define energies for LUMO and HOMO levels of every site.
    LUMO_HOMO_energies_donor = rand(distribution_donor,Int(N^dimension/2))
    LUMO_HOMO_energies_acceptor = rand(distribution_acceptor,Int(N^dimension/2))
    LUMO_HOMO_energies = hcat(LUMO_HOMO_energies_donor,LUMO_HOMO_energies_acceptor)

    return LUMO_HOMO_energies

end


"""
    is_site_pair_within_radii(dimension::Integer,exciton::Bool,electron_site_location::Vector{<:Integer},
    hole_site_location::Vector{<:Integer},current_electron_location::Vector{<:AbstractFloat},
    current_hole_location::Vector{<:AbstractFloat},current_exciton_location::Vector{<:AbstractFloat},
    electron_hamiltonian_radius::Number,hole_hamiltonian_radius::Number,exciton_hamiltonian_radius::Number)

Determines whether the location of a site-pair lies within the hamilonian radii of the current location of the charges.

# Arguments:
- `exciton`: Boolean stating whether the current state is excitonic (true) or not (false).
- `electron_site_location`: Coordinates for the electron's location in the site-pair.
- `hole_site_location`: Coordinates for the hole's location in the site-pair.
- `current_electron_location`: Coordinates for the expectation of the electron's location in the current state.
- `current_hole_location`: Coordinates for the expectation of the hole's location in the current state.
- `current_exciton_location`: oordinates for the expectation of the exciton's location in the current state.
- `electron_hamiltonian_radius`: Precalculated distance from the electron for sites that are included in Hamiltonian subset.
- `hole_hamiltonian_radius`: Precalculated distance from the hole for sites that are included in Hamiltonian subset.
- `exciton_hamiltonian_radius`: Precalculated distance from the exciton for sites that are included in Hamiltonian subset.

# Output:
- `outcome`: Boolean stating whether the site-pair falls within the hamiltonian radii.

"""
function is_site_pair_within_radii(dimension::Integer,exciton::Bool,electron_site_location::Vector{<:Integer},hole_site_location::Vector{<:Integer},current_electron_location::Vector{<:AbstractFloat},current_hole_location::Vector{<:AbstractFloat},current_exciton_location::Vector{<:AbstractFloat},electron_hamiltonian_radius::Number,hole_hamiltonian_radius::Number,exciton_hamiltonian_radius::Number)

    if exciton == true
        if sqrt(separation_squared(dimension,electron_site_location,current_electron_location))/electron_hamiltonian_radius + sqrt(separation_squared(dimension,hole_site_location,current_hole_location))/hole_hamiltonian_radius <= 1 || (separation_squared(dimension,electron_site_location,hole_site_location) == 0 && separation_squared(dimension,electron_site_location,current_electron_location)  <= exciton_hamiltonian_radius^2)
            outcome = true
        else
            outcome = false
        end
    else 
        if sqrt(separation_squared(dimension,electron_site_location,current_electron_location))/electron_hamiltonian_radius + sqrt(separation_squared(dimension,hole_site_location,current_hole_location))/hole_hamiltonian_radius <= 1 				
            outcome = true
        else
            outcome = false
        end
    end

    return outcome

end


"""
    charge_generation_coulomb_interaction(dimension::Integer,N::Integer,electron_site_location::Vector{<:Integer},
    hole_site_location::Vector{<:Integer},exciton_binding_energies::Vector{<:Number},epsilon_r::Number,
    site_spacing::Number)

Calculates the coulomb interaction between the electron and hole in a site-pair of the full charge charge generation 
hamiltonian.

# Arguments:
- `N`: Length of system, i.e., number of sites in each direction.
- `electron_site_location`: Coordinates for the electron's location in the site-pair.
- `hole_site_location`: Coordinates for the hole's location in the site-pair.
- `exciton_binding_energies`: Vector of exciton binding energies for the donor and acceptor (in meV).
- `epsilon_r`: Dielectric constant.
- `site_spacing`: The distance between sites of the cubic lattice (in m).

# Output:
- `U`: The coulomb interaction strength (in meV).
"""
function charge_generation_coulomb_interaction(dimension::Integer,N::Integer,electron_site_location::Vector{<:Integer},hole_site_location::Vector{<:Integer},exciton_binding_energies::Vector{<:Number},epsilon_r::Number,site_spacing::Number)

    #Calculate the seperation between the electron and hole in the site pair.
    r = sqrt(separation_squared(dimension,electron_site_location,hole_site_location))
    
    #If separation is zero return the exciton binding energy, if not return the coulomb interaction strength.
    if r == 0 
        if electron_site_location[1] <= N/2
            U = -exciton_binding_energies[1]
        else
            U = -exciton_binding_energies[2]
        end
    else
        U = constants.J_to_meV(-constants.q_si^2/(4*pi*constants.epsilon_0*epsilon_r*(r*site_spacing)))
    end

    return U

end


"""
    current_charge_generation_hamiltonian(dimension::Integer,N::Integer,exciton::Bool,
    LUMO_HOMO_energies::Matrix{<:AbstractFloat},exciton_binding_energies::Vector{<:Number},
    electronic_couplings::Matrix{<:Number},transition_dipole_moments::Vector{<:Number},
    dipole_orientations::Matrix{<:AbstractFloat},epsilon_r::Number,site_spacing::Number,
    bath_reorganisation_energies::Vector{<:Number},kappas::Vector{<:AbstractFloat},
    hamiltonian_radii::Matrix{<:AbstractFloat},current_electron_location::Vector{<:AbstractFloat},
    current_hole_location::Vector{<:AbstractFloat},current_exciton_location::Vector{<:AbstractFloat})

Creates a hamiltonian for a subset of site-pairs within hamiltonian radii of the current locations of the electron, 
hole, and exciton.

# Arguments:
- `dimension`: Dimension of the system (1, 2, or 3).
- `N`: Length of system, i.e., number of sites in each direction.
- `LUMO_HOMO_energies`: Matrix of LUMO (row 1) and HOMO energies (row 2) on every site (in meV).
- `exciton_binding_energies`: Vector of exciton binding energies for the donor and acceptor (in meV).
- `electronic_couplings`: Matrix of nearest neighbour electronic couplings (in meV) for electrons (row 1) and holes (row 2) in the donor (column 1), at the interface (column 2), and in the acceptor (column 3).
- `transition_dipole_moments`: Vector of magnitude of the transition dipole moments on donor and acceptor sites (in D).
- `dipole_orientations`: List of unit vectors describing the orientation of the transition dipole moment of an exciton on every site in the system. 
- `epsilon_r`: Vector of dielectric constants for the donor and acceptor materials.
- `site_spacing`: The distance between sites of the cubic lattice (in m).
- `bath_reorganisation_energies`: Vector of reorganisation energies of the bath for electrons and holes (in meV).  
- `kappas`: Vector containing renormalisation constants for electron, hole, and exciton couplings.
- `hamiltonian_radii`: Vector of precalculated distance for sites that are included in Hamiltonian subset.
- `current_electron_location`: Vector of length dimension containing the current position of the electron.
- `current_hole_location`: Vector of length dimension containing the current position of the hole.
- `current_exciton_location`: Vector of length dimension containing the current position of the exciton.

# Output:
- `H`: The system hamiltonian for the subset of sites.  
- `Ht`: The polaron-transformed system hamiltonian for the subset of sites. 
- `electron_r`: Matrix containing the electron's position for every site-pair in the hamiltonian.
- `hole_r`: Matrix containing the hole's position for every site-pair in the hamiltonian.
- `dipoles`: Matrix containing the orientation of the transition dipole moment on any exciton site-pairs in the hamiltonian.
- `transformed_coupling`: Matrix containing only polaron-transformed coupling terms, used for dKMC rate calculations.
- `electron_index`: Vector of indices to only the electron position of every site-pair.
- `hole_index`: Vector of indices to only the hole position of every site-pair.
- `exciton_index`: Vector of indices to only the exciton position of every site-pair.
- `bath_index`: Matrix of references to which particle type, and therefore which bath, corresponds to the coupling elements in Hamiltonian (1 for electrons, 2 for holes, 3 for excitons).
- `site_pair_indexes`: Vector containing vectors of the indexes of the electron and hole sites for the site-pairs included in the Hamiltonian subset.

"""
function current_charge_generation_hamiltonian(dimension::Integer,N::Integer,exciton::Bool,LUMO_HOMO_energies::Matrix{<:AbstractFloat},exciton_binding_energies::Vector{<:Number},electronic_couplings::Matrix{<:Number},transition_dipole_moments::Vector{<:Number},dipole_orientations::Matrix{<:AbstractFloat},epsilon_r::Number,site_spacing::Number,bath_reorganisation_energies::Vector{<:Number},kappas::Vector{<:AbstractFloat},hamiltonian_radii::Matrix{<:AbstractFloat},current_electron_location::Vector{<:AbstractFloat},current_hole_location::Vector{<:AbstractFloat},current_exciton_location::Vector{<:AbstractFloat})

    #Finding which hamiltonian_radii to use for each particle.
    if current_electron_location[1] <= N/2
        electron_hamiltonian_radius = hamiltonian_radii[1,1]
        exciton_hamiltonian_radius = hamiltonian_radii[3,1]
    else
        electron_hamiltonian_radius = hamiltonian_radii[1,2]
        exciton_hamiltonian_radius = hamiltonian_radii[3,2]
    end

    if current_hole_location[1] <= N/2
        hole_hamiltonian_radius = hamiltonian_radii[2,1]
    else
        hole_hamiltonian_radius = hamiltonian_radii[2,2]
    end

    #Find independent co-ordinates that are within hamiltonian radius of each charge.
    accessible_electron_coordinates = [Int.(collect(ceil(current_electron_location[i]-electron_hamiltonian_radius):floor(current_electron_location[i]+electron_hamiltonian_radius))) for i in 1:dimension]
    accessible_hole_coordinates = [Int.(collect(ceil(current_hole_location[i]-hole_hamiltonian_radius):floor(current_hole_location[i]+hole_hamiltonian_radius))) for i in 1:dimension]

    #If we are in an exciton state, we also include exciton sites within the exciton Hamiltonian radius.
    if exciton == true 			
        accessible_exciton_coordinates = [Int.(collect(ceil(current_exciton_location[i]-exciton_hamiltonian_radius):floor(current_exciton_location[i]+exciton_hamiltonian_radius))) for i in 1:dimension]
        accessible_electron_coordinates = [sort(unique(vcat(accessible_electron_coordinates[i],accessible_exciton_coordinates[i]))) for i in 1:dimension]
        accessible_hole_coordinates = [sort(unique(vcat(accessible_hole_coordinates[i],accessible_exciton_coordinates[i]))) for i in 1:dimension]
    end

    #Create a list of the locations of the electron and hole for all site-pairs within the hamiltonian radii.
    site_locations = Vector{Vector{Integer}}[]
    if dimension == 1
        for e_x in accessible_electron_coordinates[1], h_x in accessible_hole_coordinates[1]
            if is_site_pair_within_radii(dimension,exciton,[e_x],[h_x],current_electron_location,current_hole_location,current_exciton_location,electron_hamiltonian_radius,hole_hamiltonian_radius,exciton_hamiltonian_radius)
                push!(site_locations,[[e_x],[h_x]])
            end
        end
    elseif dimension == 2
        for e_x in accessible_electron_coordinates[1], e_y in accessible_electron_coordinates[2], h_x in accessible_hole_coordinates[1], h_y in accessible_hole_coordinates[2]
            if is_site_pair_within_radii(dimension,exciton,[e_x,e_y],[h_x,h_y],current_electron_location,current_hole_location,current_exciton_location,electron_hamiltonian_radius,hole_hamiltonian_radius,exciton_hamiltonian_radius)
                push!(site_locations,[[e_x,e_y],[h_x,h_y]])
            end
        end				
    elseif dimension == 3
        for e_x in accessible_electron_coordinates[1], e_y in accessible_electron_coordinates[2], e_z in accessible_electron_coordinates[3], h_x in accessible_hole_coordinates[1], h_y in accessible_hole_coordinates[2], h_z in accessible_hole_coordinates[3]
            if is_site_pair_within_radii(dimension,exciton,[e_x,e_y,e_z],[h_x,h_y,h_z],current_electron_location,current_hole_location,current_exciton_location,electron_hamiltonian_radius,hole_hamiltonian_radius,exciton_hamiltonian_radius)
                push!(site_locations,[[e_x,e_y,e_z],[h_x,h_y,h_z]])
            end
        end
    end

    #Create a list of vectors of indexes to electron and hole sites for site-pairs within the hamiltonian subset.
    site_pair_indexes = [[location_to_index(dimension,N,site_locations[i][1]), location_to_index(dimension,N,site_locations[i][2])] for i in eachindex(site_locations)]

    #Creating empty hamiltonian matrices, charge position reference matrices, dipole orientation matrices, charge index vectors, and bath index vector. 
    hamiltonian_size = length(site_locations)
    H = zeros(hamiltonian_size,hamiltonian_size)
    Ht = zeros(hamiltonian_size,hamiltonian_size)
    electron_r = Matrix{Integer}(zeros(hamiltonian_size,dimension))
    hole_r = Matrix{Integer}(zeros(hamiltonian_size,dimension))
    dipoles = zeros(hamiltonian_size,3)
    electron_index = Vector{Integer}(undef,hamiltonian_size)
    hole_index = Vector{Integer}(undef,hamiltonian_size)
    exciton_index = Vector{Integer}(zeros(hamiltonian_size))
    bath_index = Matrix{Integer}(zeros(hamiltonian_size,hamiltonian_size))

    #Loop over all site-pairs to fill the diagonal elements of the Hamiltonian with energies, and to record the electron, hole, and exciton coordinates and electron, hole, and exciton indexes for each site-pair.
    for (i,(electron_site_location,hole_site_location)) in enumerate(site_locations)

        #Fill the charge position and index reference vectors.
        electron_r[i,:] = electron_site_location
        hole_r[i,:] = hole_site_location
        electron_index[i] = site_pair_indexes[i][1]
        hole_index[i] = site_pair_indexes[i][2]

        #If the site-pair represents an exciton, save the dipole orientation to dipole vector.
        separation = sqrt(separation_squared(dimension,electron_site_location,hole_site_location))
        if separation == 0
            exciton_index[i] = site_pair_indexes[i][1]
            dipoles[i,:] = dipole_orientations[exciton_index[i],:]
        end

        #Assign energy of pair site to diagonal element in Hamiltonian.
        E = LUMO_HOMO_energies[1,site_pair_indexes[i][1]] - LUMO_HOMO_energies[2,site_pair_indexes[i][2]] + charge_generation_coulomb_interaction(dimension,N,electron_site_location,hole_site_location,exciton_binding_energies,epsilon_r,site_spacing)
        H[i,i] = E

        #Calculate the energy of the site-pair following polaron transformation and assign to the diagonal element of Ht.
        Ht[i,i] = E
        if separation == 0
            Ht[i,i] = Ht[i,i] - bath_reorganisation_energies[3]
        else
            Ht[i,i] = Ht[i,i] - bath_reorganisation_energies[1] - bath_reorganisation_energies[2]
        end

        #Loop over all site-pairs again to assign electronic couplings between sites-pairs to off diagonal elements of the hamiltionian.
        for (i_2,(electron_site_location_2,hole_site_location_2)) in enumerate(site_locations)

            #Only pairs of site-pairs where both are excitons or either the electron or hole are nearest neighbours are coupled.
            if i_2 > i
                separation_2 = sqrt(separation_squared(dimension,electron_site_location_2,hole_site_location_2))
                if separation_squared(dimension,electron_site_location,electron_site_location_2) == 1 && separation_squared(dimension,hole_site_location,hole_site_location_2) == 0 #If only electron is adjacent.
                    bath_index[i,i_2] = bath_index[i_2,i] = 1
                    if electron_site_location[1] <= N/2 && electron_site_location_2[1] <= N/2
                        H[i,i_2] = H[i_2,i] =  electronic_couplings[1,1]
                    elseif electron_site_location[1] > N/2 && electron_site_location_2[1] > N/2
                        H[i,i_2] = H[i_2,i] =  electronic_couplings[1,3]
                    else
                        H[i,i_2] = H[i_2,i] =  electronic_couplings[1,2]
                    end
                    Ht[i,i_2] = Ht[i_2,i] =  kappas[1]*H[i,i_2]
                elseif separation_squared(dimension,electron_site_location,electron_site_location_2) == 0 && separation_squared(dimension,hole_site_location,hole_site_location_2) == 1 #If only hole is adjacent.
                    bath_index[i,i_2] = bath_index[i_2,i] = 2
                    if hole_site_location[1] <= N/2 && hole_site_location_2[1] <= N/2
                        H[i,i_2] = H[i_2,i] =  electronic_couplings[2,1]
                    elseif hole_site_location[1] > N/2 && hole_site_location_2[1] > N/2
                        H[i,i_2] = H[i_2,i] =  electronic_couplings[2,3]
                    else
                        H[i,i_2] = H[i_2,i] =  electronic_couplings[2,2]
                    end
                    Ht[i,i_2] = Ht[i_2,i] =  kappas[2]*H[i,i_2]
                elseif separation == 0 && separation_2 == 0 #If both site-pairs are excitons.
                    bath_index[i,i_2] = bath_index[i_2,i] = 3
                    if electron_site_location[1] <= N/2
                        transition_dipole_moment_1 = transition_dipole_moments[1]
                    else
                        transition_dipole_moment_1 = transition_dipole_moments[2]
                    end
                    if electron_site_location_2[1] <= N/2
                        transition_dipole_moment_2 = transition_dipole_moments[1]
                    else
                        transition_dipole_moment_2 = transition_dipole_moments[2]
                    end
                    H[i,i_2] = H[i_2,i] =  dipole_coupling(transition_dipole_moment_1,transition_dipole_moment_2,dipole_orientations[site_pair_indexes[i][1],:],dipole_orientations[site_pair_indexes[i_2][1],:],electron_site_location,electron_site_location_2,site_spacing,epsilon_r)
                    Ht[i,i_2] = Ht[i_2,i] =  kappas[3]*H[i,i_2]
                end
            end
        end
    
    end 

    #Create a matrix containing only polaron-transformed electronic couplings.
    transformed_coupling = copy(Ht)
    for i=eachindex(Ht[1,:])
        transformed_coupling[i,i] = 0
    end

    return H,Ht,electron_r,hole_r,dipoles,transformed_coupling,electron_index,hole_index,exciton_index,bath_index,site_pair_indexes
    
end


"""
    compute_centres(evecs::AbstractMatrix,r::AbstractMatrix)

Calculates the centres, expectation value of position, of the eigenstates.

# Arguments:
- `evecs`: Energy eigenvectors of the polaron trasformed system hamiltonian.
- `r`: Matrix containing the position of every site in the hamiltonian.

# Output:
- `centres`: Matrix containing the coordinates of the centres of eigenstates, the expectation of the eigenstates position.

"""
function compute_centres(dimension::Integer,evecs::AbstractMatrix,r::AbstractMatrix)

    n_sites, n_states = size(evecs)
    centres = Matrix{Float64}(undef, n_states, dimension)
    @inbounds for j in 1:dimension
        for i in 1:n_states
            acc = 0.0
            @simd for k in 1:n_sites
                acc += r[k, j] * abs2(evecs[k, i])
            end
            centres[i, j] = acc
        end
    end
    
    return centres

end

end

