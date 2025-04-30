module sPTRE

#Loading the required packages.
include("../shared_functions/package_loading.jl")
include("constants.jl")

"""
    calculate_kappa(J::Function, T::Number)

Takes the spectral density as a function of ħω (in units of meV) and calculates kappa (κ), the renormalisation
constant for electronic couplings following the polaron transformation.

# Arguments:
- `J`: Spectral density as a function of energy (in meV).
- `T`: Temperature (in K).

# Output:
- `kappa`: Renormalisation constant for electronic couplings following the polaron transformation.

"""
function calculate_kappa(J::Function, T::Number)

    #Set up integral.
    F(E) = (J(E)/(E^2)) * coth(E/(2*constants.kb*T))

    #Calculate the integral.
    (a,b) = quadgk(F,0,Inf)

    #Calculate kappa.
    kappa = exp(-a)

    return kappa

end


"""
    collect_K_results(J::Function,E_step::Number,E_limit::Number,phi_step::Number,phi_limit::Number,
    bath_reorganisation_energy::Number,bath_cutoff_energy::Number,T::Number)

Collects the K results, which is a matrix containing values for K, the fourier transform of the bath correlation 
function, at different energies and for different values of the lambda function. If they are already calculated, it 
reads them from their saved location. If not, they are calculated by running calculate_K() for each of the required 
energies and values of the lambda function, before saving the full K_results in a text file for future use.

# Arguments:
- `J`: Spectral density as a function of energy (in meV).
- `E_step`: The energy step used to determine at which energies K is calculated (in meV).
- `E_limit`: The energy limit used to determine at which energies K is calculated (in meV).
- `phi_step`: The time step used in calculating phi (in fs).
- `phi_limit`: The time limit used in calculating phi (in fs).
- `bath_reorganisation_energy`: Reorganisation energy of the bath (in meV).  
- `bath_cutoff_energy`: Cutoff energy of the bath (in meV).
- `T`: Temperature (in K).
    
# Output:
- `K_results`: Matrix containing columns of values of K over energy range specified by E_limit and E_step, where each column is calculated for each of the possible values of lambda (-2 -1 0 1 2). 

"""
function collect_K_results(J::Function,E_step::Number,E_limit::Number,phi_step::Number,phi_limit::Number,bath_reorganisation_energy,bath_cutoff_energy,T::Number)

    #Create path for K data to be stored.
    DATA_K_PATH = string(dirname(dirname(@__DIR__)),"/raw_data","/K_data")
    if !isdir(dirname(DATA_K_PATH))
        mkdir(dirname(DATA_K_PATH))
    end
    if !isdir(DATA_K_PATH)
        mkdir(DATA_K_PATH)
    end

    #Create filename for the K data to be saved to.
    K_filename = "$DATA_K_PATH/E_step" * string(E_step) * "_E_limit" * string(E_limit) *"_phi_step" * string(phi_step) * "_phi_limit" * string(phi_limit) * "_λ" * string(bath_reorganisation_energy) * "_ωc" * string(bath_cutoff_energy) * "_T" * string(T) * ".csv"

    #If the file already exists, read K data from file, if not calculate it.
    if isfile(K_filename)
    
        #Read the K data from file.
        K_results = readdlm(K_filename,',',ComplexF64)
    
    else
            
        #Load phi database if already calulcated, or calculate it if not.
        phi_results = collect_phi_results(J,phi_step,phi_limit,bath_reorganisation_energy,bath_cutoff_energy,T)

        #Energies that K is calculated for.
        energies = -E_limit:E_step:E_limit #meV
        
        #Create an empty matrix for K_values at each lambda to be stored.
        K_results = complex(zeros(length(energies),5))

        #Fill the K_results matrix with values of K at given energies and lambdas.
        for lambda in [-2, -1, 1, 2]
            for i = eachindex(energies)
                K_results[i,lambda+3] = calculate_K(phi_limit,phi_step,energies[i],lambda,phi_results)
            end
        end

        #Save K_results to file.
        writedlm(K_filename,K_results,',')
    
    end

    return K_results

end


"""
    calculate_K(phi_limit::Number,phi_step::Number,E::Number,lambda::Integer,phi_results::Vector{ComplexF64})

Calculates K, the fourier transform of the bath correlation function, used as part of the sPTRE. Requires 
precalculated phi_results.

# Arguments:
- `phi_step`: The time step used in calculating phi (in fs).
- `phi_limit`: The time limit used in calculating phi (in fs).
- `E`: Energy (in meV).
- `lambda`: Integer value between -2 and 2 resulting from the lambda function in sPTRE.
- `phi_results`: Vector containing complex values of phi at times in the range specified by phi_limit and phi_step.
        
# Output:
- `K:`: The fourier transform of the bath correlation function at energy E.

"""
function calculate_K(phi_limit::Number,phi_step::Number,E::Number,lambda::Integer,phi_results::Vector{ComplexF64})
    
    #Specifies the timerange we integrate over. Phi_limit needs to be large enough that K decays to zero over this timerange, and phi_limit small enough for accurate numerical integration.
    t = collect(0:phi_step:phi_limit) 

    #Set up the integral.
    y = exp.(im * E .* t ./ constants.hbar) .*  (exp.(-lambda .* phi_results) .- 1)

    #Integrate the expression using the trapezoidal rule.
    K = NumericalIntegration.integrate(t,y)
    
    return K

end


"""
    collect_phi_results(J::Function,phi_step::Number,phi_limit::Number,bath_reorganisation_energy:::Number,
    bath_cutoff_energy::Number,T::Number)

Collects the phi_results, which is a vector containing phi values for times in the range specified by phi_limit and 
phi_step. If they are already calculated, it reads them from their saved location. If not, they are calculated by 
running calculate_phi() for each of the required times, before saving the full phi_results in a text file.

# Arguments:
- `J`: Spectral density as a function of energy (in meV).
- `phi_step`: The time step used in calculating phi (in fs).
- `phi_limit`: The time limit used in calculating phi (in fs).
- `bath_reorganisation_energy`: Reorganisation energy of the bath (in meV).  
- `bath_cutoff_energy`: Cutoff energy of the bath (in meV).
- `T`: Temperature (in K).
        
# Output:
- `phi_results`: Vector containing complex values of phi at times in the range specified by phi_limit and phi_step.

"""
function collect_phi_results(J::Function,phi_step::Number,phi_limit::Number,bath_reorganisation_energy,bath_cutoff_energy,T::Number)

    #Create path for phi data to be stored.
    DATA_PHI_PATH = string(dirname(dirname(@__DIR__)),"/raw_data","/phi_data")
    if !isdir(DATA_PHI_PATH)
        mkdir(DATA_PHI_PATH)
    end

    #Create filename for the phi data to be saved to.
    phi_filename = "$DATA_PHI_PATH/phi_step" * string(phi_step) * "_phi_limit" * string(phi_limit)  * "_λ" * string(bath_reorganisation_energy) * "_ωc" * string(bath_cutoff_energy) * "_T" * string(T) * ".csv"

    #If the file already exists, read K data from file, if not calculate it.
    if isfile(phi_filename)
    
        #Read the phi_results from file.
        phi_results = vec(readdlm(phi_filename,',',ComplexF64))
    
    else
        
        #Set time range at which phi is evaluated.
        times = collect(0:phi_step:phi_limit)

        #Calculate phi at the chosen times.
        phi_results = calculate_phi.(J,times,T)
        
        #Save phi results to file.
        writedlm(phi_filename,phi_results,',')

    end

    return phi_results

end


"""
    calculate_phi(J::Function, t::Number, T::Number)

Calculates the ϕ(t) function used as part of the sPTRE at a chosen temperature and time. 

# Arguments:
- `J`: Spectral density as a function of energy (in meV).
- `t`: Time (in fs).
- `T`: Temperature (in K).
    
# Output:
- `phi`: The phi value at time t required for calculating sPTRE rates for dKMC.

"""
function calculate_phi(J::Function, t::Number, T::Number)

    #Set up integral.
    F(E) = (J(E)./E.^2) .* (cos.((E.*t)./constants.hbar) .* coth.(E./(2*constants.kb*T)) - im .* sin.((E.*t)./constants.hbar))

    #Calculate the integral numerically using QuadGK.
    phi,phi_error = quadgk(F,0,Inf)

    return phi

end

end