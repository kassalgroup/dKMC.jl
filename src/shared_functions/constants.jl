
module constants

#This module includes a collection of universal constants required for dKMC simulations.

#Set constants
const hbar_si = 1.0545718e-34  #Js
const kb_si = 1.3806485e-23  #J/K
const q_si = 1.6021766208e-19 #C
const epsilon_0 = 8.8541878128e-12 #F.m^-1
const debye = 3.33564e-30 #C.m

#Set conversion functions
to_milli(x) = x * 1e+3
to_micro(x) = x * 1e+6
to_nano(x) = x * 1e+9
to_pico(x) = x * 1e+12
to_femto(x) = x * 1e+15
to_kilo(x) = x * 1e-3
to_mega(x) = x * 1e-6
to_giga(x) = x * 1e-9
from_milli(x) = x * 1e-3
from_micro(x) = x * 1e-6
from_nano(x) = x * 1e-9
from_pico(x) = x * 1e-12
from_femto(x) = x * 1e-15
from_kilo(x) = x * 1e+3
from_mega(x) = x * 1e+6
from_giga(x) = x * 1e+9
J_to_eV(x) = x / q_si
eV_to_J(x) = x * q_si
J_to_meV(x) = to_milli(J_to_eV(x))
meV_to_J(x) = from_milli(eV_to_J(x))
m2_to_cm2(x) = x * 1e+4
m2_to_nm2(x) = x * 1e+18
nm2_to_cm2(x) = x * 1e-14 
inv_cm_to_meV(x) = x * 0.123983

#Constants in non-SI units 
const kb = J_to_meV(kb_si)  #meV/K
const hbar = to_femto(J_to_meV(hbar_si)) #meV⋅fs

end