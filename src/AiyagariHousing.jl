module AiyagariHousing

using Parameters, QuantEcon, Distributions, NLopt, CompEcon

import LinearAlgebra.mul!

include("Types.jl")
include("Functions.jl")

export Pars, u, u_ind, Psi, myfunc, myconstraint, @unpack_Pars, lom_housing, vmax, myfunc2, hello_Irene
    #lom_assets, lom_housing, find_EV, find_c, Val_of_c, EV_of_s, find_c_constr
end # module
