module AiyagariHousing

using Parameters, QuantEcon, Distributions, CompEcon, Printf, JuMP, Setfield, AutomaticDocstrings, Roots, NLsolve
import LinearAlgebra.mul!, Ipopt

include("Types.jl")
include("Functions.jl")

export Pars, u, u_ind, Psi, @unpack_Pars, lom_housing, myfunc2, solve_model, setup_Q, Prob_a_next, calc_step, delta_low, delta_high, marketclearing_housing, find_rental_rate, solve_GE, stst_distr, comp_statics, risk_premium
#lom_assets, lom_housing, find_EV, find_c, Val_of_c, EV_of_s, find_c_constr
end # module