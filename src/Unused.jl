
# function myfunc3(x::Vector, grad::Vector, a, z_i, par::Pars, h_next, a_next, V_low, V_high, V_low_w, V_high_w, util)
#     # We optimize by choosing two types of investment: Housing and Capital.
#     # Consumption is then simply the residual variable.
#     h = x[1] # new housing
#     K = x[2] # new capital
#     # a includes already undepreciated capital and housing
#     z = par.z_chain.state_values[z_i]
#     c = a + par.w * z - (par.P_h - par.P_l) * h - K
#     #! Assume for a moment that income is constant
#     h_next .= lom_housing.(h, par.δs)
#     a_next .= par.R .* K .+ h_next

#     V_low .= funeval(par.c_low, par.basis, a_next)[1][:]
#     V_high .= funeval(par.c_high, par.basis, a_next)[1][:]
    
#     mul!(V_low_w, par.w_δt, V_low, par.z_chain.p[z_i, 1], 0)
#     mul!(V_high_w, par.w_δt, V_high, par.z_chain.p[z_i, 2], 0)

#     util .= -(u_ind(c; par) .+ par.β .* (V_low_w .+ V_high_w))
#     #! Maximize expected utlity.
#     #return -(u_ind(c; par) + β * sum((z_chain.p[z_i, 1] .* w_δ .* funeval(c_low, basis, a_next)[1][:] .+ z_chain.p[z_i, 2] .* w_δ .* funeval(c_high, basis, a_next)[1][:])))
#     return util
# end

# function myfunc2(x::Vector, grad::Vector, a, z_i, par::Pars, h_next, a_next)
#     @unpack_Pars par
#     # We optimize by choosing two types of investment: Housing and Capital.
#     # Consumption is then simply the residual variable.
#     h = x[1] # new housing
#     K = x[2] # new capital
#     # a includes already undepreciated capital and housing
#     z = z_chain.state_values[z_i]
#     c = a + w * z - (P_h - P_l) * h - K
#     #! Assume for a moment that income is constant
#     h_next .= lom_housing.(h, δs)
#     a_next .= R .* K .+ h_next

#     #! Maximize expected utlity.
#     return -(u_ind(c; par) + β * sum((z_chain.p[z_i, 1] .* w_δ .* funeval(c_low, basis, a_next)[1][:] .+ z_chain.p[z_i, 2] .* w_δ .* funeval(c_high, basis, a_next)[1][:])))
# end