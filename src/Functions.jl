"""
    CRRA utility function over consumption and housing
"""
u(c, h; par) = (c^par.θ * h^(1 - par.θ))^(1 - par.σ) / (1 - par.σ)

"""
    Indirect utility function as in Jeske, Krueger, Mitman
    where c = c_tilde + r*h is equal to total expenditure
"""
function u_ind(c; par)
    (Psi(par) * c^(1 - par.σ) - 1) / (1 - par.σ)
end

function Psi(par)
    @unpack_Pars par
    (θ^θ * (1 - θ)^(1 - θ) * r^(θ - 1))^(1 - σ)
end

"""
    Law of motion for housing investments
"""
lom_housing(h, δ_shock) = (1 - δ_shock) * h

lom_assets(b, g, δ_shock, y) = b + lom_housing(g, δ_shock) + y


"""
    The function to be optimized.
"""
function myfunc(x::Vector, grad::Vector, a, z_i, par::Pars, h_next, a_next, V_low_w, V_high_w)
    @unpack_Pars par
    # We optimize by choosing two types of investment: Housing and Capital.
    # Consumption is then simply the residual variable.

    h = x[1] # new housing
    K = x[2] # new capital

    # a includes already undepreciated capital and housing

    z = z_chain.state_values[z_i]
    c = a + w * z - (P_h - P_l) * h - K

    #! Assume for a moment that income is constant
    h_next .= lom_housing.(h, δs)
    a_next .= R .* K .+ h_next
    
    mul!(V_low_w, w_δt, funeval(c_low, basis, a_next)[1][:], z_chain.p[z_i, 1], 0)
    mul!(V_high_w, w_δt, funeval(c_high, basis, a_next)[1][:], z_chain.p[z_i, 2], 0)

    #! Maximize expected utlity.
    return -(u_ind(c; par) .+ β .* (V_low_w .+ V_high_w))[1]
end

"""
    Budget constraint.
"""
function myconstraint(x::Vector, grad::Vector, a, z_i, par::Pars)
    @unpack_Pars par
    z = z_chain.state_values[z_i]
    h = x[1]
    K = x[2]
    #! Basically: Investment needs to be less than total assets
    return (P_h - P_l) * h + K - a - w * z # ≤ 0
end

function vmax(a, z_i; par, h_next, a_next, V_low_w, V_high_w, x0 = [0.04, 0.04])
    # Two variables, choose gradient-free algorithm
    opt = Opt(:LN_COBYLA, 2)
    opt.lower_bounds = [0.0, 0.0]
    opt.xtol_rel = 1e-4
    # define objective
    myfinalfunc(x, g) = myfunc(x, g, a, z_i, par, h_next, a_next, V_low_w, V_high_w)
    opt.min_objective = myfinalfunc
    #define constraint
    inequality_constraint!(opt, (x, g) -> myconstraint(x, g, a, z_i, par), 1e-8)
    (minf, minx, ret) = optimize(opt, x0)
    
    return minf, minx, ret
end

function vmax(a, z_i; par, h_next, a_next, V_low_w, V_high_w, x0 = [0.04, 0.04])
    # Two variables, choose gradient-free algorithm
    opt = Opt(:LN_COBYLA, 2)
    opt.lower_bounds = [0.0, 0.0]
    opt.xtol_rel = 1e-4
    # define objective
    myfinalfunc(x, g) = myfunc(x, g, a, z_i, par, h_next, a_next, V_low_w, V_high_w)
    opt.min_objective = myfinalfunc
    #define constraint
    inequality_constraint!(opt, (x, g) -> myconstraint(x, g, a, z_i, par), 1e-8)
    (minf, minx, ret) = optimize(opt, x0)
    
    return minf, minx, ret
end

#! Additional constraint to add: Minimum buying requirement.
"""
    The function to be optimized without gradient.
"""
function myfunc2(h, K, a, z_i, par::Pars)
    @unpack_Pars par
    # We optimize by choosing two types of investment: Housing and Capital.
    # Consumption is then simply the residual variable.

    # a includes already undepreciated capital and housing

    z = z_chain.state_values[z_i]
    if housing_beginning == :yes
        c = a + w * z - (P_h - P_l) * h - K
        h_next = lom_housing.(h, δs)
    elseif housing_beginning == :no
        c = a + w * z - P_h * h - K
        h_next = lom_housing.(h, δs) .+ h .* P_l
    end

    a_next = R .* K .+ h_next

    V_low_w = w_δt * funeval(c_low, basis, a_next)[1][:] * z_chain.p[z_i, 1]
    V_high_w = w_δt * funeval(c_high, basis, a_next)[1][:] * z_chain.p[z_i, 2]

    #! Maximize expected utlity.
    return (u_ind(c; par).+β.*(V_low_w.+V_high_w))[1]
end


riky_return(K_risky, R) = Rs 

"""
    The function to be optimized without gradient.
"""
function reward_risky(K_safe, K_risky, a, z_i, par::Pars)
    @unpack_Pars par
    # We optimize by choosing two types of investment: Housing and Capital.
    # Consumption is then simply the residual variable.

    # a includes already undepreciated capital and housing
    z = z_chain.state_values[z_i]
    
    c      = a + w * z - (P_h - P_l) * h - K

    a_next = R .* K .+ Rs .* K_risky

    V_low_w = w_Rt * funeval(c_low, basis, a_next)[1][:] * z_chain.p[z_i, 1]
    V_high_w = w_Rt * funeval(c_high, basis, a_next)[1][:] * z_chain.p[z_i, 2]

    #! Maximize expected utlity.
    return (u_ind(c; par).+β.*(V_low_w.+V_high_w))[1]
end

function hello_Irene()
    println("Hello Irene")
end