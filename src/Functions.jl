"""
    u(c, h; par)

Direct CRRA utility function with a share θ for the normal consumption good and 1-θ for housing. σ is the coefficient of relative risk-aversion.
# Arguments:
- `c`: Consumption.
- `h`: Housing.
- `par`: Struct of Parameters.
"""
u(c, h; par) = (c^par.θ * h^(1 - par.θ))^(1 - par.σ) / (1 - par.σ)

"""
    Indirect utility function as in Jeske, Krueger, Mitman
    where c = c_tilde + P_l*h is equal to total expenditure.
"""
function u_ind(c; par)
    (Psi(par) * c^(1 - par.σ) - 1) / (1 - par.σ)
end

"""
    Psi(par::Pars)

Auxiliary function for the indirect utility function.

# Arguments:
- `par::Pars`: Struct of Parameters
"""
function Psi(par::Pars)
    @unpack_Pars par
    (θ^θ * (1 - θ)^(1 - θ) * r^(θ - 1))^(1 - σ)
end

"""
    Law of motion for housing investments
"""
lom_housing(h, δ_shock) = max((1 - δ_shock) * h, 0)

lom_assets(b, g, δ_shock, y) = b + lom_housing(g, δ_shock) + y


# """
#     The function to be optimized.
# """
# function myfunc(x::Vector, grad::Vector, a, z_i, par::Pars, h_next, a_next, V_low_w, V_high_w)
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

#     mul!(V_low_w, w_δt, funeval(c_low, basis, a_next)[1][:], z_chain.p[z_i, 1], 0)
#     mul!(V_high_w, w_δt, funeval(c_high, basis, a_next)[1][:], z_chain.p[z_i, 2], 0)

#     #! Maximize expected utlity.
#     return -(u_ind(c; par).+β.*(V_low_w.+V_high_w))[1]
# end

# """
#     Budget constraint.
# """
# function myconstraint(x::Vector, grad::Vector, a, z_i, par::Pars)
#     @unpack_Pars par
#     z = z_chain.state_values[z_i]
#     h = x[1]
#     K = x[2]
#     #! Basically: Investment needs to be less than total assets
#     return (P_h - P_l) * h + K - a - w * z # ≤ 0
# end

# function vmax(a, z_i; par, h_next, a_next, V_low_w, V_high_w, x0=[0.04, 0.04])
#     # Two variables, choose gradient-free algorithm
#     opt = Opt(:LN_COBYLA, 2)
#     opt.lower_bounds = [0.0, 0.0]
#     opt.xtol_rel = 1e-4
#     # define objective
#     myfinalfunc(x, g) = myfunc(x, g, a, z_i, par, h_next, a_next, V_low_w, V_high_w)
#     opt.min_objective = myfinalfunc
#     #define constraint
#     inequality_constraint!(opt, (x, g) -> myconstraint(x, g, a, z_i, par), 1e-8)
#     (minf, minx, ret) = optimize(opt, x0)

#     return minf, minx, ret
# end


# function vmax(a, z_i; par, h_next, a_next, V_low_w, V_high_w, x0=[0.04, 0.04])
#     # Two variables, choose gradient-free algorithm
#     opt = Opt(:LN_COBYLA, 2)
#     opt.lower_bounds = [0.0, 0.0]
#     opt.xtol_rel = 1e-4
#     # define objective
#     myfinalfunc(x, g) = myfunc(x, g, a, z_i, par, h_next, a_next, V_low_w, V_high_w)
#     opt.min_objective = myfinalfunc
#     #define constraint
#     inequality_constraint!(opt, (x, g) -> myconstraint(x, g, a, z_i, par), 1e-8)
#     (minf, minx, ret) = optimize(opt, x0)

#     return minf, minx, ret
# end

#! Additional constraint to add: Minimum buying requirement.

"""
    myfunc2(h, K, a, z_i, par::Pars)

DOCSTRING

# Arguments:
- `h`: Housing investment.
- `K`: Capital investment.
- `a`: Assets.
- `z_i`: Index of the labor productivity shock.
- `par`: Struct of Parameters.
"""
function myfunc2(h, K, a, z_i, par::Pars)
    #@unpack_Pars par
    # We optimize by choosing two types of investment: Housing and Capital.
    # Consumption is then simply the residual variable.

    # a includes already undepreciated capital and housing

    z = par.z_chain.state_values[z_i]
    if par.housing_beginning == :yes
        c = a + par.w * z - (par.P_h - par.P_l) * h - K
        h_next = lom_housing.(h, par.δs)
    elseif par.housing_beginning == :no
        c = a + par.w * z - par.P_h * h - K
        h_next = lom_housing.(h, par.δs) .+ h .* par.P_l
    end

    a_next = par.R .* K .+ h_next

    V_low_w = par.w_δt * funeval(par.c_low, par.basis, a_next)[1][:] * par.z_chain.p[z_i, 1]
    V_high_w = par.w_δt * funeval(par.c_high, par.basis, a_next)[1][:] * par.z_chain.p[z_i, 2]

    #! Maximize expected utlity.
    return (u_ind(c; par).+par.β.*(V_low_w.+V_high_w))[1]
end

"""
    calc_step(a_i, par::Pars)

Calculate the step size between a given node with index a_i and the next/previous node.

# Arguments:
- `a_i`: Index for an asset node.
- `par`: Struct of Parameters.
"""
function calc_step(a_i, par::Pars)
    if a_i == 1
        step_up = par.nodes_a[a_i+1] - par.nodes_a[a_i]
        step_down = 0
    elseif a_i == par.n_a
        step_up = 0
        step_down = par.nodes_a[a_i] - par.nodes_a[a_i-1]
    else
        step_up = par.nodes_a[a_i+1] - par.nodes_a[a_i]
        step_down = par.nodes_a[a_i] - par.nodes_a[a_i-1]
    end
    return step_up, step_down
end

"""
    delta_low(a_next, K, h, step, par::Pars)

Finds the minimal delta that would still lead the next realization of assets to be attributed to the node a_next.

# Arguments:
- `a_next`: Amount of assets in the next period.
- `K`: Capital investment.
- `h`: Housing investment.
- `step`: Step size between a_next and the next larger node.
- `par`: Struct of Parameters.
"""
function delta_low(a_next, K, h, step, par::Pars)
    if iszero(h)
        -Inf
    else
        1 - (a_next + step / 2 - par.R * K) / h #Irene: you do not have P_{t+1}^h (ie P_h) on the denominator? Or is it because you already assume that the price housing P_h is 1.0?
    end
end

"""
    delta_high(a_next, K, h, step, par::Pars)

Finds the maximum delta that would still lead the next realization of assets to be attributed to the node a_next.

# Arguments:
- `a_next`: Amount of assets in the next period.
- `K`: Capital investment.
- `h`: Housing investment.
- `step`: Step size between a_next and the next lower node.
- `par`: Struct of Parameters.
"""
function delta_high(a_next, K, h, step, par::Pars)
    if iszero(h)
        -Inf
    else
        1 - (a_next - step / 2 - par.R * K) / h
    end
end

"""
    Prob_a_next(K, h, a_i_next, par::Pars)

Calculates the probability to transition to a gridpoint with index a_i_next

# Arguments:
- `K`: Capital investment.
- `h`: Housing investment.
- `a_i_next`: Index for the next asset state.
- `par`: Struct of Parameters
"""
function Prob_a_next(K, h, a_i_next, par::Pars)
    @unpack_Pars par
    step_up, step_down = calc_step(a_i_next, par)
    a_next = par.nodes_a[a_i_next]
    if h < 1e-6
        if a_i_next == 1
            if R * K <= a_next + step_up / 2
                return 1.0
            else
                return 0.0
            end
        elseif a_i_next == length(nodes_a)
            if a_next - step_down / 2 <= R * K
                return 1.0
            else
                return 0.0
            end
        else
            if a_next - step_down / 2 <= R * K <= a_next + step_up / 2
                return 1.0
            else
                return 0.0
            end
        end

    elseif h > 1e-6
        δ_low = delta_low(a_next, K, h, step_up, par)
        δ_high = delta_high(a_next, K, h, step_down, par)
        @assert δ_low < δ_high
        if a_i_next == 1
            if δ_high > 1 # need to check because we effectively lump all probabilities with δ > 1 to δ = 1
                return 0
            else
                return 1 - cdf(GeneralizedExtremeValue(δ_min, σ_δ, k), δ_high) # delta needs to be at least equal to delta_low
            end
            #! Seems ok
        elseif a_i_next == length(nodes_a)
            cdf(GeneralizedExtremeValue(δ_min, σ_δ, k), δ_high)
        else
            if δ_high > 1 && δ_low < 1
                return 1 - cdf(GeneralizedExtremeValue(δ_min, σ_δ, k), δ_low) # in this case, all values larger than δ_low count
            elseif δ_low > 1
                # delta cannot take values larger than one
                return 0.0
            else
                # normal case, take the area between
                return cdf(GeneralizedExtremeValue(δ_min, σ_δ, k), δ_high) - cdf(GeneralizedExtremeValue(δ_min, σ_δ, k), δ_low)
            end
        end
    end
end

"""
    setup_Q(sol_K, sol_h, par::Pars)

Populates the transition matrix Q based on the policy functions for capital and housing.

# Arguments:
- `sol_K`: Policy function for capital investment.
- `sol_h`: Policy function for housing investment.
- `par`: Struct of Parameters.
"""
function setup_Q(sol_K, sol_h, par::Pars)
    @unpack_Pars par
    Q = zeros(n_a * length(par.z_chain.state_values), n_a * length(par.z_chain.state_values))
    for a_i in eachindex(par.nodes_a) # action
        for z_i in eachindex(par.z_chain.state_values) # current state
            #z = par.z_chain.state_values[z_i] # z state
            K = sol_K[a_i, z_i] # action 1, capital
            h = sol_h[a_i, z_i] # action 2, housing shares
            for next_a_i in eachindex(par.nodes_a)
                for next_z_i in eachindex(eachindex(par.z_chain.state_values))
                    Q[a_i+(z_i-1)*n_a, next_a_i+(next_z_i-1)*n_a] = z_chain.p[z_i, next_z_i] * Prob_a_next(K, h, next_a_i, par)
                end
            end
        end
    end
    #! Need to check that without the normalization the rows of Q sum up to 1
    for j in 1:size(Q, 2)
        for i in 1:size(Q, 1)
            Q[i, :] = Q[i, :] ./ sum(Q[i, :])
        end
    end
    return Q
end

# function myfunc2(h, K, a, z_i, par::Pars)
#     @unpack_Pars par
#     # We optimize by choosing two types of investment: Housing and Capital.
#     # Consumption is then simply the residual variable.

#     # a includes already undepreciated capital and housing

#     z = z_chain.state_values[z_i]
#     if housing_beginning == :yes
#         c = a + w * z - (P_h - P_l) * h - K
#         h_next = lom_housing.(h, δs)
#     elseif housing_beginning == :no
#         c = a + w * z - P_h * h - K
#         h_next = lom_housing.(h, δs) .+ h .* P_l
#     end

#     a_next = R .* K .+ h_next

#     V_low_w = w_δt * funeval(c_low, basis, a_next)[1][:] * z_chain.p[z_i, 1]
#     V_high_w = w_δt * funeval(c_high, basis, a_next)[1][:] * z_chain.p[z_i, 2]

#     #! Maximize expected utlity.
#     return (u_ind(c; par).+β.*(V_low_w.+V_high_w))[1]
# end
# """
#     reward_risky(K_safe, K_risky, a, z_i, par::Pars)

# DOCSTRING

# # Arguments:
# - `K_safe`: DESCRIPTION
# - `K_risky`: DESCRIPTION
# - `a`: DESCRIPTION
# - `z_i`: DESCRIPTION
# - `par`: DESCRIPTION
# """
# function reward_risky(K_safe, K_risky, a, z_i, par::Pars)
#     @unpack_Pars par
#     # We optimize by choosing two types of investment: Housing and Capital.
#     # Consumption is then simply the residual variable.

#     # a includes already undepreciated capital and housing
#     z = z_chain.state_values[z_i]

#     c = a .+ w .* z .- K_risky .- K_safe

#     a_next = R .* K_safe .+ Rs .* K_risky

#     V_low_w = w_Rt * funeval(c_low, basis, a_next)[1][:] * z_chain.p[z_i, 1]
#     V_high_w = w_Rt * funeval(c_high, basis, a_next)[1][:] * z_chain.p[z_i, 2]

#     #! Maximize expected utlity.
#     return (u_ind(c; par).+β.*(V_low_w.+V_high_w))[1]
# end

# """
#     u_static(c; par::Pars)

# DOCSTRING
# """
# u_static(c; par::Pars) = (c)^(1 - par.σ) / (1 - par.σ)



# """
#     reward_risky_static(K_risky, a, par::Pars)

#     Try to understand CRRA utility a little bit better. Optimize the allocation between risky and safe asset in a static setting.

# # Arguments:
# - `K_risky`: DESCRIPTION
# - `a`: DESCRIPTION
# - `par`: DESCRIPTION
# """
# function reward_risky_static(K_risky, a, par::Pars)
#     @unpack_Pars par

#     K_safe = a .- K_risky

#     a_next = R .* K_safe .+ Rs .* K_risky

#     #! Maximize expected utlity.
#     return (w_Rt*u_static.(a_next; par))[1]
# end

"""
    solve_model(par; tol=1e-8, maxit=500)

Performs value function iteration, solving for the optimal policy functions for housing and capital investment.

# Arguments:
- `par`: Struct of Parameters.
- `tol`: Tolerance for convergence of the value function. Default value 1e-8.
- `maxit`: Maximum number of iterations. Default value 500.
"""

function solve_model(par; tol=1e-8, maxit=500)
    i = 1
    diff = 1000
    V_a_guess = zeros(length(par.nodes_a), length(par.z_chain.state_values))
    V_a_sol = zeros(length(par.nodes_a), length(par.z_chain.state_values))
    sol_K = zeros(length(par.nodes_a), length(par.z_chain.state_values))
    sol_h = zeros(length(par.nodes_a), length(par.z_chain.state_values))
    sol_c = zeros(length(par.nodes_a), length(par.z_chain.state_values))
    sol_y = zeros(length(par.nodes_a), length(par.z_chain.state_values))

    while diff > tol && i < maxit
        for z_i in eachindex(par.z_chain.state_values)
            for a_i in eachindex(par.nodes_a)

                a = par.nodes_a[a_i]
                z = par.z_chain.state_values[z_i]

                #model = Model(NLopt.Optimizer)
                #set_optimizer_attribute(model, "algorithm", :LD_MMA)
                model = Model(Ipopt.Optimizer)
                set_silent(model)
                @variable(model, a + par.w * z >= h >= 0)
                @variable(model, a + par.w * z >= K >= 0)

                f(h, K) = myfunc2(h, K, a, z_i, par)
                register(model, :f, 2, f, autodiff=true)

                @NLobjective(model, Max, f(h, K))
                if par.housing_beginning == :yes
                    @NLconstraint(model, (par.P_h - par.P_l) * h + K - a - par.w * z <= 0)
                elseif par.housing_beginning == :no
                    @NLconstraint(model, par.P_h * h + K - a - par.w * z <= 0)
                end

                if i > 1
                    set_start_value(h, sol_h[a_i, z_i])
                    set_start_value(K, sol_K[a_i, z_i])
                end

                JuMP.optimize!(model)
                if termination_status(model) != OPTIMAL && termination_status(model) != LOCALLY_SOLVED && termination_status(model) != ALMOST_LOCALLY_SOLVED
                    @show termination_status(model)
                    error("No solution was found!")
                end
                V_a_sol[a_i, z_i] = objective_value(model)
                sol_h[a_i, z_i] = value(h)
                sol_K[a_i, z_i] = value(K)
                sol_y[a_i, z_i] = a + par.w * z
                if par.housing_beginning == :yes
                    sol_c[a_i, z_i] = a + par.w * z - (par.P_h - par.P_l) * sol_h[a_i, z_i] - sol_K[a_i, z_i]
                elseif par.housing_beginning == :no
                    sol_c[a_i, z_i] = a + par.w * z - par.P_h * sol_h[a_i, z_i] - sol_K[a_i, z_i]
                end
            end
        end
        diff = copy(sum((V_a_guess .- V_a_sol) .^ 2))
        @printf "Iteration: %i. Residual: %.3e\n" i diff
        if diff < tol
            break
        else
            @set! par.c_low = par.phi \ V_a_sol[:, 1]
            @set! par.c_high = par.phi \ V_a_sol[:, 2]
            V_a_guess = copy(V_a_sol)
            i += 1
        end
    end
    return V_a_sol, sol_K, sol_h, sol_c, sol_y, par
end

"""
    stst_distr(sol_K, sol_h, par::Pars; tol = 1.0e-8)

Derives the steady state distribution through iteration using the transition matrix.

# Arguments:
- `sol_K`: Policy function for capital investment.
- `sol_h`: Policy function for housing investment.
- `par`: Struct of Parameters.
- `tol`: Tolerance for convergence of the distribution. Default value 1e-8.
"""
function stst_distr(sol_K, sol_h, par::Pars; tol=1e-8)
    Q = setup_Q(sol_K, sol_h, par)
    distr = ones(par.n_a * 2) .* 1 / (par.n_a * 2)
    dif = 1000
    while dif > tol
        new_distr = Q' * distr
        dif = sum((new_distr - distr) .^ 2)
        if dif > tol
            distr = copy(new_distr)
        elseif dif <= tol
            distr = copy(new_distr)
            break
        else
            error("I should never be here!")
        end
    end

    # Normalize distribution, inconsequential
    distr = distr ./ sum(distr)

    return distr
end

"""
    marketclearing_housing(sol_K, sol_h, par::Pars; fixed_supply::Bool = false)

Return total housing ownership minus total housing demand using policy functions derived from the solution of the model and the steady state distribution to measure aggregate demand and supply.

Takes the keyword argument fixed_supply, where true means that total housing supply is fixed to one and false means that total housing supply can adjust.

# Arguments:
- `sol_K`: Policy function for capital investment.
- `sol_h`: Policy function for housing investment.
- `par`: Struct of Parameters.
- `fixed_supply`: true means that total housing supply is fixed to one. false for elastic housing supply.
"""
function marketclearing_housing(sol_K, sol_h, par::Pars; fixed_supply::Bool=false)

    distr = stst_distr(sol_K, sol_h, par)
    # Reshape parameters
    sol_h_long = reshape(sol_h, (par.n_a * 2, 1))
    sol_c_long = reshape(sol_c, (par.n_a * 2, 1))

    total_housing_ownership = sum(sol_h_long .* distr)
    # Since preferences are cobb-douglas over housing and the consumption good, households will spend a fraction θ of total consumption expenditure on the consumption good and a fraction 1-θ on housing. To get the total amount of rented housing, need to also divide by the price of housing, which is <1.
    total_housing_demand = sum((1 - par.θ) .* sol_c_long ./ par.P_l .* distr)

    #! Give one or two residuals depending on whether the supply of housing is fixed to one or not.
    if fixed_supply
        1 - total_housing_demand, 1 - total_housing_ownership
    else
        total_housing_ownership - total_housing_demand, 0.0
    end

end