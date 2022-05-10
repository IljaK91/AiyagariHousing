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
lom_housing(h, δ_shock) = max((1 - δ_shock) * h, 0)

lom_assets(b, h, δ_shock, y; par) = b + lom_housing(h, δ_shock) + y - par.B

function EV_of_s(h, s_guess; par)
    bonds = s_guess - h * (1 - par.r)
    a_next = lom_assets.(bonds, h, par.shock_states[:, 2], par.shock_states[:, 1]; par)
    V_next_big = par.interp_linear.(a_next)
    EV_next = sum(par.p_states_grid .* V_next_big)

    return EV_next
end

function find_EV(s_guess; par)
    f(h) = -EV_of_s(h, s_guess; par)
    sol = optimize(f, 0.0, s_guess / (1 - par.r) - 1e-6)
    return sol.minimizer, -sol.minimum
end

"""
    Saving decision.
"""
function find_EV_h_bar(s_guess; par)
    f(h) = -EV_of_s(h, s_guess; par)
    @set! par.h_bar_mort = max(par.h_bar, par.h_min_mort)
    #@show par.h_bar_mort
    #@show s_guess / (1 - par.r)
    if par.h_bar_mort >= s_guess / (1 - par.r) - 1e-6
        return 0, EV_of_s(0, s_guess; par)
    elseif par.h_bar_mort < s_guess / (1 - par.r) - 1e-6
        sol = optimize(f, par.h_bar_mort, s_guess / (1 - par.r) - 1e-6)
        return sol.minimizer, -sol.minimum
    else
        error("I should never be here")
    end
end

"""
    This one will be depreciated in favor of Val_of_c_h_bar with h_bar = 0.
"""
function Val_of_c(c, a; par)
    s = a - c
    h_sol, EV_next = find_EV(s; par)
    u_ind(c; par) + par.β * EV_next
end

function Val_of_c_h_bar(c, a; par)
    s = a - c
    h_sol, EV_next = find_EV_h_bar(s; par)
    u_ind(c; par) + par.β * EV_next
end

function Val_of_c_h_zero(c, a; par)
    s = a - c
    u_ind(c; par) + par.β * EV_of_s(0, s; par)
end

function h_constr(c, a; par)
    Val_h_bar = Val_of_c_h_bar(c, a; par)
    Val_h_zero = Val_of_c_h_zero(c, a; par)
    s = a - c
    if Val_h_bar > Val_h_zero
        find_EV_h_bar(s; par)[1]
    elseif Val_h_bar <= Val_h_zero
        0.0
    else
        @show Val_h_zero
        @show Val_h_bar
        error("I should never be here")
    end
end

function find_c(a; par)
    f(c) = -Val_of_c(c, a; par)
    sol = optimize(f, 1e-6, a - 1e-6)
    return sol.minimizer, -sol.minimum
end

function find_c_h_bar(a; par)
    f(c) = -Val_of_c_h_bar(c, a; par)
    sol = optimize(f, 1e-6, a - 1e-6)
    return sol.minimizer, -sol.minimum
end

"""
    Two steps in the solution algorithm:

    What is the optimal decision if we constrain ourselves to no housing?
    What is the optimal decision if we consume at least h_bar housing?

    Need to extend by assuming that borrowing is possible.
"""
function find_c_constr(a; par)
    if par.B > 0
        c_h_bar, V_h_bar = find_c_h_bar(a; par)
        return c_h_bar, V_h_bar
    elseif par.B == 0.0
        c_h_zero, V_h_zero = find_c_h_zero(a; par)
        c_h_bar, V_h_bar = find_c_h_bar(a; par)
    
        if V_h_zero > V_h_bar
            return c_h_zero, V_h_zero
        elseif V_h_zero <= V_h_bar
            return c_h_bar, V_h_bar
        else
            error("This cannot happen!")
        end
    end
end

"""
    Two steps in the solution algorithm:

    What is the optimal decision if we constrain ourselves to no housing?
    What is the optimal decision if we consume at least h_bar housing?

    Need to extend by assuming that borrowing is possible.
"""
function find_c_constr_mortgages(a, B; par)
    a_with_mort = a + B
    sol_c[i], V_a_sol[i] = find_c_constr(a_with_mort; par)

end

"""
    Optimal borrowing decision.
"""
function find_optimal_mortage(a; par, max_B = effective_B_max(a; par))
    min_B = 0
    #@set! h_min_mort = h_min(B; par)
    f(B) = begin
        a_new = a + B
        @set! par.B = B
        @set! par.h_min_mort = h_min(B; par)
        c_h_bar, V_h_bar = find_c_h_bar(a_new; par)

        return -V_h_bar
    end

    sol = optimize(f, min_B, max_B * 0.99)

    return sol.minimizer, -sol.minimum
end


"""
    Maximal possible borrowing.
"""
B_max(a; par) = (1 - par.Eδ) * a / par.R

"""
    Maximal possible borrowing.
"""
function effective_B_max(a; par)
    max_B = B_max(a; par)
    if (a + max_B) / (1 - par.r) < par.h_bar
        0.0
    elseif (a + max_B) / (1 - par.r) > par.h_bar
        max_B
    end
end 

"""
    Minimum necessary amount of housing that needs to be bought to back up mortgage.
"""
h_min(B; par) = par.R * B / (1 - par.Eδ)


function find_c_h_zero(a; par)
    f(c) = -Val_of_c_h_zero(c, a; par)
    sol = optimize(f, 1e-6, a - 1e-6)
    return sol.minimizer, -sol.minimum
end

function create_grid(y_min, init_step, scale, max)
    grid = Vector()
    push!(grid, y_min)
    step = copy(init_step)
    push!(grid, y_min + step)
    while maximum(grid) < max
        step = copy(scale * step)
        push!(grid, y_min + step)
    end
    return grid
end

function solve_model(par; max_B = nothing)
    V_a_sol = zeros(length(par.a_states))
    sol_s = zeros(length(par.a_states))
    sol_c = zeros(length(par.a_states))

    ##! First guess just uses the log 
    V_guess = u_ind.(par.a_states; par)
    tol = 1e-8
    dif = 10
    counter = 0
    @set! par.interp_linear = LinearInterpolation(par.a_states, V_guess, extrapolation_bc = Line())

    while dif > tol
        counter += 1
        @show counter
        interp_linear = LinearInterpolation(par.a_states, V_guess, extrapolation_bc = Line())
        for i in eachindex(V_guess)
            @show i
            #! Check first whether household takes up a mortgage
            if isnothing(max_B)
                sol_B, sol_V_h_bar = find_optimal_mortage(par.a_states[i]; par)
            else
                sol_B, sol_V_h_bar = find_optimal_mortage(par.a_states[i]; par, max_B = max_B)
            end
            #! This is all a bit inefficient, as I solve the same thing over and over. Especially the lower part, as the upper part may already 
            @set! par.B = sol_B
            #! Then check what the next solution is
            if par.h_bar > 0
                sol_c[i], V_a_sol[i] = find_c_constr(par.a_states[i] + par.B; par)
            elseif par.h_bar == 0
                sol_c[i], V_a_sol[i] = find_c(par.a_states[i] + par.B; par)
            end
        end
        dif = sum((V_a_sol - V_guess) .^ 2)
        @show dif
        if counter >= 100
            break
        elseif dif > tol
            V_guess = copy(V_a_sol)
            @set! par.interp_linear = LinearInterpolation(par.a_states, V_guess, extrapolation_bc = Line())
        else
            break
        end
    end

    return V_a_sol
end
