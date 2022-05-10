"""
    The parameters of our model.
"""
@with_kw struct Pars
    β = 0.919 # Discount Rate
    # Probably the interest rate is too high, it is too attractive to accumulate savings endlessly.
    R = 1 / β - 0.06 # Interest Rate
    P_h = 1 # House price
    P_b = 1 / R # Bond price
    r = (R - 1) + 0.006 # Rental Rate
    σ = 3.911 # CRRA coefficient
    θ = 0.8590 # Consumption Share, taken from Jensen Krueger Mitman

    h_bar = 0.5
    h_bar_mort = 0.0 # minimal buying requirement due to mortgages

    y_min = 0.1 # minimal income shock
    y_max = 2 # maximal income shock
    y_steps = 20 # number of steps for income shock
    y_states = collect(range(y_min, y_max, y_steps)) # grid income shock

    p_y_states = 1 / y_steps .* ones(length(y_states)) # equal probability income shocks

    δ_min = -0.5 # minimal depreciation shock
    δ_max = 1 # maximal depreciation shock
    δ_steps = 20 # number of steps for depreciation shock
    δ_states = collect(range(δ_min, δ_max, δ_steps)) # grid depreciation shock

    p_δ_states = pdf.(Normal(0, 0.3), δ_states) ./ sum(pdf.(Normal(0, 0.3), δ_states)) # probabilities depreciation shock, borrowed from normal distribution
    Eδ = sum(δ_states.*p_δ_states)

    #! Mortgages
    λ = 0.8 # financial friction
    h_min_mort = 0. # how much housing needs to be bought at least to back up the mortgage.
    B = 0.0 # Borrowing

    shock_states = gridmake(y_states, δ_states) # grid both shocks
    p_states_intermediate = gridmake(p_y_states, p_δ_states) # probability of all states
    p_states_grid = p_states_intermediate[:, 1] .* p_states_intermediate[:, 2] # grid probabilities

    a_max = 10 # maximal wealth
    #a_states = create_grid(0.2, 0.001, 1.2, a_max) # minimal to maximal wealth
    a_states = collect(y_min:0.05:a_max) # minimal to maximal wealth
    #log.(y_min:0.05:a_max) 
    μ_a_states = zeros(length(a_states)) # mass distribution over states, equilibrium outcome.

    interp_linear = LinearInterpolation(a_states, zeros(length(a_states)), extrapolation_bc = Line())
end