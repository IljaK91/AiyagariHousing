#c_low_init = [-304.66965824049703, -177.35499303773508, -84.58338346687913, -41.340961862539466, -22.845372703036276, -13.447504500542841, -8.093914878313798, -4.783424083486419, -2.6071557837281145, -1.1076307310363753, -0.034089767880941674, 0.7580566987710144, 1.3580565841183616, 1.6671626065072889, 1.8054585850430716]

#c_high_init = [-53.134057382850045, -37.75789860395679, -23.14944720411788, -12.758361189165873, -7.540498261124606, -4.241911969640917, -2.2119501908602426, -0.7750781912859815, 0.21922965572759567, 0.9655139711592133, 1.5224815691066174, 1.9572369894920274, 2.297862827970178, 2.4804819820713995, 2.562861344205294]

@with_kw mutable struct Pars
    β = 0.919 # Discount Rate
    δ_K = 0.04 # Capital Depreciation Rate
    r = 1 / β - δ_K - 0.02 - 1 # Net Interest Rate on Capital (0.02 for saving motive)
    R = 1 + r # Gross Interest Rate on Capital
    P_h = 1 # House price, normalized
    P_b = 1 / R # Bond price
    premium = 0.3
    P_l = (R - 1) + premium # Rental Rate for Housing, plus premium
    σ = 3.911 # CRRA coefficient
    θ = 0.8590 # Consumption Share, taken from Jensen Krueger Mitman
    w = 1.0 # wage rate

    #! Minimum buying requirement
    h_bar = 1.0
    
    #! Depreciation shocks
    δ_min = -0.0082 #minimum value of depreciation shock for housing
    δ_max = 1.0
    k = 0.7304
    σ_δ = 0.0077 #standard distribution
    N = 8 # Number of nodes

    #! This is a different formulation that distributes weight more equally
    #dist = truncated(GeneralizedExtremeValue(δ_min, σ_δ, k), δ_min, δ_max)

    dist = GeneralizedExtremeValue(δ_min, σ_δ, k)
    q0 = 0.001
    qN = cdf(dist, 1)
    δs = qnwdist(dist, N, q0, qN)[1]
    w_δ = qnwdist(dist, N, q0, qN)[2]
    w_δt = qnwdist(dist, N, q0, qN)[2]'

    #! States of wealth
    a_min = 0.1 #smin
    a_max = 20.0 #smax
    n_a = 10 # order of nodes
    type = :spli
    basis = fundefn(type, n_a, a_min, a_max) # define basis
    nodes_a = funnode(basis)[1] # Basis nodes / Collocation nodes
    phi = funbase(basis)
    c_low = phi \ log.(nodes_a) # Collocation coefficients low income state, first guess
    c_high = phi \ log.(nodes_a) # Collocation coefficients high income state, first guess

    #! States of income
    z_low = 0.1
    z_high = 1.0
    z_chain = MarkovChain([0.9 0.1; 0.1 0.9], [z_low; z_high])
    #evaluate function with funeval(c, basis, [x])[1][1]

    #! Housing is bought at the beginning or end of the period
    housing_beginning = :yes
end