using Plots, Distributions, PDMats

function calc_mixed_concentration(V0, C0, V1, C1)
    return (V0 * C0 + V1 * C1)/(V0 + V1)
end

function calc_concentrations(x, C0, B0, N0)
    ka = 0.55
    kc = 0.35
    kn = 0.25
    Cs = 10
    U = 6

    α1 = exp(-ka*x / U)
    α2 = kc * (exp(-kc*x / U) - exp(-ka*x / U)) / (ka - kc)
    α3 = kn * (exp(-kn*x / U) - exp(-ka*x / U)) / (ka - kn)

    C = Cs*(1-α1) + C0*α1 - B0*α2 - N0*α3
    B = B0*exp(-kc*(x/U))
    N = N0*exp(-kn*(x/U))

    return C, B, N
end

function DO_simulation(E1 = 0, E2 = 0, max_dist = 50)
    C = zeros(max_dist + 1)
    B = zeros(max_dist + 1)
    N = zeros(max_dist + 1)

    C[1] = calc_mixed_concentration(1e5, 7.5, 1e4, 5)
    B[1] = calc_mixed_concentration(1e5, 5, 1e4, (1 - E1)*50)
    N[1] = calc_mixed_concentration(1e5, 5, 1e4, (1 - E1)*35)
    
    for x = 1:(max_dist)
        if x < 15
            C[x + 1], B[x + 1], N[x + 1] = calc_concentrations(x, C[1], B[1], N[1])
        elseif x == 15
            C[16] = calc_mixed_concentration(1.1e5, C[15], 1.5e4, 5)
            B[16] = calc_mixed_concentration(1.1e5, B[15], 1.5e4, (1 - E2)*45)
            N[16] = calc_mixed_concentration(1.1e5, N[15], 1.5e4, (1 - E2)*35)
        else
            C[x + 1], B[x + 1], N[x + 1] = calc_concentrations(x - 15, C[16], B[16], N[16])
        end
    end

    return C, B, N
end

function stochastic_DO_simulation1(E1 = 0, E2 = 0, max_dist = 50)
    C = zeros(max_dist + 1)
    B = zeros(max_dist + 1)
    N = zeros(max_dist + 1)

    C[1] = calc_mixed_concentration(1e5, 7.5, 1e4, 5)
    B[1] = calc_mixed_concentration(1e5, rand(Uniform(4, 7)), 1e4, (1 - E1)*50)
    N[1] = calc_mixed_concentration(1e5, rand(Uniform(3, 8)), 1e4, (1 - E1)*35)
    
    for x = 1:(max_dist)
        if x < 15
            C[x + 1], B[x + 1], N[x + 1] = calc_concentrations(x, C[1], B[1], N[1])
        elseif x == 15
            C[16] = calc_mixed_concentration(1.1e5, C[15], 1.5e4, 5)
            B[16] = calc_mixed_concentration(1.1e5, B[15], 1.5e4, (1 - E2)*45)
            N[16] = calc_mixed_concentration(1.1e5, N[15], 1.5e4, (1 - E2)*35)
        else
            C[x + 1], B[x + 1], N[x + 1] = calc_concentrations(x - 15, C[16], B[16], N[16])
        end
    end

    return C, B, N
end

min_DO_conc = zeros(1000)
for n in 1:1000
    DO_conc, CBOD_conc, NBOD_conc = stochastic_DO_simulation1(0.15, 0.15, 50)
    min_DO_conc[n] = minimum(DO_conc)
end
sum(min_DO_conc .>= 4.0)/10

function sample_correlated_uniform(n, a, b, corr_coef=0.7)
    mvnorm = MvNormal([0, 0], PDMat([1 corr_coef; corr_coef 1])) # set up a multivariate normal with each marginal variance of 1 and the right correlation
    norm_samples = rand(mvnorm, n)' # sample from the multivariate normal, the marginal distributions are a standard normal
    unif_samples = cdf.(Normal(0, 1), norm_samples) # convert samples to a uniform distribution using the pdf of a standard Normal
    samples = (unif_samples .* [a[2] - a[1] b[2] - b[1]]) .+ [a[1] b[1]]
    return samples
end

function stochastic_DO_simulation2(E1 = 0, E2 = 0, max_dist = 50)
    C = zeros(max_dist + 1)
    B = zeros(max_dist + 1)
    N = zeros(max_dist + 1)

    inflow_org_matter = sample_correlated_uniform(1, [4, 7], [3, 8])

    C[1] = calc_mixed_concentration(1e5, 7.5, 1e4, 5)
    B[1] = calc_mixed_concentration(1e5, inflow_org_matter[1], 1e4, (1 - E1)*50)
    N[1] = calc_mixed_concentration(1e5, inflow_org_matter[2], 1e4, (1 - E1)*35)
    
    for x = 1:(max_dist)
        if x < 15
            C[x + 1], B[x + 1], N[x + 1] = calc_concentrations(x, C[1], B[1], N[1])
        elseif x == 15
            C[16] = calc_mixed_concentration(1.1e5, C[15], 1.5e4, 5)
            B[16] = calc_mixed_concentration(1.1e5, B[15], 1.5e4, (1 - E2)*45)
            N[16] = calc_mixed_concentration(1.1e5, N[15], 1.5e4, (1 - E2)*35)
        else
            C[x + 1], B[x + 1], N[x + 1] = calc_concentrations(x - 15, C[16], B[16], N[16])
        end
    end

    return C, B, N
end

min_DO_conc = zeros(1000)
for n in 1:1000
    DO_conc, CBOD_conc, NBOD_conc = stochastic_DO_simulation2(0.2, 0.2, 50)
    min_DO_conc[n] = minimum(DO_conc)
end
sum(min_DO_conc .>= 4.0)/10