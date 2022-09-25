using Plots

function calc_DO_concentration(x, C0, B0, N0)
    ka = 0.55
    kc = 0.35
    kn = 0.25
    Cs = 10
    U = 6

    α1 = exp(-ka*x / U)
    α2 = kc * (exp(-kc*x / U) - exp(-ka*x / U)) / (ka - kc)
    α3 = kn * (exp(-kn*x / U) - exp(-ka*x / U)) / (ka - kn)

    conc = Cs*(1-α1) + C0*α1 - B0*α2 -N0*α3

    return conc
end

function DO_concentration(E1 = 0, E2 = 0, max_dist = 50)
    conc = zeros(max_dist + 1)

    kc = 0.35
    kn = 0.25
    U = 6

    C0 = (750 + 50)/110
    B0 = (500 + (1 - E1)*500)/110
    N0 = (500 + (1 - E1)*350)/110

    C15 = (110*calc_DO_concentration(15, C0, B0, N0) + 75)/125
    B15 = (110*B0*exp(-kc*(15/U)) + (1 - E2)*675)/125
    N15 = (110*N0*exp(-kn*(15/U)) + (1 - E2)*525)/125
    
    for x = 0:(max_dist)
        if x < 15
            conc[x + 1] = calc_DO_concentration(x, C0, B0, N0)
        else
            conc[x + 1] = calc_DO_concentration(x - 15, C15, B15, N15)
        end
    end

    return conc
end

conc = DO_concentration(0, 0, 50)
plot(0:50, conc, yaxis=[3, 7.5], label=false, xlabel="Distance (km)", ylabel="DO concentration (mg/L)", linewidth=2)
vline!([0, 15], label="Waste discharge", linewidth=2)

# findall(x -> x >= 6.0, conc)[findfirst(findall(x -> x >= 6.0, conc) .>= 16)] -1

test_E2s = 0:0.01:1.0
min_DO_conc = zeros(length(test_E2s))
i = 1
for E2 in test_E2s
    # conc = DO_concentration(0, E2, 50)
    conc = DO_concentration(E2, E2, 50)
    min_DO_conc[i] = minimum(conc[16:51])
    i += 1
end
plot(test_E2s, min_DO_conc, label=false, xlabel="Treatment level (% removal)", ylabel="Minimum DO concentration (mg/L)", linewidth=2)
hline!([4.0], label="Target", linewidth=2)