function calc_DO_concentration(x, C0, B0, N0, E1, E2)
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

function DO_concentration(max_dist = 50, E1, E2)
    conc = zeros(max_dist)

    kc = 0.35
    kn = 0.25
    U = 6

    C0 = (750 + 50)/110
    B0 = (500 + 500)/110
    N0 = (500 + 350)/110
    
    for x = 0:(max_dist - 1)
        if x < 15
            conc[x + 1] = calc_DO_concentration(x, C0, B0, N0, E1, E2)
        else
            C15 = (110*calc_D0_concentration(x, C0, B0, N0, E1, E2) + 75)/125
            B15 = B0*exp(-kc*(15/U)) + 5.4
            N15 = N0*exp(-kn*(15/U)) + 4.2

            conc[x + 1] = calc_D0_concentration(x - 15, C15, B15, N15, E1, E2)
        end
    end
end

