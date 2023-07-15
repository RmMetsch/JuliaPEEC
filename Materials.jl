using PlotlyJS, BenchmarkTools

## Tape dimensions, all values are crossectional surface area's in [m^2]
TapeWidth = 0.01
CopperA   = 50e-6 * TapeWidth
HastTA    = 50e-6 * TapeWidth
SilverA   = 4e-6  * TapeWidth
RebcoA    = 2e-6  * TapeWidth


## Rebco mat properties
function Rebco_Ic(T)
    # Scaling law input parameters
    # General parameters
    w = 0.01
    B = 0
    angle = 0
    T_c0 = 93  # [K]

    # Parameters for ab-plane
    alpha_ab = 3.1e13  # [AT/m^2]
    gamma_ab = 1.8
    a = 0.18
    n_ab = 0.72
    n1 = 1.7
    n2 = 6.9
    Bi0_ab = 250  # [T]
    p_ab = 0.88
    q_ab = 1.6

    # Parameters for c-plane
    alpha_c = 8.1e12  # [AT/m^2]
    gamma_c = 2.1
    Bi0_c = 140  # [T]
    n_c = 0.4
    p_c = 0.64
    q_c = 2.3

    # Parameters for anisotropy
    g0 = 0.045
    g1 = 0.23
    g2 = 0.014
    g3 = 0.1
    nu = 2

    # Field relations
    # prevent zero field at which an anomalous behavior occurs
    B = max.(B, 0.1)

    # dimensionless temperature
    T = T ./ T_c0

    Bi_ab = Bi0_ab * ((1 .- complex(T) .^ n1) .^ n2 + a * (1 .- complex(T) .^ n_ab))
    Bi_c = Bi0_c * (1 .- complex(T) .^ n_c)

    b_ab = B ./ Bi_ab
    b_c = B ./ Bi_c

    Jc_ab = ifelse.((T .< 1.0) .& (real(b_ab) .< 1.0), (alpha_ab ./ B) .* (b_ab .^ p_ab) .* (1 .- b_ab) .^ q_ab .* ((1 .- complex(T) .^ n1) .^ n2 + a * (1 .- complex(T) .^ n_ab)) .^ gamma_ab, 0)
    Jc_c = ifelse.((T .< 1.0) .& (real(b_c) .< 1), (alpha_c ./ B) .* (b_c .^ p_c) .* (1 .- b_c) .^ q_c .* (1 .- complex(T) .^ n_c) .^ gamma_c, 0)

    g = g0 .+ g1 .* exp.(-g2 .* B .* exp.(g3 .* T))

    Jc_ab = real(Jc_ab)
    Jc_c  = real(Jc_c)
        

    Jc = min.(Jc_ab, Jc_c) .+ max.(0, Jc_ab .- Jc_c) ./ (1 .+ ((angle .- pi / 2) ./ g) .^ nu)

    # Convert units
    # thickness of superconducting layer
    dsup = 2e-6  # from J. Fleiter report 
    return Jc .* w .* dsup   # crittical current [A]

end

function Rebco_C(T)
    # From j.Nughteren
    # output in J/m^3K
    a = -7.567485538209158e-10
    b = 6.351452642016898e-7
    c = -1.947975786547597e-4
    d = 0.023616673974415
    e = 0.239331954284042
    f = -1.096191721280114

    return 6.3e3 .* (a .* T .^ 5 .+ b .* T .^ 4 .+ c .* T .^ 3 .+ d .* T .^ 2 .+ e .* T .+ f)
end

function Rebco_k(T)
    # from bergen
    # takes a vector of Temperatures in K
    # outputs a vector of specific heat conductivities in W/(m*K)
    a = [-1.266103106492942e-05, -1.316704893540544e-09]
    b = [0.002670105477219, 2.636151594117440e-06]
    c = [-0.197035302542769, -0.001601689632073]
    d = [4.933962659604384, 0.428760641312538]
    e = [29.651536939501670, -53.306643333287260]
    f = [66.578192505447330, 3.682252343338599e+03]

    n = ifelse.(T .<= 70, 1, 2)

    thermal_conductivity_ybco = a[n] .* T .^ 5 .+ b[n] .* T .^ 4 .+ c[n] .* T .^ 3 .+ d[n] .* T .^ 2 .+ e[n] .* T .+ f[n]

    return thermal_conductivity_ybco
end

## Copper matieral properties
function Copper_C(T)
    @. 8.94e3 .* (
        (T .>= 1.00) .& (T .< 17.50) .* (2.24642893e-05 .* T.^4 + 2.84703334e-04 .* T.^3 + 3.44121866e-03 .* T.^2 + 1.04457033e-03 .* T + 8.16805501e-03) .+
        (T .>= 17.50) .& (T .< 62.00) .* (3.01020641e-07 .* T.^4 - 1.04836396e-03 .* T.^3 + 1.54053918e-01 .* T.^2 - 3.76716858e+00 .* T + 2.90597210e+01) .+
        (T .>= 62.00) .& (T .< 300.00) .* (-1.35703145e-07 .* T.^4 + 1.29111169e-04 .* T.^3 - 4.73210818e-02 .* T.^2 + 8.23639228e+00 .* T - 2.15281402e+02) .+
        (T .>= 300.00) .* (1.14074710e-10 .* T.^4 - 1.97122089e-07 .* T.^3 + 5.53525209e-05 .* T.^2 + 1.33834821e-01 .* T + 3.42764033e+02)
    )
end

function Copper_k(T, B, RRR)
    # from Jeroen v.N.
    # in units of W/(Km)
    a = 2.3797
    b = -0.4918
    c = -0.98615
    d = 0.13942
    e = 0.30475
    f = -0.019713
    g = -0.046897
    h = 0.0011969
    i = 0.0029988
    j = (@. (a + c * sqrt(T) + e * T + g * T^1.5 + i * T^2) / (1 + b * sqrt(T) + d * T + f * T^1.5 + h * T^2))
    k = (@. Copper_ρ(T, 0, RRR) / Copper_ρ(T, B, RRR))
        
    return @. 10^j * k
end

function Copper_ρ(T, B, RRR)
    # from Jeroen v.N.
    # ohm*m
    return @. ((1.467 ./ RRR + 1 ./ (2.32547e9 ./ T .^ 5 + 9.57137e5 ./ T .^ 3 + 1.62735e2 ./ T)) .* 1e-8 + B .* (0.37 + 0.0005 .* RRR) .* 1e-10)
end

# Hasteloy Material properties
function Hast_C(T)
    # from Bergen
    # in J/m^3-1K^-1
    
    a = @. ifelse(T < 35, -8.659744297212208e-7, 2.150134339872814e-10)
    b = @. ifelse(T < 35, 7.209308577841222e-5, -1.458398200879766e-7)
    c = @. ifelse(T < 35, -0.001346146407585, 4.376675190884510e-5)
    d = @. ifelse(T < 35, 0.0179929464336959, -0.01504610350977)
    e = @. ifelse(T < 35, 0.084834978274794, 4.491152093383801)
    f = @. ifelse(T < 35, 0.642508075273653, -110.7126478142060)
    
    return @. 8.89e3 * (a * T^5 + b * T^4 + c * T^3 + d * T^2 + e * T + f)
end

function Hast_k(T)
    # from Bergen
    # in W/m^-1K^-1   
    
    a = @. ifelse(T < 55, -1.976150097497916e-8, 3.544514666284345e-11)
    b = @. ifelse(T < 55, 1.692069531515161e-6, -3.085569423939295e-8)
    c = @. ifelse(T < 55, 3.428435748706821e-5, 9.862456406179634e-6)
    d = @. ifelse(T < 55, -0.008289421572708, -0.001448613557764)
    e = @. ifelse(T < 55, 0.399017559651135, 0.120321620954427)
    f = @. ifelse(T < 55, -0.827062339102484, 3.632700117022619)
    
    return @. a * T^5 + b * T^4 + c * T^3 + d * T^2 + e * T + f
end

function Hast_ρ(T)
    p1 = 0.0002667
    p2 = 1.234
    p3 = 0.0001469
    p4 = 1.229
    
    return @. ifelse(T <= 12, (p1 * T + p2) * 1e-6, (p3 * T + p4) * 1e-6)
end

## Silver Material Properties
function Silver_C(T)
    # B.J. McBride, S. Gordon and M.A. Reno, Thermodynamic Data for 
    # Fifty Reference Elements, NASA Technical Paper 3287 (1993) available 
    # online at https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20010021116.pdf 
    # and D.R. Smith, F.R. Fickett, Journal of Research National Institute 
    # Standards Technology, v100, p119 (1995)
    # Tmp near 962C (1763F), unique ID = 881|9|0|-1
    # data is in units of J/(kg-K)
    # T must be in degrees Kelvin for these equations
    
    result = similar(T, Float64)
    
    mask1 = (1.00 .<= T) .& (T .< 12.30)
    result[mask1] .= -3.51234917e-06 .* T[mask1].^5 .+ 1.66168020e-04 .* T[mask1].^4 .- 4.95765926e-04 .* T[mask1].^3 .+ 1.00959334e-02 .* T[mask1].^2 .- 1.39436904e-02 .* T[mask1] .+ 1.22422710e-02
    
    mask2 = (12.30 .<= T) .& (T .< 75.00)
    result[mask2] .= -8.43901503e-08 .* T[mask2].^5 .+ 3.38022515e-05 .* T[mask2].^4 .- 4.85354460e-03 .* T[mask2].^3 .+ 2.95644369e-01 .* T[mask2].^2 .- 4.66981261e+00 .* T[mask2] .+ 2.42484682e+01
    
    mask3 = (75.00 .<= T) .& (T .< 300.00)
    result[mask3] .= 2.41143501e-10 .* T[mask3].^5 .- 3.10537526e-07 .* T[mask3].^4 .+ 1.57045352e-04 .* T[mask3].^3 .- 3.96147773e-02 .* T[mask3].^2 .+ 5.17726486e+00 .* T[mask3] .- 6.38588439e+01
    
    mask4 = (300.00 .<= T) .& (T .<= 1235.00)
    result[mask4] .= -1.76849772e-08 .* T[mask4].^3 .+ 5.00714342e-05 .* T[mask4].^2 .+ 1.70570216e-02 .* T[mask4] .+ 2.25706549e+02
    
    result[.!(mask1 .| mask2 .| mask3 .| mask4)] .= 1.000000e+100
    
    return result .* 10.49e3
end

function Silver_k(T)
    # data is in units of W/(m-K)
    # taken from MPDB
    tc = @. ifelse((T >= 0.00) & (T < 12.00),
                    -5.30182898e-02 * T^6 + 1.83852143e+00 * T^5 - 2.01076316e+01 * T^4 + 5.51775003e+01 * T^3 - 7.96943362e+01 * T^2 + 3.98415514e+03 * T,
                    ifelse((T >= 12.00) & (T < 35.00),
                            -1.93171328e-04 * T^6 + 3.06222302e-02 * T^5 - 1.94183418e+00 * T^4 + 6.16763958e+01 * T^3 - 9.76813559e+02 * T^2 + 5.78952713e+03 * T + 1.16851851e+04,
                            ifelse((T >= 35.00) & (T < 100.00),
                                    1.62424109e-08 * T^6 - 8.86687517e-06 * T^5 + 2.03489646e-03 * T^4 - 2.52172418e-01 * T^3 + 1.78523078e+01 * T^2 - 6.87054377e+02 * T + 1.17412623e+04,
                                    ifelse((T >= 100.00) & (T < 173.00),
                                            -3.82464365e-05 * T^3 + 2.03715000e-02 * T^2 - 3.63423730e+00 * T + 6.47744954e+02,
                                            ifelse((T >= 173.00) & (T < 300.00),
                                                    -2.32630945e-08 * T^4 + 2.22554027e-05 * T^3 - 7.74648962e-03 * T^2 + 1.14373432e+00 * T + 3.70275575e+02,
                                                    ifelse((T >= 300.00) & (T <= 1235.00),
                                                            3.63030887e-17 * T^6 - 1.06327174e-14 * T^5 - 3.37202514e-10 * T^4 + 7.57664035e-07 * T^3 - 7.07574291e-04 * T^2 + 2.36992453e-01 * T + 4.03534347e+02,
                                                            1.000000e+100))))))
    return tc
end

function Silver_ρ(T)
    ## electrical resistivity in ohm*m
    # taken from MPDB
    rho = @. ifelse((T >= 1.00) & (T < 15.80),
                    6.14418300e-15 * T^3 - 6.69009400e-14 * T^2 + 2.25956700e-13 * T + 9.82204800e-12,
                    ifelse((T >= 15.80) & (T < 27.50),
                            2.02666700e-14 * T^3 - 6.16000000e-13 * T^2 + 7.47333300e-12 * T - 2.33000000e-11,
                            ifelse((T >= 27.50) & (T < 60.00),
                                    -1.20000000e-14 * T^3 + 2.21095200e-12 * T^2 - 7.58642900e-11 * T + 8.01547600e-10,
                                    ifelse((T >= 60.00) & (T < 200.00),
                                            6.84118400e-17 * T^3 - 4.44702800e-14 * T^2 + 6.97450800e-11 * T - 2.42874100e-09,
                                            ifelse((T >= 200.00) & (T <= 1235.00),
                                                    8.26904500e-18 * T^3 - 3.07705900e-15 * T^2 + 6.07474200e-11 * T - 1.81275200e-09,
                                                    1.000000e+100)))))
    return rho
end

## Tape Material Properties
function Tape_C(T)
    # returns Heat capacitace per meter [J/mK]
    return (HastTA.*Hast_C(T) + CopperA.*Copper_C(T)+ SilverA.*Silver_C(T)  +RebcoA.*Rebco_C(T))
end

function Tape_k(T)
    # returns conductance meter, [Wm/k]
    return (HastTA.*Hast_k(T) + CopperA.*Copper_k(T,0,100) + SilverA.*Silver_k(T)  +RebcoA.*Rebco_k(T))
end

function Tape_ρ(T)
    # returns resistance per meter tape [Ω/m]
    return inv.(HastTA./Hast_ρ(T) + CopperA./Copper_ρ(T,0,100) + SilverA./Silver_ρ(T))
end



# T = collect(1:1:300)

# ### Plot for the thermal conductance 
# traceTape   = scatter(x = T, y = Tape_k(T), name = "k_Tape")
# traceCopper = scatter(x = T, y = Copper_k(T,0,100).*CopperA , name = "k_Copper") 
# traceSilver = scatter(x = T, y = Silver_k(T).*SilverA , name = "k_Silver") 
# traceHast   = scatter(x = T, y = Hast_k(T).*HastTA , name = "k_Hasteloy") 
# traceRebco  = scatter(x = T, y = Rebco_k(T).*RebcoA, name = "k_ReBCO") 

# layout = Layout(
#     title = "Thermal conductance of various layers of the tape",
#     yaxis_title = "Thermal conductance [Wm/K]",
#     xaxis_title = "Temperature [K]",
#     yaxis_type = "log",
#     xaxis_type = "log",
#     font = attr(size = 25),
#     legend = attr(x = 0.7, y = 0.2)
# )

# plot([traceTape,traceCopper,traceHast,traceSilver,traceRebco],layout)


# ### Plot for the heat capacitace
# traceTape   = scatter(x = T, y = Tape_C(T), name = "C_Tape")
# traceCopper = scatter(x = T, y = Copper_C(T).*CopperA , name = "C_Copper") 
# traceSilver = scatter(x = T, y = Silver_C(T).*SilverA , name = "C_Silver") 
# traceHast   = scatter(x = T, y = Hast_C(T).*HastTA , name = "C_Hasteloy") 
# traceRebco  = scatter(x = T, y = Rebco_C(T).*RebcoA, name = "C_ReBCO") 

# layout = Layout(
#     title = "Heat capacitace of various layers of the tape",
#     yaxis_title = "Heat Capacity [J/mK]",
#     xaxis_title = "Temperature [K]",
#     yaxis_type = "log",
#     font = attr(size = 25),
#     legend = attr(x = 0.7, y = 0.2)
# )

# plot([traceTape,traceCopper,traceHast,traceSilver,traceRebco],layout)


### Plot for the resistivity
# traceTape = scatter(x = T, y = Tape_ρ(T), name = "ρ_Tape")
# traceCopper = scatter(x = T, y = Copper_ρ(T,0,100)./CopperA , name = "ρ_Copper") 
# traceSilver = scatter(x = T, y = Silver_ρ(T)./SilverA , name = "ρ_Silver") 
# traceHast   = scatter(x = T, y = Hast_ρ(T)./HastTA , name = "ρ_Hasteloy") 
# traceRebco  = 0

# layout = Layout(
#     title = "Resistance of various layers of the tape",
#     yaxis_title = "Resistivity [Ω/m]",
#     xaxis_title = "Temperature [K]",
#     yaxis_type = "log",
#     font = attr(size = 25),
#     legend = attr(x = 0, y = 0.9)
# )

# plot([traceTape,traceCopper,traceHast,traceSilver],layout)

