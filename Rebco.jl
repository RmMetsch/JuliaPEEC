module Rebco

    function Ic(T)
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
        t = T ./ T_c0

        Bi_ab = Bi0_ab * ((1 .- t .^ n1) .^ n2 + a * (1 .- t .^ n_ab))
        Bi_c = Bi0_c * (1 .- t .^ n_c)

        b_ab = B ./ Bi_ab
        b_c = B ./ Bi_c

        Jc_ab = ifelse.((t .< 1) .& (b_ab .< 1), (alpha_ab ./ B) .* (b_ab .^ p_ab) .* (1 .- b_ab) .^ q_ab .* ((1 .- t .^ n1) .^ n2 + a * (1 .- t .^ n_ab)) .^ gamma_ab, 0)
        Jc_c = ifelse.((t .< 1) .& (b_c .< 1), (alpha_c ./ B) .* (b_c .^ p_c) .* (1 .- b_c) .^ q_c .* (1 .- t .^ n_c) .^ gamma_c, 0)

        g = g0 .+ g1 .* exp.(-g2 .* B .* exp.(g3 .* T))

        Jc = min.(Jc_ab, Jc_c) .+ max.(0, Jc_ab .- Jc_c) ./ (1 .+ ((angle .- pi / 2) ./ g) .^ nu)

        # Convert units
        # thickness of superconducting layer
        dsup = 2e-6  # from J. Fleiter report 
        return Jc .* w .* dsup   # crittical current [A]

    end

    function C(T)
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

    function k(T)
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

end