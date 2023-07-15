# import packages
using DifferentialEquations, PlotlyJS, LinearAlgebra, Symbolics,SparseArrays,BenchmarkTools, Sundials, NLsolve, LaTeXStrings
# import files 


N  =  101

include("Geometry.jl")
include("Matrices.jl")
include("Functions.jl")
include("Materials.jl")


# constants
p = 1
r = 1*dx                                                     # factor is in [ohm/m], imes dx gives [ohm]
Iapp = 4633.718019360891
tend = 1
tspan = (0, tend)                                            # Time
ext = [zeros(N-2); -1*Iapp; zeros(2*N-1)]                    # External sources


## event handling
# CenterTempIndex = round(Int,(2*N-1)+ (N-1)/2)
# condition(u,t,integrator) = t == 1
# affect!(integrator) = integrator.u[CenterTempIndex] += 10
# discreteCall = DiscreteCallback(condition,affect!)

# #Define the problem
function HeatEq!(res,du, u, p, t)
    ### unpack the solution vector
    # V = u[1:N-1]                                            # voltage vector
    u_I = u[N:2*N-2]                                          # current vector
    u_T = u[(2*N-1): end]                                     # temperature vector

    ### add all contributions to the residual
    # external contribution
    t <= tend ? res[:] = ext.*t : res[:] = ext   
    
    # linear contribution
    res[:] += linMatrix*u                                       # Add linear contribution    

    # Joule heating
    Tcell           = Tn2c*u_T                                            # Convert to cells
    In              = FindIn(u_I,Tcell)                                   # Find In in cells       
    P_nl            = Tc2n*(r.*In.^2)                                     # Find power generated inside cells due currentsharing
    P_nl[[1,end]]  .= 0                                         # Set power generated at boundaries to 0   
    res[2*N-1:end] += P_nl                                     # Add contribution to residual
    
    # Nonlinear contribution
    res[N:2*N-2] += r.*In                                       # Non linear cont of sc elements
    
    # time rate of change
    res[:] += -MassMatrix*du                                    # Multiply by mass matrix
end

# setting intial condition
# u0Elec  = vcat(collect(LinRange(0,-1,N)),Iapp*ones(N-1))
# u0Elec  = u0Elec[2:end]
u0Elec  = zeros(2*(N-1))
u0Therm = 20*ones(N)
u0      = vcat(u0Elec,u0Therm)
du0     = [zeros(N-1); Iapp*ones(N-1); zeros(N)]


# solving the problem
difvars = Bool[zeros(N-1);zeros(N-1); ones(N)]
prob    = DAEProblem(HeatEq!, du0, u0, tspan, p, differential_vars = difvars) 
sol = solve(prob, IDA(),tstops = tend)


## current ramping plot
Current = Vector{Float64}()
Voltage = Vector{Float64}()
tspan = 0:0.001:tend

for t in tspan
    push!(Current, sol(t)[N])
    push!(Voltage, abs(sol(t)[N-1]))
end

TraceVoltage = scatter(x = tspan, y= Voltage, mode = "lines", name = "Voltage", yaxis ="y2")
TraceCurrent = scatter(x = tspan, y= Current, mode = "lines", name = "Current")
layout  	 = Layout(title = "     Current and Voltage over time       Ic = 4633 A, Ec = 100 μ/m, R = 1",
                xaxis_title = "Time [sec]",
                yaxis_title = "Current [A]",
                yaxis2 = attr(title    = "Voltage Drop [V]",
                            overlaying = "y",
                            side       =  "right"),
                font = attr(size = 20,),
                legend = attr(x=0,y=1)) 

p = plot([TraceCurrent,TraceVoltage],layout) 


## temperature plot
# TraceArray = GenericTrace[]

# for t in 1:0.05:tend
#     Trace = scatter(x = x, y = sol(t)[2*N-1:end], mode = "lines", name = "t = $t")
#     push!(TraceArray, Trace)
# end
# layout = Layout(title = "Qench initiation",
#                 xaxis_title = "Distance [m]",
#                 yaxis_title = "Temperature [K]",
#                 font = attr(size = 15,),
#                 legend = attr(x=0.85,y=1))


# plot(TraceArray,layout)    

