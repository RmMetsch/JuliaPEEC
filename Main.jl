# import packages
using DifferentialEquations, Plots, LinearAlgebra, Symbolics,SparseArrays,BenchmarkTools, Sundials, NLsolve
# import files 
include("Matrices.jl")
include("Functions.jl")
import .Rebco as Re
plotlyjs()

# geometry
x_begin =  0                # length in [m]
x_end   =  1
xrange  = x_end - x_begin
N  = 5                    
dx = xrange/(N-1)
x  = x_begin:dx:x_end
tapeWidth  = 0.01           # Take width in [m]

# constants
p = 1
r = 1*dx                    # factor is in [ohm/m], imes dx gives [ohm]

#  time independent linear Matrix
# submatrices
Mkvl =  BuildMkvl(N) 
Mkcl =  BuildMkcl(N)
G    =  zeros(N,N)
R    =  zeros(N-1,N-1)
# R    =  BuildResistanceMatrix(N)

# Thermal/electrical linear matrix
linThermMat = sparse( BuildHeatTransMatrix(N,dx))
LinElecMat  = sparse(hcat(vcat(G,Mkvl),vcat(Mkcl,R)))
# full Linear Matrix
linMatrix   = blockdiag(LinElecMat,linThermMat)
linMatrix   = linMatrix[2:end, 2:end]
# Mass Matrix
MassThermal    = sparse(BuildThermalMassMatrix(N))
MassElectrical = sparse(zeros(2N-1,2N-1)) 
MassMatrix = blockdiag(MassElectrical,MassThermal) 
MassMatrix = MassMatrix[2:end, 2:end]
# node to cell Matrix for Temperatures
Tn2c = Node2Cell(N)
Tc2n = Cell2Node(N)
Tc2n[[1,end]] .= 0


ext = [zeros(N-2); -1; zeros(2*N-1)]                        # external sources
tspan = (0, 1)                                             # Time

# #Define the problem
function HeatEq!(res,du, u, p, t)
    ### unpack the solution vector
    # V = u[1:N]
    T = u[2*N:end]
    I = u[N:2*N-2]
    ### add all contributions to the residual
    
    # linear contribution
    res[:]  = linMatrix*u .+ ext                            # Add external contribution    
    
    # Joule heating
    Tcell = Tn2c*T                                          # Convert to cells
    In    = FindIn(I,Tcell)                               # Determine In in cells
    
    P_nl  = Tc2n*(r.*In.^2)                                  # Convert back to nodes                                   
    res[:] += [zeros(N-1); zeros(N-1); P_nl]                # Add contribution to residual
    # time rate of change
    res[:] += -MassMatrix*du                                # Multiply by mass matrix
end

# setting intial condition
u0Elec  = vcat(collect(LinRange(0,-1,N)),ones(N-1))
u0Elec  = u0Elec[2:end]
u0Therm = ones(N)


u0 = vcat(u0Elec,u0Therm)
du0 = HeatEq!(zeros(length(u0)),zeros(length(u0)), u0, p,0)

# du0 = [zeros(length(u0Elec));ones(length(u0Therm))]
# solving the problem
difvars = Bool[zeros(N-1);zeros(N-1);[0;ones(N-2);0]]

prob = DAEProblem(HeatEq!, du0, u0, tspan, p, differential_vars = difvars) 
sol = solve(prob, IDA())

plot(sol.u)

