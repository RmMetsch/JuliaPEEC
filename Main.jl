# import packages
using DifferentialEquations, Plots, LinearAlgebra, Symbolics,SparseArrays,BenchmarkTools, Sundials
# import files 
include("Matrices.jl")
plotlyjs()

# geometry
x_begin =  0
x_end   =  1
xrange  = x_end - x_begin
N  = 101
dx = xrange/(N-1)
x  = x_begin:dx:x_end

# constants
p = 1
r = 1*dx    #factor is in [ohm/m], imes dx gives [ohm]

#  time independent linear Matrix
# submatrices
Mkvl =  BuildMkvl(N) 
Mkcl =  BuildMkcl(N)
G    =  zeros(N,N)
R    =  BuildResistanceMatrix(N)
# Thermal/electrical linear matrix
linThermMat = sparse( BuildHeatTransMatrix(N,dx))
LinElecMat  = sparse(hcat(vcat(G,Mkvl),vcat(Mkcl,R)))
# full Linear Matrix
linMatrix   = blockdiag(LinElecMat,linThermMat)
linMatrix   = linMatrix[2:end, 2:end]
# Mass Matrix
MassThermal    = sparse( BuildThermalMassMatrix(N))
MassElectrical = sparse(zeros(2N-1,2N-1)) 
MassMatrix = blockdiag(MassElectrical,MassThermal) 
MassMatrix = MassMatrix[2:end, 2:end]
# # joule conduction matrix
a = [0;1/2*ones(N-2)]
JouleHeatFlowMatrix = spdiagm(N,N-1,0=>a,-1=>reverse(a) )
# external sources
ext = [zeros(N-2); -1; zeros(2*N-1)] 
# time
tspan = (0, 10)

# #Define the problem
function HeatEq!(res,du, u, p, t)
    ### unpack the solution vector
    # V = u[1:N]
    # T = u[2*N:end]
    I = u[N:2*N-2]
    ### add all contributions to the residual
    # #linear contribution
    res[:]  = linMatrix*u .+ ext
    # time rate of change
    res[:] += -MassMatrix*du
    # nonlinear contribution'
    res[:] += [zeros(N-1);zeros(N-1); 100*r.*JouleHeatFlowMatrix*(I.^2)]
    
end

# setting intial condition
u0Elec  = vcat(collect(LinRange(0,-1,N)),ones(N-1))
u0Elec  = u0Elec[2:end]
u0Therm = zeros(N)
# u0Therm = [0.25-(x-0.5)^2 for x in x] 

u0 = vcat(u0Elec,u0Therm)
du0 = HeatEq!(zeros(length(u0)),zeros(length(u0)), u0, p,0)

# du0 = [zeros(length(u0Elec));ones(length(u0Therm))]
# solving the problem
difvars = Bool[zeros(N-1);zeros(N-1);[0;ones(N-2);0]]

prob = DAEProblem(HeatEq!, du0, u0, tspan, p, differential_vars = difvars) 
sol = solve(prob, IDA())

plot(sol.u)


# jac_sparsity = Symbolics.jacobian_sparsity((du,u) -> HeatEq!(1, du, u, p, 0.0), du0, u0)
# # problem with jacobian