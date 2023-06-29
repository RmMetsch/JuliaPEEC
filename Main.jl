# import packages
using DifferentialEquations, Plots, LinearAlgebra, Symbolics,SparseArrays,BenchmarkTools
# import files
using .Matrices


# geometry
x_begin =  0
x_end   =  1
xrange  = x_end - x_begin
nodes = 11
dx = xrange/(nodes-1)
x = x_begin:dx:x_end


#  Matrices
M = Diagonal(ones(nodes))
K = buildHeatTransMat()

# time
tspan = (0.0,1)

# constants
p = 1

# inital conditon
T0 = ones(nodes)

# #Define the problem
function HeatEq(du, u, p, t)


    du[:] = K*u 
    du[1] = 0
    du[end] = 0


end


# # defining Jaconian
du0 = copy(T0)
u0  = copy(T0)


jac_sparsity = Symbolics.jacobian_sparsity((du,u) -> HeatEq(du, u, p, 0.0), du0, u0)

# problem with jacobian
f   = ODEFunction(HeatEq,
                  jac_prototype = float.(jac_sparsity),
                  mass_matrix = M) 


prob_jacobian  = ODEProblem(f,T0,tspan,p,)
sol = solve(prob_jacobian, TRBDF2())

plot(sol.u)