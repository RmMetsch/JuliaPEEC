module Matrices
using LinearAlgebra


# Heat transfer matrix
function BuildHeatTransMatrix(N::Int64,dx::Float64)

    # make base heattransfer matrix
    k_1 = ones(N-1)
    K = diagm(1 => k_1)
    K += K'
    K = K./dx^2

    for i in 1: 1 : length(K[:,1])
        K[i,i] = -sum(K[i,:])
    end
    
    K[1,:] .= 0
    K[end,:] .= 0
    
    return K
end

function BuildConductanceMatrix(N::Int64)
    g0 = N-1
    G = diagm(1=> -g0*ones(N-1))
    G += G'
    G += Diagonal([g0; 2*g0*ones(N-2);g0])
    return G
end

function BuildResistanceMatrix(N::Int64)
    R = Diagonal(ones(N-1)/(N-1))
    return R
end

function BuildThermalMassMatrix(N)
    return Diagonal([0; ones(N-2) ;0])
end


function BuildMkcl(N)
    Mkcl = zeros(N,N-1)
    # Mkcl[1:end,1:end] .= 1
    for i in 1:1:N-1
        Mkcl[i,i]   = -1
        Mkcl[i+1,i] = 1
    end

    return Mkcl
end


function BuildMkvl(N)
    return BuildMkcl(N)'
end


end