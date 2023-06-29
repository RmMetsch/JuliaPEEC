module Matrices
using LinearAlgebra



# Heat transfer matrix
function buildHeatTransMat(nodes,dx)

    # make base heattransfer matrix
    k_1 = ones(nodes-1)
    K = diagm(1 => k_1)
    K += K'
    K = K./dx^2

    for i in 1: 1 : length(K[:,1])
        K[i,i] = -sum(K[i,:])
    end

    return K
end


function buildMkcl(N)
    Mkcl = zeros(N,N-1)
    # Mkcl[1:end,1:end] .= 1
    for i in 1:1:N-1
        Mkcl[i,i]   = -1
        Mkcl[i+1,i] = 1
    end

    return Mkcl
end


function BuildMkvl(N)
    return buildMkcl(N)'
end

end