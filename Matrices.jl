# Heat transfer matrix
function BuildHeatTransMatrix(T::Vector{Float64},dx::Float64)
    # make base heattransfer matrix
    k_1 = ones(N-1)
    K = diagm(1 => k_1)
    K += K'
    K = K./dx^2

    # # make base heattransfer matrix
    # kOffDiagonal = Tape_k(Node2Cell(T))./dx^2  
    # K = diagm(1 => kOffDiagonal, -1=> kOffDiagonal)
    
    # kDiagonal = -(pushfirst!(kOffDiagonal,0) + popfirst!(push!(kOffDiagonal,0)))
    # K += diagm(0 => kDiagonal)

    for i in 1: 1 : length(K[:,1])
        K[i,i] = -sum(K[i,:])
    end

    K[1,:] .=  0
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

function BuildThermalMassMatrix(N,dx::Float64)
    return Diagonal(ones(N))
    # return Diagonal(Tape_C(T))
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

function Node2Cell(N)
    return spdiagm(N-1,N, 0=> 1/2*ones(N-1),1=>1/2*ones(N-1))
end

function Cell2Node(N)
    return Node2Cell(N)' 
end

# Time independent linear Matrix'
Mkvl =  BuildMkvl(N) 
Mkcl =  BuildMkcl(N)
G    =  zeros(N,N)
R    =  zeros(N-1,N-1)

# Thermal/electrical linear matrix
linThermMat = sparse(BuildHeatTransMatrix(N,dx))
LinElecMat  = sparse(hcat(vcat(G,Mkvl),vcat(Mkcl,R)))

# full Linear Matrix
linMatrix   = blockdiag(LinElecMat,linThermMat)
linMatrix   = linMatrix[2:end, 2:end]

# Mass Matrix
MassThermal    = 1/dx.* sparse(BuildThermalMassMatrix(N))
MassElectrical = sparse(zeros(2N-1,2N-1)) 
MassMatrix = blockdiag(MassElectrical,MassThermal) 
MassMatrix = MassMatrix[2:end, 2:end]

# node to cell Matrix for Temperatures
Tn2c = Node2Cell(N)
Tc2n = Cell2Node(N)


