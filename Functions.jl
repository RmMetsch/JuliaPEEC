# using NLsolve
# using Roots
using BenchmarkTools
import .Rebco as Re

function FindIn(I,T)

    function ImplicitEq(x,I,Ic)
        return I - x - (x/Ic)^30 
    end
    
    Ic = Re.Ic(T)
    In = zero(Ic)
    println(Ic)

    for i in eachindex(In)

        Il = 0
        Ih = I[i]
        while true
            # bisection
            if ImplicitEq((Ih+Il)/2, I[i], Ic[i]) > 0
                Il = (Ih+Il)/2
            else
                Ih = (Ih+Il)/2
            end
            # break condition
            if abs(Ih-Il) < 10^(-7)
                break
            end
        end
        # return 
        In[i] = I[i] - (Ih+Il)/2
    end

    return In

end







 
