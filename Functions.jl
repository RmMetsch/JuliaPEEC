# using NLsolve
# using Roots

function FindIn(I,T)

    function ImplicitEq(x,I,Ic)
        return I - x - 100e-6*(x/Ic)^30 
    end
    
    Ic = Rebco_Ic(T)
    In = zero(Ic)

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
            if abs(Ih-Il) < 10e-7
                break
            end
        end
        # return 
        In[i] = I[i] - (Ih+Il)/2
    end

    return In

end
