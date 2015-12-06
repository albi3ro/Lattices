module HAtom
    export Ψf
    
    type Orbital
        r::Array{Float64}
        θ::Array{Float64}
        ϕ::Array{Float64}
        n::Int64
        l::Int64
        m::Int64
        Ψ:Array{Float64}
    end
    
    using GSL

    function Yml(m,l,θ,ϕ)
        return (-1)^m*sf_legendre_Plm(l,m,cos(θ))*e^(im*m*ϕ)
    end
    
    function R(n,l,ρ)
        return sf_laguerre_n(n-l-1,2*l+1,ρ)*e^(-ρ/2)*ρ^l
    end

    function Ψf(r,θ,ϕ,n::Int64,l::Int64,m::Int64)
        if(l>(n-1) || l<0)
            error("angular momentum l not in bounds 0<l<n-1")
        end
        if(abs(m)>l)
            error("magnetic moment not in bounds abs(m)<l")
        end
        
        ρ=(2r)/(n);
        
        Norm=sqrt((2/n)^3 * factorial(n-l-1)/(2n*factorial(n+l)));
        
        Ψn=Norm*R(n,l,ρ)*Yml(m,l,θ,ϕ);
        
        return Ψn
    end

    function 
        r=0:.01:10;
        θ=0:.01:2π;
        ϕ=0:.01:π;
        Ψ=Array{Float64}(length(r),length(θ),length(ϕ))
        
        for i in 1:length(r)
            for j in 1:length(θ)
                for k in 1:length(ϕ)
                    Ψ
            
                end
            end
        end
        
    
    end

    function ΨSurf              
    
    
    end 

end
