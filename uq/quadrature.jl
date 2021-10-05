__precompile__
import FastTransforms
using FastGaussQuadrature

struct Quadrature
    Nq::Int;
    w::Array{Float64,1};
    wTens::Array{Float64,1};
    xi::Array{Float64,1};

   function Quadrature(NqVal,quadType)
        if quadType == "Gauss"
            #xi,w = QuadNodesWeights(NqVal,-1.0,1.0)
            xi,w = gausslegendre(NqVal)
        elseif quadType == "ClenshawCurtis"
            xi,w = QuadNodesWeightsClCu(NqVal,-1.0,1.0)
        end
        xi = xi[end:-1:1];
        w = w[end:-1:1];

        wTens = zeros(NqVal^2);

        for k = 1:NqVal
            for q = 1:NqVal
                wTens[(q-1)*NqVal+k] = w[k]*w[q];
            end
        end


        new(NqVal,w,wTens,xi)
   end

end

function QuadNodesWeightsClCu(Nq::Int64,a::Float64,b::Float64) # transformation to a,b not tested!
    xiQuad,wQuad = FastTransforms.clenshawcurtis(Nq,0.,0.)
    xiQuad=(a*(1-xiQuad)+b*(1+xiQuad))/2.0;
    wQuad = wQuad*(2.0/(b-a));
    return xiQuad,wQuad;
end

function Integral(q::Quadrature, f::Function)
    dot( q.w, f.(q.xi) );
end

function Integral(q::Quadrature, f::Function, a::Float64, b::Float64)
    w = q.w*(b-a)*0.5;
    xi = (a*(ones(Float64,length(w)).-q.xi)+b*(ones(Float64,length(w))+q.xi))/2.0;
    #return q.w'f.(xi);
    return sum(w'f.(xi))
end

function IntegralVec(q::Quadrature, fVec::Array{Float64,1}, a::Float64, b::Float64)
    w = q.w*(b-a)*0.5;
    return sum(w'fVec)
end

function Integral(q::Quadrature, fVec::Array{Float64,1})
    return sum(q.wTens'fVec)
end

function IntegralMat(q::Quadrature, fVec::Array{Float64,2}, a::Float64, b::Float64)
    w = q.w*(b-a)*0.5;
    return fVec*w;
end

function IntegralMat(q::Quadrature, fVec::Array{Float64,2})
    return fVec*q.w;
end
