__precompile__
import FastTransforms
using FastGaussQuadrature

struct Quadrature
    Nq::Int;
    w::Array{Float64,1};
    xi::Array{Float64,1};

   function Quadrature(NqVal,quadType)
        if quadType == "Gauss"
            xi,w = gausslegendre(NqVal)
            #xi=(a*(1-xi)+b*(1+xi))/2.0;FastTransforms.clenshawcurtis
            #w = w*(2.0/(b-a));\
        elseif quadType == "ClenshawCurtis"
            xi,w = QuadNodesWeightsClCu(NqVal,-1.0,1.0)
        elseif quadType == "Grid"
            dxi = 2.0/(NqVal-1);
            w = dxi*ones(NqVal-1);
            xi = 1.0:-dxi:-1.0
        end
        xi = xi[end:-1:1];
        w = w[end:-1:1];

        new(NqVal,collect(w),collect(xi))
   end

end

function QuadNodesWeightsClCu(Nq::Int64,a::Float64,b::Float64) # transformation to a,b not tested!
    xiQuad,wQuad = FastTransforms.clenshawcurtis(Nq,0.,0.)
    xiQuad=(a*(1-xiQuad)+b*(1+xiQuad))/2.0;
    wQuad = wQuad*(2.0/(b-a));
    return xiQuad,wQuad;
end

function Integral(q::Quadrature, f)
    dot( q.w, f.(q.xi) );
end

function Integral(q::Quadrature, f, a::Float64, b::Float64)
    w = q.w*(b-a)*0.5;
    xi = (a*(1.0.-q.xi)+b*(1.0.+q.xi))./2.0;
    #return q.w'f.(xi);
    return sum(w'f.(xi))
end

function IntegralVec(q::Quadrature, fVec::Array{Float64,1}, a::Float64, b::Float64)
    w = q.w*(b-a)*0.5;
    return sum(w'fVec)
end

function IntegralMat(q::Quadrature, fVec::Array{Float64,2}, a::Float64, b::Float64)
    w = q.w*(b-a)*0.5;
    return fVec*w;
end

function IntegralMat(q::Quadrature, fVec::Array{Float64,2})
    return fVec*q.w;
end
