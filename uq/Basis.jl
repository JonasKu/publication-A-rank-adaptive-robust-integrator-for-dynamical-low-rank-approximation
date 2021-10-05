__precompile__
import GSL

struct Basis
    # number of moments
    N::Int64;
    Nq::Int64;
    states::Int64;

    # precomputed Legendre Polynomials at quadrature points
    PhiQuad::Array{Float64,2};
    PhiQuadW::Array{Float64,2};

    function Basis(quadrature::Quadrature,settings::Settings)
        N = settings.N;
        Nq = quadrature.Nq;

        # precompute Legendre basis functions at quad points
        PhiQuad = zeros(quadrature.Nq^2,N^2);
        PhiQuadW = zeros(N^2,quadrature.Nq^2);
        for i = 1:N
            for j = 1:N
                for k = 1:Nq
                    for q = 1:Nq
                        PhiQuad[(q-1)*Nq+k,(j-1)*N+i] = Phi.(i-1,quadrature.xi[k])*Phi.(j-1,quadrature.xi[q]);
                        PhiQuadW[(j-1)*N+i,(q-1)*Nq+k] = Phi.(i-1,quadrature.xi[k])*Phi.(j-1,quadrature.xi[q])*quadrature.w[k]*quadrature.w[q];
                    end
                end
            end
        end

        new(N,Nq,1,PhiQuad,PhiQuadW);
    end

end

# Legendre Polynomials on [-1,1]
function Phi(n::Int64,xi::Float64) ## can I use this with input vectors xVal?
    return sqrt(2.0*n+1.0)*GSL.sf_legendre_Pl.(n,xi);
end

# evaluate polynomial with moments u at spatial position x
function Eval(obj::Basis,u::Array{Float64,1},xi)
    y=zeros(length(xi))
    for i = 1:obj.N
        y = y+u[i]*Phi.(i-1,xi);
    end
    return y;
end

# evalueates states at quadrature points. returns solution[quadpoints,states]
function Eval(obj::Basis,u::Array{Float64,2},xi)
    y=zeros(length(xi),obj.states)
    for s = 1:obj.states
        for i = 1:obj.N
            y[:,s] = y[:,s]+u[i,s]*Phi.(i-1,xi);
        end
    end
    return y;
end

function EvalAtQuad(obj::Basis,u::Array{Float64,1})
    return obj.PhiQuad*u;
end

# evalueates states at quadrature points. returns solution[quadpoints,states]
function EvalAtQuad(obj::Basis,u::Array{Float64,2})
    return obj.PhiQuad*u;
end

# returns N moments of a function evaluated at the quadrature points and stored in uQ
function ComputeMoments(obj::Basis,uQ::Array{Float64,1})
    return obj.PhiQuadW*uQ;
end

# returns N moments of all states evaluated at the quadrature points and stored in uQ
# uQ in [Nq,states], returns [N,states]
function ComputeMoments(obj::Basis,uQ::Array{Float64,2})
    return obj.PhiQuadW*uQ;
end

# evaluate polynomial with moments u at spatial position x
function Eval(obj::Basis,u::Array{Float64,1},xi,eta)
    y=zeros(length(xi))
    for i = 1:obj.N
        for j = 1:obj.N
            y = y.+u[(j-1)*obj.N+i]*Phi.(i-1,xi).*Phi.(j-1,eta);
        end
    end
    return y;
end