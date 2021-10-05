__precompile__
include("quadrature.jl")
include("basis.jl")

using ProgressMeter
using LinearAlgebra

struct DLRSolver
    # spatial grid of cell interfaces
    x::Array{Float64,1};

    # Solver settings
    settings::Settings;

    # preallocate memory for performance
    outRhs::Array{Float64,2};
    
    # squared L2 norms of Legendre coeffs
    gamma::Array{Float64,1};
    # flux matrix PN system
    A::Array{Float64,2};

    # physical parameters
    sigmaT::Float64;
    sigmaS::Float64;

    # constructor
    function DLRSolver(settings)
        x = settings.x;

        outRhs = zeros(settings.NCells,settings.nPN);

        # setup gamma vector
        gamma = zeros(settings.nPN);
        for i = 1:settings.nPN
            n = i-1;
            gamma[i] = 2/(2*n+1);
        end
        
        # setup flux matrix
        A = zeros(settings.nPN,settings.nPN)
        if settings.problem == "UQ" # UQ for advection equation
            N = settings.nPN;
            q = Quadrature(2*N,"Gauss");
            b = Basis(q,settings);
            aXi = (1.0 - settings.sigmaS) .+settings.sigmaS*q.xi;
            for i = 1:N
                for j = 1:N
                    A[i,j] = IntegralVec(q, aXi.*b.PhiQuad[:,i].*b.PhiQuad[:,j]*0.5,-1.0,1.0)
                end
            end
        else # radiative transfer
            for i = 1:(settings.nPN-1)
                n = i-1;
                A[i,i+1] = (n+1)/(2*n+1)*sqrt(gamma[i+1])/sqrt(gamma[i]);
            end

            for i = 2:settings.nPN
                n = i-1;
                A[i,i-1] = n/(2*n+1)*sqrt(gamma[i-1])/sqrt(gamma[i]);
            end
        end

        # update dt with correct maximal speed lmax
        lmax = maximum(abs.(eigvals(A)));
        settings.dt = settings.dx*settings.cfl/lmax;

        new(x,settings,outRhs,gamma,A,settings.sigmaT,settings.sigmaS);
    end
end

function SetupIC(obj::DLRSolver)
    u = zeros(obj.settings.NCells,obj.settings.nPN); # Nx interfaces, means we have Nx - 1 spatial cells
    u[:,1] = 2.0/sqrt(obj.gamma[1])*IC(obj.settings,obj.settings.xMid);
    return u;
end

function SolveUnconventionalEfficient(obj::DLRSolver)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;
    dx = obj.x[2]-obj.x[1];
    Nx = obj.settings.NCells;
    r = obj.settings.r; # DLR rank
    N = obj.settings.nPN; # here, N is the number of quadrature points

    rMaxTotal = Int(floor(obj.settings.r/2));
    r = 2;#15;


    # Set up initial condition
    u = SetupIC(obj);

    # Low-rank approx of init data:
    X,S,W = svd(u); 
    
    # rank-r truncation:
    X = X[:,1:r]; 
    W = W[:,1:r];
    S = Diagonal(S);
    S = S[1:r, 1:r]; 

    K = zeros(Nx,r);
    KNew = X*S;
    L = zeros(N,r);
    
    Nt = Integer(round(tEnd/dt));

    rankInTime = zeros(2,Nt);
    NormInTime = zeros(2,Nt);

    prog = Progress(Nt,1)
    for n = 1:Nt
        rankInTime[1,n] = t;
        rankInTime[2,n] = r;
        NormInTime[1,n] = t;
        NormInTime[2,n] = norm(S,2);

        ############## Streaming ##############

        ###### K-step ######
        K = X*S;

        y = zeros(Nx,r);
        dt = obj.settings.dt;
        dx = obj.settings.dx;
        WAW = W'*obj.A*W;
        for j = 2:(obj.settings.NCells-1) # leave out ghost cells
            y[j,:] = (K[j+1,:]-2*K[j,:]+K[j-1,:])/dt/2- (WAW*(K[j+1,:]-K[j-1,:])/dx/2);
            # add source term and scattering
            y[j,:] += - obj.sigmaT*K[j,:];
            y[j,:] += (obj.sigmaS*K[j,:]*W[1,:]')'*W[1,:]
        end

        K1 = K + dt*y;
        K1 = [K1 X];

        XNew,STmp = qr(K1);
        XNew = Matrix(XNew);
        XNew = XNew[:,1:2*r];

        MUp = XNew' * X;
        ###### L-step ######
        L = W*S';

        dt = obj.settings.dt;
        dx = obj.settings.dx;

        X1 = zeros(r,r);
        X2 = zeros(r,r);
        
        for k = 1:r
            for l = 1:r
                for j = 2:(obj.settings.NCells-1)
                    X1[k,l] = X1[k,l] + X[j,k].*(X[j+1,l]-2*X[j,l]+X[j-1,l])/dt/2;
                    X2[k,l] = X2[k,l] + X[j,k].*(X[j+1,l]-X[j-1,l])/dx/2;
                end
            end
        end

        y = L*X1' - obj.A*L*X2' - obj.sigmaT*L;
        y[1,:] += obj.sigmaS*L[1,:]

        L1 = L .+ dt*y;

        L1 = [L1 W];
        WNew,STmp = qr(L1);
        WNew = Matrix(WNew);
        WNew = WNew[:,1:2*r];

        NUp = WNew' * W;

        ################## S-step ##################
        # use Lemma 1 to note that we can use rank-r solution inside right-hand side to reduce computational costs

        X1 = zeros(2*r,r);
        X2 = zeros(2*r,r);
        
        for k = 1:2*r
            for l = 1:r
                for j = 2:(obj.settings.NCells-1)
                    X1[k,l] = X1[k,l] + XNew[j,k].*(X[j+1,l]-2*X[j,l]+X[j-1,l])/dt/2;
                    X2[k,l] = X2[k,l] + XNew[j,k].*(X[j+1,l]-X[j-1,l])/dx/2;
                end
            end
        end

        dt = obj.settings.dt;

        ySV = obj.sigmaS*MUp*S*(WNew[1,:]'.*W[1,:])

        STilde = MUp*S*(NUp');

        S = STilde .+ dt*X1*S*NUp' - dt*X2*L'*obj.A*WNew .- dt*obj.sigmaT*STilde .+ dt*ySV

        ################## truncate ##################

        # Compute singular values of S1 and decide how to truncate:
        U,D,V = svd(S);
        rmax = -1;
        S .= zeros(size(S));


        tmp = 0.0;
        tol = obj.settings.epsAdapt*norm(D);
        
        rmax = Int(floor(size(D,1)/2));
        
        for j=1:2*rmax
            tmp = sqrt(sum(D[j:2*rmax]).^2);
            if(tmp<tol)
                rmax = j;
                break;
            end
        end
        
        rmax = min(rmax,rMaxTotal);
        rmax = max(rmax,2);

        for l = 1:rmax
            S[l,l] = D[l];
        end

        # if 2*r was actually not enough move to highest possible rank
        if rmax == -1
            rmax = rMaxTotal;
        end

        # update solution with new rank
        XNew = XNew*U;
        WNew = WNew*V;

        # update solution with new rank
        S = S[1:rmax,1:rmax];
        X = XNew[:,1:rmax];
        W = WNew[:,1:rmax];

        # update rank
        r = rmax;
      
        next!(prog) # update progress bar

        t = t+dt;
    end

    # return end time and solution
    return t,rankInTime,NormInTime, 0.5*sqrt(obj.gamma[1])*X*S*W';

end
