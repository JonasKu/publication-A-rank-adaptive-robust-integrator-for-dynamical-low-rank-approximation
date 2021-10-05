__precompile__
include("quadrature.jl")
include("Basis.jl")
include("TimeSolver.jl")

using ProgressMeter
using LinearAlgebra

struct DLRSolver
    # spatial grid of cell interfaces
    x;

    # quadrature
    q::Quadrature;

    # DLRSolver settings
    settings::Settings;

    # spatial basis functions
    basis::Basis;

    # low-rank solution matrices
    X::Array{Float64,2}
    W::Array{Float64,2}
    S::Array{Float64,2}

    # preallocated matrices for Rhs
    A::Array{Float64,3}
    B::Array{Float64,3}
    Y::Array{Float64,3}
    Y1::Array{Float64,2}
    ACons::Array{Float64,3}

    fluxS::Array{Float64,2}
    fluxL::Array{Float64,2}

    yL::Array{Float64,2}
    yS::Array{Float64,2}
    yK::Array{Float64,2}

    # Dirichlet BCs
    uL::Array{Float64,1}
    uR::Array{Float64,1}

    # time solver
    rkUpdate::TimeSolver

    # constructor
    function DLRSolver(settings)
        x = settings.x;
        r = settings.r;
        q = Quadrature(settings.Nq,"Gauss");
        basis = Basis(q,settings);

        # note that these are actually the hat variables
        X = zeros(settings.Nx,r)
        S = zeros(r,r)
        W = zeros(settings.N,r) 

        A = zeros(r,r,r);
        B = zeros(r,r,settings.Nq^2);
        Y = zeros(r,r,r);
        Y1 = zeros(r,r);

        fluxS = zeros(r,r);
        fluxL = zeros(r,settings.N^2);
        yL = zeros(r,settings.N^2);
        yS = zeros(r,r);
        yK = zeros(settings.Nx,r);

        uL = zeros(settings.Nq^2);
        uR = zeros(settings.Nq^2);

        ACons = zeros(settings.NCons,settings.N^2,settings.N^2);

        rkUpdate = TimeSolver(settings);

        new(x,q,settings,basis,X,S,W,A,B,Y,Y1,ACons,fluxS,fluxL,yL,yS,yK,uL,uR,rkUpdate);
    end
end


function RhsK(obj::DLRSolver,K::Array{Float64,2},W::Array{Float64,2},r=obj.settings.r)
    Nx = obj.settings.Nx;
    Nq = obj.settings.Nq;
    N = obj.settings.N;
    fXi = 0.25;
    dt = obj.settings.dt;
    flux = zeros(r);
    dx = obj.settings.dx;

    WQuad = EvalAtQuad(obj.basis,W)';

    # Compute A_{i,j,m} = E[W_i W_j W_m]
    WQuad = EvalAtQuad(obj.basis,W)';
    for i = 1:r
        for j = 1:r
            for m = 1:r
                obj.A[i,j,m] = Integral(obj.q,WQuad[i,:].*WQuad[j,:].*WQuad[m,:].*fXi);
            end
        end
    end


    yK = zeros(Nx,r);

    for j = 2:(Nx-1)
        for p = 1:r
            flux[p] = 0.0;
            for l = 1:r
                for m = 1:r
                    flux[p] += 1/4/dx * (K[j+1,l]*K[j+1,m]-K[j-1,l]*K[j-1,m])*obj.A[l,m,p];
                end
            end
        end
        for p = 1:r
            yK[j,p] = 1/2/dt .* (K[j+1,p]-2*K[j,p]+K[j-1,p]) - flux[p];
        end
    end
    return yK;
end

function RhsS(obj::DLRSolver,X::Array{Float64,2},W::Array{Float64,2},L::Array{Float64,2},r=obj.settings.r)
    Nx = obj.settings.Nx;
    Nq = obj.settings.Nq^2;
    N = obj.settings.N^2;
    fXi = 0.25;
    dt = obj.settings.dt;
    # Compute A_{i,j,m} = E[L_i L_j phi_m]
    LQuad = EvalAtQuad(obj.basis,L)';
    WQuad = EvalAtQuad(obj.basis,W)';
    for i = 1:r
        for j = 1:r
            for m = 1:r
                obj.A[i,j,m] = Integral(obj.q,LQuad[i,:].*LQuad[j,:].*WQuad[m,:].*fXi);
            end
        end
    end

    for p = 1:r
        for l = 1:r
            obj.Y1[p,l] = 0.0;
            for j = 2:(Nx-1)
                if obj.settings.stabilization != 1
                    if obj.settings.stabilization == 2
                        obj.Y1[p,l] += -1/2/dt * X[j,p]*(X[j+1,l]-2*X[j,l]+X[j-1,l]);# DLR first, discretize second; stable version for projector splitting
                    else
                        obj.Y1[p,l] += 1/2/dt * X[j,p]*(X[j+1,l]-2*X[j,l]+X[j-1,l]); # discretize first, DLR second
                    end
                end
            end
        end
    end

    for p = 1:r
        for m = 1:r
            for l = 1:r
                obj.Y[m,l,p] = 0;
                for j = 2:(Nx-1) 
                    obj.Y[m,l,p] += 1/4/obj.settings.dx * X[j,p]*(X[j+1,l]*X[j+1,m]-X[j-1,l]*X[j-1,m]);
                end
            end
        end
    end

    yS = zeros(r,r);

    for q = 1:r
        for l = 1:r
            obj.fluxS[q,l] = 0.0;
            for p = 1:r
                for m = 1:r
                    obj.fluxS[q,l] += obj.Y[m,p,q]*obj.A[m,p,l];
                end
            end
        end
        for l = 1:r
            obj.yS[q,l] = 0.0
            for m = 1:r
                for i = 1:N
                    yS[q,l] += obj.Y1[q,m]*L[i,m]*W[i,l];
                end
            end
            yS[q,l] -= obj.fluxS[q,l];
        end
    end
    return yS;

end

function RhsL(obj::DLRSolver,X::Array{Float64,2},L::Array{Float64,2},recompute::Bool=true,r=obj.settings.r)
    Nx = obj.settings.Nx;
    Nq = obj.settings.Nq^2;
    N = obj.settings.N^2;
    fXi = 0.25;
    # Compute A_{i,j,m} = E[L_i L_j phi_m]
    LQuad = EvalAtQuad(obj.basis,L)';
    for i = 1:r
        for j = 1:r
            for m = 1:N
                obj.B[i,j,m] = Integral(obj.q,LQuad[i,:].*LQuad[j,:].*obj.basis.PhiQuad[:,m].*fXi);
            end
        end
    end

    if recompute
        for p = 1:r
            for m = 1:r
                for l = 1:r
                    obj.Y[m,l,p] = 0;
                    for j = 2:(Nx-1) 
                        obj.Y[m,l,p] += 1/4/obj.settings.dx * X[j,p]*(X[j+1,l]*X[j+1,m]-X[j-1,l]*X[j-1,m]);
                    end
                end
            end
        end

        for p = 1:r
            for l = 1:r
                obj.Y1[p,l] = 0.0;
                for j = 2:(Nx-1)
                    if obj.settings.stabilization != 1
                        obj.Y1[p,l] += 1/2/obj.settings.dt * X[j,p]*(X[j+1,l]-2*X[j,l]+X[j-1,l]);
                    end
                end
            end
        end
    end
    
    for p = 1:r
        for i = 1:N
            obj.fluxL[p,i] = 0.0;
            for l = 1:r
                for m = 1:r
                    obj.fluxL[p,i] += obj.Y[m,l,p]*obj.B[m,l,i];
                end
            end
        end
    end

    yL = zeros(r,N);

    for p = 1:r
        for i = 1:N
            yL[p,i] = 0.0;
            for l = 1:r
                yL[p,i] += obj.Y1[p,l]*L[i,l];
            end
            yL[p,i] -= obj.fluxL[p,i];
        end
    end

    return yL;
end

function SetupIC(obj::DLRSolver)
    u = zeros(obj.settings.N*obj.settings.N,obj.settings.Nx);
    uVals = zeros(obj.settings.Nq^2)
    for j = 1:obj.settings.Nx
        for q = 1:obj.settings.Nq
            for k = 1:obj.settings.Nq
                uVals[(k-1)*obj.settings.Nq+q] = obj.settings.IC(obj.x[j],obj.q.xi[k],obj.q.xi[q])[1];
            end
        end
        u[:,j] = ComputeMoments(obj.basis,uVals*0.25);
    end
    return u;
end

function SolveForward(obj::DLRSolver)
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.x[2]-obj.x[1];
    tEnd = obj.settings.tEnd;
    Nx = obj.settings.Nx;
    r = obj.settings.r; # DLR rank
    Nq = obj.settings.Nq^2;
    N = obj.settings.N^2; # here, N is the number of quadrature points

    fXi = 0.25;

    # Set up initial condition
    u = SetupIC(obj);

    # Low-rank approx of init data:
    X,S,W = svd(u'); 
    
    # rank-r truncation:
    X = X[:,1:r]; 
    W = W[:,1:r];
    S = Array(Diagonal(S));
    S = S[1:r, 1:r]; 

    K = zeros(Nx,r);
    K1 = zeros(Nx,r);
    KNew = X*S;
    L = zeros(N,r);
    L1 = zeros(N,r);
    S1 = zeros(r,r);

    PhiQuad = obj.basis.PhiQuad;
    
    Nt = Integer(round(tEnd/dt));
    
    # compute Dirichlet values if they are independent of time
    for k = 1:obj.settings.Nq
        for q = 1:obj.settings.Nq
            obj.uL[(k-1)*obj.settings.Nq+q] = obj.settings.IC(obj.x[1],obj.q.xi[k],obj.q.xi[q])[1];
            obj.uR[(k-1)*obj.settings.Nq+q] = obj.settings.IC(obj.x[end],obj.q.xi[k],obj.q.xi[q])[1];
        end
    end

    prog = Progress(Nt,1)
    
    # time loop
    #@gif 
    for n = 1:Nt

        ################## K-step ##################
        K .= X*S;

        # impose BCs
        WQuad = EvalAtQuad(obj.basis,W)';
        for i = 1:r
            K[1,i] = Integral(obj.q,obj.uL.*WQuad[i,:].*fXi);
            K[end,i] = Integral(obj.q,obj.uR.*WQuad[i,:].*fXi);
        end
        
        if obj.settings.rkType == "Heun"
            K1 .= K .+ dt*RhsK(obj,K,W);
            K1 .= K1 .+ dt*RhsK(obj,K1,W);
            K .= 0.5.*(K.+K1);
        elseif obj.settings.rkType == "Euler"
            K .= K .+ dt*RhsK(obj,K,W);
        elseif obj.settings.rkType == "SSP"
            K .= UpdateK(obj.rkUpdate,obj,K,W);
        end

        XNew,STmp = qr(K);
        XNew = Matrix(XNew)
        XNew = XNew[:,1:r];

        MUp = XNew' * X;

        ################## L-step ##################
        L .= W*S';

        if obj.settings.rkType == "Heun"
            L1 .= L .+ dt*RhsL(obj,X,L)';
            L1 .= L1 .+ dt*RhsL(obj,X,L1)';
            L .= 0.5.*(L.+L1);
        elseif obj.settings.rkType == "Euler"
            L .= L .+ dt*RhsL(obj,X,L)';
        elseif obj.settings.rkType == "SSP"
            L .= UpdateL(obj.rkUpdate,obj,X,L);
        end
                
        WNew,STmp = qr(L);
        WNew = Matrix(WNew)
        WNew = WNew[:,1:r];

        NUp = WNew' * W;

        W .= WNew;
        X .= XNew;

        ################## S-step ##################
        S .= MUp*S*(NUp')

        L = W*S';

        if obj.settings.rkType == "Heun"
            S1 .= S .+ dt.*RhsS(obj,X,W,L);
            L .= W*S';
            S1 .= S1 .+ dt.*RhsS(obj,X,W,L);
            S .= 0.5.*(S.+S1);
        elseif obj.settings.rkType == "Euler"
            S .= S .+ dt.*RhsS(obj,X,W,L);
        elseif obj.settings.rkType == "SSP"
            S .= UpdateS(obj.rkUpdate,obj,X,S,W,false);
        end
        
        next!(prog) # update progress bar

        t = t+dt;
    end

    # return end time and solution
    return t, X,S,W;
end

function SolveForwardAdaptive(obj::DLRSolver)
    truncateChristian = true;
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.x[2]-obj.x[1];
    tEnd = obj.settings.tEnd;
    Nx = obj.settings.Nx;
    r = obj.settings.r; # DLR rank
    Nq = obj.settings.Nq^2;
    N = obj.settings.N^2; # here, N is the number of quadrature points

    # start with max rank, all objects are preallocated with rank obj.settings.r
    rMaxTotal = Int(floor(obj.settings.r/2));
    r = rMaxTotal;#15;

    fXi = 0.25;

    # Set up initial condition
    u = SetupIC(obj);

    # Low-rank approx of init data:
    X,S,W = svd(u'); 
    
    # rank-r truncation:
    X = X[:,1:r]; 
    W = W[:,1:r];
    S = Array(Diagonal(S));
    S = S[1:r, 1:r]; 

    K = zeros(Nx,rMaxTotal);
    K1 = zeros(Nx,rMaxTotal);
    KNew = X*S;
    L = zeros(N,rMaxTotal);
    L1 = zeros(N,rMaxTotal);
    S1 = zeros(r,rMaxTotal);

    PhiQuad = obj.basis.PhiQuad;
    
    Nt = Integer(round(tEnd/dt));
    
    # compute Dirichlet values if they are independent of time
    for k = 1:obj.settings.Nq
        for q = 1:obj.settings.Nq
            obj.uL[(k-1)*obj.settings.Nq+q] = obj.settings.IC(obj.x[1],obj.q.xi[k],obj.q.xi[q])[1];
            obj.uR[(k-1)*obj.settings.Nq+q] = obj.settings.IC(obj.x[end],obj.q.xi[k],obj.q.xi[q])[1];
        end
    end

    prog = Progress(Nt,1)

    rankInTime = zeros(2,Nt);
    
    # time loop
    for n = 1:Nt
        rankInTime[1,n] = t;
        rankInTime[2,n] = r;

        ################## K-step ##################
        K = X*S;

        # impose BCs
        WQuad = EvalAtQuad(obj.basis,W)';
        for i = 1:r
            K[1,i] = Integral(obj.q,obj.uL.*WQuad[i,:].*fXi);
            K[end,i] = Integral(obj.q,obj.uR.*WQuad[i,:].*fXi);
        end
        
        K1 = K + dt*RhsK(obj,K,W,r);
        K1 = [K1 X];
        XNew,STmp = qr(K1);
        XNew = Matrix(XNew);

        XNew = XNew[:,1:2*r];

        MUp = XNew' * X;

        ################## L-step ##################
        L = W*S';

        L1 = L .+ dt*RhsL(obj,X,L,true,r)';
        L1 = [L1 W];
        WNew,STmp = qr(L1);
        WNew = Matrix(WNew);
       
        WNew = WNew[:,1:2*r];

        NUp = WNew' * W;

        ################## S-step ##################
        S = MUp*S*(NUp')

        L = WNew*S';

        S .= S .+ dt.*RhsS(obj,XNew,WNew,L,2*r);

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
        rmax = max(rmax,3);

        for l = 1:rmax
            S[l,l] = D[l];
        end

        # if 2*r was actually not enough move to highest possible rank
        if rmax == -1
            rmax = rMaxTotal;
        end

        # put SVD onto basis
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
    return t,rankInTime, X,S,W;
end

# SSP Update function for K-step
function UpdateK(obj::TimeSolver,solver::DLRSolver,K::Array{Float64,2},W::Array{Float64,2})
    NCells = solver.settings.Nx;
    obj.KRK[1,:,:] .= K;
    for s = 1:obj.rkStages
        obj.Krhs[s,:,:] .= RhsK(solver,obj.KRK[s,:,:],W);;
        obj.KRK[s+1,:,:] .= zeros(NCells,solver.settings.r);
        for j = 1:s
            obj.KRK[s+1,:,:] = obj.KRK[s+1,:,:]+obj.alpha[s,j].*obj.KRK[j,:,:]+obj.dt*obj.beta[s,j].*obj.Krhs[j,:,:];
        end
    end

    return obj.KRK[obj.rkStages+1,:,:];
end

# SSP Update function for S-step
function UpdateS(obj::TimeSolver,solver::DLRSolver,X::Array{Float64,2},S::Array{Float64,2},W::Array{Float64,2},backward::Bool=true)
    obj.SRK[1,:,:] .= S;
    for s = 1:obj.rkStages
        obj.Srhs[s,:,:] .= RhsS(solver,X,W,W*obj.SRK[s,:,:]');
        if backward # if projector splitting integrator is used, sign in rhsS must be changed
            obj.Srhs[s,:,:] .= -obj.Srhs[s,:,:];
        end
        obj.SRK[s+1,:,:] .= zeros(solver.settings.r,solver.settings.r);
        for j = 1:s
            obj.SRK[s+1,:,:] = obj.SRK[s+1,:,:]+obj.alpha[s,j].*obj.SRK[j,:,:]+obj.dt*obj.beta[s,j].*obj.Srhs[j,:,:];
        end
    end

    return obj.SRK[obj.rkStages+1,:,:];
end

# SSP Update function for L-step
function UpdateL(obj::TimeSolver,solver::DLRSolver,X::Array{Float64,2},L::Array{Float64,2},recompute::Bool=true)
    N = solver.settings.N^2;
    obj.LRK[1,:,:] .= L;
    for s = 1:obj.rkStages
        obj.Lrhs[s,:,:] .= RhsL(solver,X,obj.LRK[s,:,:],recompute)';
        obj.LRK[s+1,:,:] .= zeros(N,solver.settings.r);
        for j = 1:s
            obj.LRK[s+1,:,:] = obj.LRK[s+1,:,:]+obj.alpha[s,j].*obj.LRK[j,:,:]+obj.dt*obj.beta[s,j].*obj.Lrhs[j,:,:];
        end
    end

    return obj.LRK[obj.rkStages+1,:,:];
end