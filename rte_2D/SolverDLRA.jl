__precompile__

using ProgressMeter
using LinearAlgebra
using LegendrePolynomials
using QuadGK
using SparseArrays
using SphericalHarmonicExpansions,SphericalHarmonics,TypedPolynomials,GSL
using MultivariatePolynomials
using Einsum

include("PNSystem.jl")

struct SolverDLRA
    # spatial grid of cell interfaces
    x::Array{Float64};
    y::Array{Float64};

    # Solver settings
    settings::Settings;

    # preallocate memory for performance
    outRhs::Array{Float64,3};
    
    # squared L2 norms of Legendre coeffs
    gamma::Array{Float64,1};
    # Roe matrix
    AbsAx::Array{Float64,2};
    AbsAz::Array{Float64,2};
    # normalized Legendre Polynomials
    P::Array{Float64,2};
    # quadrature points
    mu::Array{Float64,1};
    w::Array{Float64,1};

    # functionalities of the PN system
    pn::PNSystem;

    L1x::SparseMatrixCSC{Float64, Int64};
    L1y::SparseMatrixCSC{Float64, Int64};
    L2x::SparseMatrixCSC{Float64, Int64};
    L2y::SparseMatrixCSC{Float64, Int64};

    # constructor
    function SolverDLRA(settings)
        x = settings.x;
        y = settings.y;

        # setup flux matrix
        gamma = zeros(settings.nPN+1);
        for i = 1:settings.nPN+1
            n = i-1;
            gamma[i] = 2/(2*n+1);
        end
        A = zeros(settings.nPN,settings.nPN);
            # setup flux matrix (alternative analytic computation)
        for i = 1:(settings.nPN-1)
            n = i-1;
            A[i,i+1] = (n+1)/(2*n+1)*sqrt(gamma[i+1])/sqrt(gamma[i]);
        end

        for i = 2:settings.nPN
            n = i-1;
            A[i,i-1] = n/(2*n+1)*sqrt(gamma[i-1])/sqrt(gamma[i]);
        end

        # construct PN system matrices
        pn = PNSystem(settings)
        SetupSystemMatrices(pn);

        outRhs = zeros(settings.NCellsX,settings.NCellsY,pn.nTotalEntries);

        # setup Roe matrix
        S = eigvals(pn.Ax)
        V = eigvecs(pn.Ax)
        AbsAx = V*abs.(diagm(S))*inv(V)

        S = eigvals(pn.Az)
        V = eigvecs(pn.Az)
        AbsAz = V*abs.(diagm(S))*inv(V)

        # compute normalized Legendre Polynomials
        Nq=200;
        (mu,w) = gauss(Nq);
        P=zeros(Nq,settings.nPN);
        for k=1:Nq
            PCurrent = collectPl(mu[k],lmax=settings.nPN-1);
            for i = 1:settings.nPN
                P[k,i] = PCurrent[i-1]/sqrt(gamma[i]);
            end
        end

        # setupt stencil matrix
        nx = settings.NCellsX;
        ny = settings.NCellsY;
        N = pn.nTotalEntries;
        L1x = spzeros(nx*ny,nx*ny);
        L1y = spzeros(nx*ny,nx*ny);
        L2x = spzeros(nx*ny,nx*ny);
        L2y = spzeros(nx*ny,nx*ny);

        # setup index arrays and values for allocation of stencil matrices
        II = zeros(3*(nx-2)*(ny-2)); J = zeros(3*(nx-2)*(ny-2)); vals = zeros(3*(nx-2)*(ny-2));
        counter = -2;

        for i = 2:nx-1
            for j = 2:ny-1
                counter = counter + 3;
                # x part
                index = vectorIndex(nx,i,j);
                indexPlus = vectorIndex(nx,i+1,j);
                indexMinus = vectorIndex(nx,i-1,j);

                II[counter+1] = index;
                J[counter+1] = index;
                vals[counter+1] = 2.0/2/settings.dx; 
                if i > 1
                    II[counter] = index;
                    J[counter] = indexMinus;
                    vals[counter] = -1/2/settings.dx;
                end
                if i < nx
                    II[counter+2] = index;
                    J[counter+2] = indexPlus;
                    vals[counter+2] = -1/2/settings.dx; 
                end
            end
        end
        L1x = sparse(II,J,vals,nx*ny,nx*ny);

        II .= zeros(3*(nx-2)*(ny-2)); J .= zeros(3*(nx-2)*(ny-2)); vals .= zeros(3*(nx-2)*(ny-2));
        counter = -2;

        for i = 2:nx-1
            for j = 2:ny-1
                counter = counter + 3;
                # y part
                index = vectorIndex(nx,i,j);
                indexPlus = vectorIndex(nx,i,j+1);
                indexMinus = vectorIndex(nx,i,j-1);

                II[counter+1] = index;
                J[counter+1] = index;
                vals[counter+1] = 2.0/2/settings.dy; 

                if j > 1
                    II[counter] = index;
                    J[counter] = indexMinus;
                    vals[counter] = -1/2/settings.dy;
                end
                if j < ny
                    II[counter+2] = index;
                    J[counter+2] = indexPlus;
                    vals[counter+2] = -1/2/settings.dy; 
                end
            end
        end
        L1y = sparse(II,J,vals,nx*ny,nx*ny);

        II = zeros(2*(nx-2)*(ny-2)); J = zeros(2*(nx-2)*(ny-2)); vals = zeros(2*(nx-2)*(ny-2));
        counter = -1;

        for i = 2:nx-1
            for j = 2:ny-1
                counter = counter + 2;
                # x part
                index = vectorIndex(nx,i,j);
                indexPlus = vectorIndex(nx,i+1,j);
                indexMinus = vectorIndex(nx,i-1,j);

                if i > 1
                    II[counter] = index;
                    J[counter] = indexMinus;
                    vals[counter] = -1/2/settings.dx;
                end
                if i < nx
                    II[counter+1] = index;
                    J[counter+1] = indexPlus;
                    vals[counter+1] = 1/2/settings.dx;
                end
            end
        end
        L2x = sparse(II,J,vals,nx*ny,nx*ny);

        II .= zeros(2*(nx-2)*(ny-2)); J .= zeros(2*(nx-2)*(ny-2)); vals .= zeros(2*(nx-2)*(ny-2));
        counter = -1;

        for i = 2:nx-1
            for j = 2:ny-1
                counter = counter + 2;
                # y part
                index = vectorIndex(nx,i,j);
                indexPlus = vectorIndex(nx,i,j+1);
                indexMinus = vectorIndex(nx,i,j-1);

                if j > 1
                    II[counter] = index;
                    J[counter] = indexMinus;
                    vals[counter] = -1/2/settings.dy;
                end
                if j < ny
                    II[counter+1] = index;
                    J[counter+1] = indexPlus;
                    vals[counter+1] = 1/2/settings.dy;
                end
            end
        end
        L2y = sparse(II,J,vals,nx*ny,nx*ny);

        new(x,y,settings,outRhs,gamma,AbsAx,AbsAz,P,mu,w,pn,L1x,L1y,L2x,L2y);
    end
end

function SetupIC(obj::SolverDLRA)
    u = zeros(obj.settings.NCellsX,obj.settings.NCellsY,obj.pn.nTotalEntries);
    u[:,:,1] = IC(obj.settings,obj.settings.xMid,obj.settings.yMid)
    return u;
end

function SolveUnconventionalAdaptive(obj::SolverDLRA)
    # Get rank
    r=50;
    rMaxTotal = Int(floor(obj.settings.r/2));
    s = obj.settings;
    # Set up initial condition and store as matrix
    v = SetupIC(obj);
    nx = obj.settings.NCellsX;
    ny = obj.settings.NCellsY;
    N = obj.pn.nTotalEntries
    u = zeros(nx*ny,N);
    for k = 1:N
        u[:,k] = vec(v[:,:,k]);
    end

    nT = Int(ceil(s.tEnd/s.dt))
    dt = s.dt

    prog = Progress(nT,1)

    # Low-rank approx of init data:
    X,S,W = svd(u);
    
    # rank-r truncation:
    X = X[:,1:r];
    W = W[:,1:r];
    S = Diagonal(S);
    S = S[1:r, 1:r];
    K = zeros(size(X));

    WAxW = zeros(r,r)
    WAzW = zeros(r,r)
    WAbsAxW = zeros(r,r)
    WAbsAzW = zeros(r,r)
    WeW = zeros(r,r)

    XL2xX = zeros(r,r)
    XL2yX = zeros(r,r)
    XL1xX = zeros(r,r)
    XL1yX = zeros(r,r)

    MUp = zeros(2*r,r)
    NUp = zeros(2*r,r)

    XNew = zeros(nx*ny,r)
    STmp = zeros(r,r)

    e1 = sparse([1],[1],[1.0],N,N); 

    rankInTime = zeros(2,nT);
    NormInTime = zeros(2,nT);

    prog = Progress(nT,1)
    t = 0.0;

    for n=1:nT
        rankInTime[1,n] = t;
        rankInTime[2,n] = r;
        NormInTime[1,n] = t;
        NormInTime[2,n] = norm(S,2);

        ################## K-step ##################
        K = X*S;

        WAzW = W'*obj.pn.Az'*W
        WAbsAzW = W'*obj.AbsAz'*W
        WAbsAxW = W'*obj.AbsAx'*W
        WAxW = W'*obj.pn.Ax'*W
        WeW = W'*e1*W;

        K = K .- dt*(obj.L2x*K*WAxW .+ obj.L2y*K*WAzW .+ obj.L1x*K*WAbsAxW .+ obj.L1y*K*WAbsAzW .+ s.sigmaT*K .- s.sigmaS*K*WeW);

        K = [K X];
        XNew,STmp = qr!(K);
        XNew = Matrix(XNew)
        XNew = XNew[:,1:2*r];

        MUp = XNew' * X;
        ################## L-step ##################
        L = W*S';

        XL2xX = X'*obj.L2x*X
        XL2yX = X'*obj.L2y*X
        XL1xX = X'*obj.L1x*X
        XL1yX = X'*obj.L1y*X

        L .= L .- dt*(obj.pn.Ax*L*XL2xX' .+ obj.pn.Az*L*XL2yX' .+ obj.AbsAx*L*XL1xX' .+ obj.AbsAz*L*XL1yX' .+ s.sigmaT*L .- s.sigmaS*e1*L);
                
        L = [L W];
        WNew,STmp = qr(L);
        WNew = Matrix(WNew)
        WNew = WNew[:,1:2*r];

        NUp = WNew' * W;
        W = WNew;
        X = XNew;
        ################## S-step ##################
        S = MUp*S*(NUp')

        XL2xX = X'*obj.L2x*X
        XL2yX = X'*obj.L2y*X
        XL1xX = X'*obj.L1x*X
        XL1yX = X'*obj.L1y*X

        WAzW = W'*obj.pn.Az'*W
        WAbsAzW = W'*obj.AbsAz'*W
        WAbsAxW = W'*obj.AbsAx'*W
        WAxW = W'*obj.pn.Ax'*W
        WeW = W'*e1*W;

        S .= S .- dt.*(XL2xX*S*WAxW .+ XL2yX*S*WAzW .+ XL1xX*S*WAbsAxW .+ XL1yX*S*WAbsAzW .+ s.sigmaT*S .- s.sigmaS*S*WeW);
        
        ################## truncate ##################

        # Compute singular values of S1 and decide how to truncate:
        U,D,V = svd(S);
        U = Matrix(U); V = Matrix(V)
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
        
        t += dt;
        next!(prog) # update progress bar
    end

    # return end time and solution
    return 0.5*sqrt(obj.gamma[1])*X*S*W',rankInTime;

end

function SolveUnconventionalAdaptiveHeun(obj::SolverDLRA)
    # Get rank
    r=50;
    rMaxTotal = Int(floor(obj.settings.r/2));
    s = obj.settings;
    # Set up initial condition and store as matrix
    v = SetupIC(obj);
    nx = obj.settings.NCellsX;
    ny = obj.settings.NCellsY;
    N = obj.pn.nTotalEntries
    u = zeros(nx*ny,N);
    for k = 1:N
        u[:,k] = vec(v[:,:,k]);
    end

    nT = Int(ceil(s.tEnd/s.dt))
    dt = s.dt

    prog = Progress(nT,1)

    # Low-rank approx of init data:
    X,S,W = svd(u);
    u = 1; v = 1;
    
    # rank-r truncation:
    X = X[:,1:r];
    W = W[:,1:r];
    S = Diagonal(S);
    S = S[1:r, 1:r];
    K = zeros(size(X));

    WAxW = zeros(r,r)
    WAzW = zeros(r,r)
    WAbsAxW = zeros(r,r)
    WAbsAzW = zeros(r,r)
    WeW = zeros(r,r)

    XL2xX = zeros(r,r)
    XL2yX = zeros(r,r)
    XL1xX = zeros(r,r)
    XL1yX = zeros(r,r)

    MUp = zeros(2*r,r)
    NUp = zeros(2*r,r)

    XNew = zeros(nx*ny,r)
    STmp = zeros(r,r)

    e1 = sparse([1],[1],[1.0],N,N); 

    rankInTime = zeros(2,nT);
    NormInTime = zeros(2,nT);

    prog = Progress(nT,1)
    t = 0.0;

    for n=1:nT
        rankInTime[1,n] = t;
        rankInTime[2,n] = r;
        NormInTime[1,n] = t;
        NormInTime[2,n] = norm(S,2);

        ################## K-step ##################
        K = X*S;

        WAzW = W'*obj.pn.Az'*W
        WAbsAzW = W'*obj.AbsAz'*W
        WAbsAxW = W'*obj.AbsAx'*W
        WAxW = W'*obj.pn.Ax'*W
        WeW = W'*e1*W;

        K1 = K .- dt*(obj.L2x*K*WAxW .+ obj.L2y*K*WAzW .+ obj.L1x*K*WAbsAxW .+ obj.L1y*K*WAbsAzW .+ s.sigmaT*K .- s.sigmaS*K*WeW);
        K2 = K1 .- dt*(obj.L2x*K1*WAxW .+ obj.L2y*K1*WAzW .+ obj.L1x*K1*WAbsAxW .+ obj.L1y*K1*WAbsAzW .+ s.sigmaT*K1 .- s.sigmaS*K1*WeW);

        K = [0.5*(K+K2) X];
        XNew,STmp = qr!(K);
        XNew = Matrix(XNew)
        XNew = XNew[:,1:2*r];

        MUp = XNew' * X;
        ################## L-step ##################
        L = W*S';

        XL2xX = X'*obj.L2x*X
        XL2yX = X'*obj.L2y*X
        XL1xX = X'*obj.L1x*X
        XL1yX = X'*obj.L1y*X

        L1 = L .- dt*(obj.pn.Ax*L*XL2xX' .+ obj.pn.Az*L*XL2yX' .+ obj.AbsAx*L*XL1xX' .+ obj.AbsAz*L*XL1yX' .+ s.sigmaT*L .- s.sigmaS*e1*L);
        L2 = L1 .- dt*(obj.pn.Ax*L1*XL2xX' .+ obj.pn.Az*L1*XL2yX' .+ obj.AbsAx*L1*XL1xX' .+ obj.AbsAz*L1*XL1yX' .+ s.sigmaT*L1 .- s.sigmaS*e1*L1);
                
        L = [0.5*(L+L2) W];
        WNew,STmp = qr(L);
        WNew = Matrix(WNew)
        WNew = WNew[:,1:2*r];

        NUp = WNew' * W;
        W = WNew;
        X = XNew;
        ################## S-step ##################
        S = MUp*S*(NUp')

        XL2xX = X'*obj.L2x*X
        XL2yX = X'*obj.L2y*X
        XL1xX = X'*obj.L1x*X
        XL1yX = X'*obj.L1y*X

        WAzW = W'*obj.pn.Az'*W
        WAbsAzW = W'*obj.AbsAz'*W
        WAbsAxW = W'*obj.AbsAx'*W
        WAxW = W'*obj.pn.Ax'*W
        WeW = W'*e1*W;

        S1 = S .- dt.*(XL2xX*S*WAxW .+ XL2yX*S*WAzW .+ XL1xX*S*WAbsAxW .+ XL1yX*S*WAbsAzW .+ s.sigmaT*S .- s.sigmaS*S*WeW);
        S2 = S1 .- dt.*(XL2xX*S1*WAxW .+ XL2yX*S1*WAzW .+ XL1xX*S1*WAbsAxW .+ XL1yX*S1*WAbsAzW .+ s.sigmaT*S1 .- s.sigmaS*S1*WeW);
        S .= 0.5*(S .+ S2);
        
        ################## truncate ##################

        # Compute singular values of S1 and decide how to truncate:
        U,D,V = svd(S);
        U = Matrix(U); V = Matrix(V)
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
        
        t += dt;
        next!(prog) # update progress bar
    end

    # return end time and solution
    return 0.5*sqrt(obj.gamma[1])*X*S*W',rankInTime;

end

function vectorIndex(nx,i,j)
    return (i-1)*nx + j;
end

function Vec2Mat(nx,ny,v)
    m = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            m[i,j] = v[(i-1)*nx + j]
        end
    end
    return m;
end