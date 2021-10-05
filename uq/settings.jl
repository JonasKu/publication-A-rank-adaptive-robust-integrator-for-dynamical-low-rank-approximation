__precompile__

mutable struct Settings
    # grid settings
    # number spatial interfaces
    Nx::Int64;
    # start and end point
    a::Float64;
    b::Float64;
    # grid cell width
    dx::Float64

    # time settings
    # end time
    tEnd::Float64;
    # time increment
    dt::Float64;
    # cfl number
    cfl::Float64;
    
    # number of quadrature points
    Nq::Int64;
    # definition of quadrature type
    quadratureType::String;

    # maximal polynomial degree
    N::Int64;

    # DLR rank
    r::Int64

    # solution values
    uL::Float64;
    uR::Float64;

    # initial shock position at xi = 0
    x0::Float64;

    # grid
    x

    sigma::Float64;

    # spatial limiter
    limiterType::String;
    # time update type
    rkType::String;
    # number of RK stages
    rkStages::Int;
    # filter
    filterType::String;
    lambda::Float64;
    filterOrder::Int64;

    useStabilizingTermsS::Bool;
    useStabilizingTermsL::Bool;
    stabilization::Int;

    # conservation settings
    NCons::Int; # number of conserved basis functions
    iCons; #::Array{Int64,2}; # indices of conserved basis funtions

    # epsilon for adaptive code
    epsAdapt::Float64;

    # initial condition
    IC::Function;
    solutionExact::Function;

    function Settings()
        # define spatial grid
        Nx = 600; #600
        a = 0.0;
        b = 1.0;
        x = range(a,stop = b,length = Nx)
        dx = (b-a)/(Nx-1.0);

        # define number of stochastic moments
        N = 20; #10

        # define DLR rank (if DLR is chosen)
        r = 80; #80

        # define time settings
        tEnd = 0.04;#0.03;
        cfl = 0.9;#0.9;

        # define test case settings
        x0 = 0.3; # 0.5
        x1 = 0.4;
        uL = 12.0;
        uM = 6.0;
        uR = 1.0;
        sigma0 = 0.2;#0.2;
        sigma1 = 5.0;#5.0;

        quadratureType = "Gauss"; # Possibilities are Gauss and ClenshawCurtis

        # compute time step size
        dt = cfl*dx/uL;

        # determine quadrature points
        Nq = ceil(1.5*N+1);

        # filter parameters
        filterType = "L2" # L2, EXP
        lambda = 0.0;#0.00001
        filterOrder = 1;

        # limiter type
        limiterType = "None";#"Minmod"

        # Runge-Kutta type Euler, Heun, SSP
        rkType = "Euler"
        rkStages = 1;

        useStabilizingTermsS = true;
        useStabilizingTermsL = true;

        # stabilization method for PS: 0 - standard, 1 - Lax Wendroff in L and S, 2 - stabilization in all steps
        stabilization = 0;

        # conservation settings
        NCons = 2;#3; # number of conserved basis functions
        iCons = [1 1; 1 2];#[1 1;1 2;1 3]; # indices of conserved basis funtions

        NCons = 0;#3; # number of conserved basis functions
        iCons = 0;

        epsAdapt = 1.0e-2;#1.5e-2;
        
        # build class 
        new(Nx,a,b,dx,tEnd,dt,cfl,Nq,quadratureType,N,r,uL,uR,x0,x,sigma0,limiterType,rkType,rkStages,filterType,lambda,filterOrder,useStabilizingTermsS,useStabilizingTermsL,stabilization,NCons,iCons,epsAdapt,
            #(xV,xi,eta)->IC3(xV,xi,eta,sigma0,sigma1,uL,uM,uR,x0,x1),
            #(t,xV,xi,eta)->IC3Exact(t,xV,xi,eta,sigma0,sigma1,uL,uM,uR,x0,x1))
            (xV,xi,eta)->IC1(xV,xi,eta,sigma0,sigma1,uL,uR,x0,x1),
            (t,xV,xi,eta)->IC1Exact(t,xV,xi,eta,sigma0,sigma1,uL,uR,x0,x1))
    end

end

function IC1(x,xi,eta,sigma0::Float64,sigma1::Float64,uL::Float64,uR::Float64,x0::Float64,x1::Float64)
    y = zeros(size(xi));
    
    for j = 1:length(y);
        uREta = uR+ sigma1*(eta[j]+1)*0.5;
        if x < x0+sigma0*xi[j]
            y[j] = uL;
        elseif x < x1+sigma0*xi[j]
            y[j] = uL + (uREta - uL)*(x-sigma0*xi[j]-x0)/(x1-x0);
        else
            y[j] = uREta;
        end
    end
    return y;
end

function IC1Exact(t::Float64,x,xi::Float64,eta::Float64,sigma0::Float64,sigma1::Float64,uL::Float64,uR::Float64,x0::Float64,x1::Float64)
    y = zeros(length(x));
    sigma = sigma0;

    for j = 1:length(y);
        uREta = uR+ sigma1*(eta[j]+1)*0.5;
        

        if t >= (x1-x0)/(uL-uREta);
            tS = (x1-x0)/(uL-uREta);
            x0BeforeShock = x0+sigma*xi + tS*uL;
            x1BeforeShock = x1+sigma*xi + tS*uREta;
            x0 = x0BeforeShock + (t-tS)*(uL+uREta)*0.5;
            x1 = x0 - 1.0;
        else
            x0 = x0+sigma*xi + t*uL;
            x1 = x1+sigma*xi + t*uREta;
        end
        
        if x[j] < x0
            y[j] = uL;
        elseif x[j] < x1
            y[j] = uL + (uREta - uL)*(x[j]-x0)/(x1-x0);
        else
            y[j] = uREta;
        end
    end

    return y;
end
