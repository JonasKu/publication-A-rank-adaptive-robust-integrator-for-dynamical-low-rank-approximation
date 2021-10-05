__precompile__
mutable struct Settings
    # grid settings
    # number spatial interfaces
    Nx::Int64;
    # number spatial cells
    NCells::Int64;
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
    # CFL number 
    cfl::Float64;
    
    # degree PN
    nPN::Int64;

    # solution values
    uL::Float64;
    uR::Float64;
    x0::Float64;
    x1::Float64;

    # spatial grid
    x
    xMid

    # problem definitions
    problem::String;

    # IC definitions
    ICType::String;
    BCType::String;

    # physical parameters
    sigmaT::Float64;
    sigmaS::Float64;

    # low rank parameters
    r::Int;
    useStabilizingTermsS::Bool;
    useStabilizingTermsL::Bool;
    stabilization::Int;

    epsAdapt::Float64;

    function Settings(Nx::Int=1002)
        # spatial grid setting
        #Nx = 100; # number of spatial grid points
        NCells = Nx - 1;
        a = -5#-1.5; # left boundary
        b = 5#1.5; # right boundary
        
        # time settings
        tEnd = 5;
        cfl = 1.0; # CFL condition

        # moment system parameters
        ICType = "LS"; # Possibilities are sin, shock, LS, ManufacturedSolution, shock

        # initial condition parameters
        uL = 12.0;
        uR = 1.0;
        x0 = 0.5;
        x1 = x0;
        
        # number PN moments
        nPN = 200; #100

        problem = "LineSource"; # Possibilities are ManufacturedSolution, ManufacturedSolutionLinear, ManufacturedSolutionSteady and LineSource and UQ

        x = collect(range(a,stop = b,length = NCells));
        dx = x[2]-x[1];
        x = [x[1]-dx;x]; # add ghost cells so that boundary cell centers lie on a and b
        x = x.+dx/2;
        println("NCells = ",NCells);
        println("size x = ",size(x));
        xMid = x[1:(end-1)].+0.5*dx
        # a lies on x(3), b lies on x(N-1). Cells 1,2 and N-1,N are ghost cells

        dt = cfl*dx;

        # stabilization method for PS: 0 - standard, 1 - Lax Wendroff in L and S, 2 - stabilization in all steps
        stabilization = 0;

        # physical parameters
        if problem =="LineSource"
            sigmaS = 1.0;#1.0;
            sigmaA = 0.0;
        else
            sigmaS = 0.4;
            sigmaA = 0.0;
        end
        sigmaT = sigmaA + sigmaS;

        BCType = "exact" # periodic,dirichlet,exact

        r = 60; #r/2 is the maximum rank;

        epsAdapt = 1e-1;

        # build class
        new(Nx,NCells,a,b,dx,tEnd,dt,cfl,nPN,uL,uR,x0,x1,x,xMid,problem,ICType,BCType,sigmaT,sigmaS,r,true,true,stabilization,epsAdapt);
    end

end

function IC(obj::Settings,x,xi=0.0)
    y = zeros(size(x));
    if obj.ICType == "shock"
        uL = 1;
        uR = 0.1;
        sigma = 0.0;
        x0 = -1.0+sigma*xi;
        x1 = 0.0+sigma*xi;
        for j = 1:length(y);
            if x[j] < x0
                y[j] = uR;
            elseif x[j] < x1
                y[j] = uL;
            else 
                y[j] = uR;
            end
        end
    elseif obj.ICType == "sin"
        for j = 1:length(y)
            y[j] = sin(2*pi*x[j]);
        end
    elseif obj.ICType == "LS"
        x0 = 0.0
        s1 = 0.03
        s2 = s1^2
        floor = 1e-4
        x0 = 0.0
        for j = 1:length(y);
            #println(1.0/(4.0*pi*s2) *exp(-((x[j]-x0)*(x[j]-x0))/4.0/s2))
            y[j] = max(floor,1.0/(sqrt(2*pi)*s1) *exp(-((x[j]-x0)*(x[j]-x0))/2.0/s2))
        end
    elseif obj.ICType == "ManufacturedSolution"
        for j = 1:length(y);
            y[j] = cos(x[j]);
        end
    end
    return y;
end

function IC0Deterministic(obj::Settings,x)
    y = zeros(size(x));
    for j = 1:length(y);
        if x[j] < obj.x0
            y[j] = obj.uL;
        else
            y[j] = obj.uR;
        end
    end
    return y;
end

function IC1Deterministic(obj::Settings,x)
    y = zeros(size(x));
    x0 = 0.0
    s2 = 0.03^2
    floor = 1e-4
    println("x = ",x)

    for j = 1:length(y);
        println(1.0/(4.0*pi*s2) *exp(-((x[j]-x0)*(x[j]-x0))/2))
        y[j] = max(floor,1.0/(sqrt(2*pi)*s2) *exp(-((x[j]-x0)*(x[j]-x0))/4.0/s2))
    end

    return y;
end

function IC3Deterministic(obj::Settings,x)
    y = zeros(size(x));
    for j = 1:length(y);
        if x[j] < obj.x0
            y[j] = obj.uL;
        elseif x[j] < obj.x1
            y[j] = 0.5*(obj.uR+obj.uL);            
        else
            y[j] = obj.uR;
        end
    end
    return y;
end

function IC1DeterministicExact(obj::Settings,t::Float64,x)
    y = zeros(size(x));

    if obj.problem == "Burgers"
        if t >= (obj.x1-obj.x0)/(obj.uL-obj.uR);
            tS = (obj.x1-obj.x0)/(obj.uL-obj.uR);
            x0BeforeShock = obj.x0 + tS*obj.uL;
            x1BeforeShock = obj.x1 + tS*obj.uR;
            x0 = x0BeforeShock + (t-tS)*(obj.uL+obj.uR)*0.5;
            x1 = x0 - 1.0;
        else
            x0 = obj.x0 + t*obj.uL;
            x1 = obj.x1 + t*obj.uR;
        end
        for j = 1:length(y);
            if x[j] < x0
                y[j] = obj.uL;
            elseif x[j] < x1
                y[j] = obj.uL + (obj.uR - obj.uL)*(x[j]-x0)/(x1-x0);
            else
                y[j] = obj.uR;
            end
        end
    elseif obj.problem == "Advection"
        if obj.ICType == "shock"
            x0 = obj.x0 + t*obj.advectionSpeed;
            x1 = obj.x1 + t*obj.advectionSpeed;
            for j = 1:length(y)
                if x[j] < x0
                    y[j] = obj.uL;
                elseif x[j] < x1
                    y[j] = obj.uL + (obj.uR - obj.uL)*(x[j]-x0)/(x1-x0);
                else
                    y[j] = obj.uR;
                end
            end
        elseif obj.ICType == "sin"
            for j = 1:length(y)
                y[j] = sin(2*pi*(x[j]-t*obj.advectionSpeed));
            end
        end
    end
    return y;
end

function ExactSolution(t::Float64,x::Float64,problem::String)
    if problem == "ManufacturedSolution"
        return exp(t)*cos(x)*2/sqrt(2);
    elseif problem == "ManufacturedSolutionLinear"
        return (t+1.0)*cos(x)*2/sqrt(2);
    elseif problem == "ManufacturedSolutionSteady"
        return cos(x)*2/sqrt(2);
    else
        return 0.0;
    end
end
