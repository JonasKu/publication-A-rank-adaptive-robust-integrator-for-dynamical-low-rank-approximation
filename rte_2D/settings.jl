__precompile__

mutable struct Settings
    # grid settings
    # number spatial interfaces
    Nx::Int64;
    Ny::Int64;
    # number spatial cells
    NCellsX::Int64;
    NCellsY::Int64;
    # start and end point
    a::Float64;
    b::Float64;
    c::Float64;
    d::Float64;
    # grid cell width
    dx::Float64
    dy::Float64

    # time settings
    # end time
    tEnd::Float64;
    # time increment
    dt::Float64;
    # CFL number 
    cfl::Float64;
    
    # degree PN
    nPN::Int64;

    # spatial grid
    x
    xMid
    y
    yMid

    # problem definitions
    problem::String;

    # physical parameters
    sigmaT::Float64;
    sigmaS::Float64;    

    # rank
    r::Int;
    epsAdapt::Float64;

    function Settings(Nx::Int=302,Ny::Int=302,r::Int=15,problem::String="LineSource")

        # spatial grid setting
        NCellsX = Nx - 1;
        NCellsY = Ny - 1;

        a = -1.5; # left boundary
        b = 1.5; # right boundary

        c = -1.5; # lower boundary
        d = 1.5; # upper boundary

        # physical parameters
        sigmaS = 1.0;
        sigmaA = 0.0;        
        sigmaT = sigmaA + sigmaS;

        # spatial grid
        x = collect(range(a,stop = b,length = NCellsX));
        dx = x[2]-x[1];
        x = [x[1]-dx;x]; # add ghost cells so that boundary cell centers lie on a and b
        x = x.+dx/2;
        xMid = x[1:(end-1)].+0.5*dx
        y = collect(range(c,stop = d,length = NCellsY));
        dy = y[2]-y[1];
        y = [y[1]-dy;y]; # add ghost cells so that boundary cell centers lie on a and b
        y = y.+dy/2;
        yMid = y[1:(end-1)].+0.5*dy

        # time settings
        tEnd = 1.0;
        cfl = 0.7#1.7 # CFL condition
        dt = cfl*dx;
        
        # number PN moments
        nPN = 39#39; # use odd number

        epsAdapt = 5e-2;

        # build class
        new(Nx,Ny,NCellsX,NCellsY,a,b,c,d,dx,dy,tEnd,dt,cfl,nPN,x,xMid,y,yMid,problem,sigmaT,sigmaS,r,epsAdapt);
    end
end

function IC(obj::Settings,x,y)
    x0 = 0.0;
    y0 = 0.0;
    out = zeros(length(x),length(y));
    s1 = 0.01
    s2 = 0.03^2
    floor = 1e-4
    for j = 1:length(x);
        for i = 1:length(y);
            out[j,i] = max(floor,1.0/(4.0*pi*s2) *exp(-((x[j]-x0)*(x[j]-x0)+(y[i]-y0)*(y[i]-y0))/4.0/s2))/4.0/pi;
        end
    end
    
    return out;
end