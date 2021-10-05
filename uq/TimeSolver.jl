__precompile__
struct TimeSolver
    alpha::Array{Float64,2};
    beta::Array{Float64,2};
    gamma::Array{Float64,1};
    rkStages::Int64;
    Nx::Int64;
    dt::Float64;

    KRK::Array{Float64,3};
    SRK::Array{Float64,3};
    LRK::Array{Float64,3};
    Krhs::Array{Float64,3};
    Srhs::Array{Float64,3};
    Lrhs::Array{Float64,3};

    # constructor
    function TimeSolver(settings::Settings)
        rkStages = settings.rkStages;
        Nx = settings.Nx;
        dt = settings.dt;

        # set up alpha and beta values for SSP RK method
        alpha, beta = SetUpRKTable(rkStages);

        # setup gamma according to https://gkyl.readthedocs.io/en/latest/dev/ssp-rk.html
        gamma = zeros(rkStages);
        if rkStages == 1
            gamma[1] = 0.0;
        elseif rkStages == 2
            gamma[1] = 0.0;
            gamma[2] = 1.0;
        elseif rkStages == 3
            gamma[1] = 0.0;
            gamma[2] = 1.0;
            gamma[3] = 0.5;
        end

        # allocate memory for matrices used in Update function
        KRK = zeros(rkStages+1,Nx,settings.r);
        SRK = zeros(rkStages+1,settings.r,settings.r);
        LRK = zeros(rkStages+1,settings.N^2,settings.r);
        Krhs = zeros(rkStages,Nx,settings.r);
        Srhs = zeros(rkStages,settings.r,settings.r);
        Lrhs = zeros(rkStages,settings.N^2,settings.r);

        new(alpha,beta,gamma,rkStages,Nx,dt,KRK,SRK,LRK,Krhs,Srhs,Lrhs)
    end
end

function SetUpRKTable(p::Int64)
    alpha = zeros(p,p);
    beta = zeros(p,p);
    if p == 1
        alpha[1,1] = 1.0;
        beta[1,1] = 1.0;
    elseif p == 2
        alpha[1,1] = 1.0;
        alpha[2,1] = 0.5;
        alpha[2,2] = 0.5;
        beta[1,1] = 1.0;
        beta[2,1] = 0.0;
        beta[2,2] = 0.5;
    elseif p == 3
        alpha[1,1] = 1.0;
        alpha[2,1] = 3.0/4.0;
        alpha[2,2] = 1.0/4.0;
        alpha[3,1] = 1.0/3.0;
        alpha[3,3] = 2.0/3.0;
        beta[1,1] = 1.0;
        beta[2,1] = 0.0;
        beta[2,2] = 1.0/4.0;   
        beta[3,3] = 2.0/3.0; 
    end
    return alpha, beta
end
