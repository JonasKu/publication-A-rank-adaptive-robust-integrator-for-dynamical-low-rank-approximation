struct PNSystem
    # system matrices
    Ax::Array{Float64,2};
    Ay::Array{Float64,2};
    Az::Array{Float64,2};

    # Solver settings
    settings::Settings;

    # total number of moments
    nTotalEntries::Int
    # degree of expansion
    N::Int

    # constructor
    function PNSystem(settings) 
        N = settings.nPN;
        nTotalEntries = GlobalIndex( N, N ) + 1;    # total number of entries for sytem matrix
        Ax = zeros(nTotalEntries,nTotalEntries);
        Ay = zeros(nTotalEntries,nTotalEntries);
        Az = zeros(nTotalEntries,nTotalEntries);

        new(Ax,Ay,Az,settings,nTotalEntries,N);
    end
end

function AParam( l::Int, k::Int )
    return sqrt( ( ( l - k + 1 ) * ( l + k + 1 ) ) / ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) );
end

function BParam( l::Int, k::Int ) 
    return sqrt( ( ( l - k ) * ( l + k ) ) / ( ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ) );
end

function CParam( l::Int, k::Int ) 
    return sqrt( ( ( l + k + 1 ) * ( l + k + 2 ) ) / ( ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ) );
end

function DParam( l::Int, k::Int ) 
    return sqrt( ( ( l - k ) * ( l - k - 1 ) ) / ( ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ) );
end

function EParam( l::Int, k::Int ) 
    return sqrt( ( ( l - k + 1 ) * ( l - k + 2 ) ) / ( ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ) );
end

function FParam( l::Int, k::Int ) 
    return sqrt( ( ( l + k ) * ( l + k - 1 ) ) / ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) );
end

function CTilde( l::Int, k::Int ) 
    if k < 0  return 0.0; end
    if k == 0 
        return sqrt( 2 ) * CParam( l, k );
    else
        return CParam( l, k );
    end
end

function DTilde( l::Int, k::Int ) 
    if k < 0  return 0.0; end
    if k == 0 
        return sqrt( 2 ) * DParam( l, k );
    else
        return DParam( l, k );
    end
end

function ETilde( l::Int, k::Int ) 
    if k == 1 
        return sqrt( 2 ) * EParam( l, k );
    else
        return EParam( l, k );
    end
end

function FTilde( l::Int, k::Int ) 
    if k == 1
        return sqrt( 2 ) * FParam( l, k );
    else
        return FParam( l, k );
    end
end

function Sgn( k::Int ) 
    if k >= 0 
        return 1;
    else
        return -1;
    end
end

function GlobalIndex( l::Int, k::Int ) 
    numIndicesPrevLevel  = l * l;    # number of previous indices untill level l-1
    prevIndicesThisLevel = k + l;    # number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
end

function kPlus( k::Int )  return k + Sgn( k ); end

function kMinus( k::Int )  return k - Sgn( k ); end

function unsigned(x::Float64)
    return Int(floor(x))
end

function int(x::Float64)
    return Int(floor(x))
end

function SetupSystemMatrices(obj::PNSystem)
    nTotalEntries = obj.nTotalEntries;    # total number of entries for sytem matrix
    N = obj.N

    # loop over columns of A
    for l = 0:N
        for k=-l:l
            i = GlobalIndex( l, k ) ;

            # flux matrix in direction x
            if k != -1
                j = GlobalIndex( l - 1, kMinus( k ) );
                if j >= 0 && j < nTotalEntries obj.Ax[i+1,j+1] = 0.5 * CTilde( l - 1, abs( k ) - 1 ); end

                j = GlobalIndex( l + 1, kMinus( k ) );
                if j >= 0 && j < nTotalEntries obj.Ax[i+1,j+1] = -0.5 * DTilde( l + 1, abs( k ) - 1 ); end
            end

            j = GlobalIndex( l - 1, kPlus( k ) );
            if  j >= 0 && j < nTotalEntries  obj.Ax[i+1,j+1] = -0.5 * ETilde( l - 1, abs( k ) + 1 ); end

            j = GlobalIndex( l + 1, kPlus( k ) );
            if  j >= 0 && j < nTotalEntries  obj.Ax[i+1,j+1] = 0.5 * FTilde( l + 1, abs( k ) + 1 ); end

            # flux matrix in direction y
            if  k != 1
                j = GlobalIndex( l - 1, -kMinus( k ) );
                if  j >= 0 && j < nTotalEntries  obj.Ay[i+1,j+1] = -0.5 * Sgn( k ) * CTilde( l - 1, abs( k ) - 1 ); end

                j = GlobalIndex( l + 1, -kMinus( k ) );
                if  j >= 0 && j < nTotalEntries  obj.Ay[i+1,j+1] = 0.5 * Sgn( k ) * DTilde( l + 1, abs( k ) - 1 ); end
            end

            j = GlobalIndex( l - 1, -kPlus( k ) );
            if  j >= 0 && j < nTotalEntries  obj.Ay[i+1,j+1] = -0.5 * Sgn( k ) * ETilde( l - 1, abs( k ) + 1 ); end

            j = GlobalIndex( l + 1, -kPlus( k ) );
            if  j >= 0 && j < nTotalEntries  obj.Ay[i+1,j+1] = 0.5 * Sgn( k ) * FTilde( l + 1, abs( k ) + 1 ); end

            # flux matrix in direction z
            j = GlobalIndex( l - 1, k );
            if  j >= 0 && j < nTotalEntries  obj.Az[i+1,j+1] = AParam( l - 1, k ); end

            j = GlobalIndex( l + 1, k );
            if  j >= 0 && j < nTotalEntries  obj.Az[i+1,j+1] = BParam( l + 1, k ); end
        end
    end
    #println("System Matrix Set UP!")
    #println("A_x =", obj.Ax)
    #println("A_y =", obj.Ay)
    #println("A_z =", obj.Az)
end

#function G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
#    unused( nUnit );
#
#    // return F( 0.5 * ( u + v ) ) * n - 0.5 * ( v - u ) * norm( n );
#    return 0.5 * ( F( u ) + F( v ) ) * n - 0.5 * _AbsAx * ( v - u ) * fabs( n[0] ) - 0.5 * _AbsAz * ( v - u ) * fabs( n[1] );
#end
