function y = Ref(t,U0, V0, S0)
    global D
    global V_cos
    global fun
    
    Y0 = U0*S0*V0';
    y = myDopri(fun, Y0, 0, t);
end

function sol = myDopri(fun, y0, t0, t1)
    N = size(y0);

    odefun = @(t,y) F(t,y, fun, N);
    
    tspan = [t0 t1];
    
    param = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

    [a,b] = ode45(odefun, tspan, y0(:), param);    
    
    b = reshape(b.', N(1), N(2),  []);

    sol = b(:,:,end);
end


function dy = F(t,y, fun, N)

    
    X = reshape(y, N);
    tmp = fun(X);
    dy = tmp(:);

end