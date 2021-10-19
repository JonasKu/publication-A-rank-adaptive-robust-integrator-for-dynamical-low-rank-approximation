%% Cleaning and set parameters
run Init

%% Reference solution:

sol = Ref(T, U0, V0, S0);

% hold off
% 
% subplot(1,2,1);
% v = svd(sol);
% v = v(1:12);
% semilogy(v);
% pause(0.01);
% 
% hold on

%% Mat proj split VS Adaptive unconventional:
    
norm_ref = norm(S0, 'fro');

tmp1 = Build({U0,V0,S0});
tmp2 = 1i*fun(tmp1);
ener_ref = abs(tmp1(:)'*tmp2(:));

tt = [10^-1 5*10^-2 2*10^-2 10^-2 5*10^-3 2*10^-3 10^-3];

index_h=1;
for h = tt
    
    index_r=1;
    for r=[4 8 12] %[4 6 8 10 12] 

        UU0 = U0(:,1:r);
        VV0 = V0(:,1:r);
        SS0 = S0(1:r, 1:r);
        
        Y0 ={UU0, VV0, SS0};
        
        Y1 = Y0;
        Y1_std = Y0; 
        %Y1_PSI = Y0;
        
        Max = T/h;
        for i=1:Max
            clc;
            fprintf('r = %d, t = %f \n', rank(Y0{end}), i*h);
            
            Y1 = Method_adaptive(Y1, (i-1)*h, i*h);
            rr(i, index_h) = rank(Y1{end});
                        
            Y1_std = Method_standard(Y1_std, (i-1)*h, i*h);
            
            %Y1_PSI = PSI(Y1_PSI, (i-1)*h, i*h);            
            
            % Adaptive
            val_norm(i,index_h) = abs(norm(Y1{end},'fro'));
            
            tmp1 = Build(Y1);
            tmp2 = 1i*fun(tmp1);
            val_energy(i,index_h) = abs(tmp1(:)'*tmp2(:));
            
            % Unconventional:
            val_norm_std(i,index_h) = abs(norm(Y1_std{end},'fro'));
            
            tmp1 = Build(Y1_std);
            tmp2 = 1i*fun(tmp1);
            val_energy_std(i,index_h) = abs(tmp1(:)'*tmp2(:));

            
            % Projector-splitting integrator:
            %val_norm_psi(i,index_h) = abs(norm(Y1_PSI{end},'fro'));
            
            %tmp1 = Build(Y1_PSI);
            %tmp2 = 1i*fun(tmp1);
            %val_energy_psi(i,index_h) = abs(tmp1(:)'*tmp2(:));
            
        end 
        
        approxSol{index_r, index_h} = Y1;

        Error(index_r,index_h) = norm(sol - Build(Y1), 'fro');
        
        Error_std(index_r,index_h) = norm(sol - Build(Y1_std), 'fro');

        
        index_r = index_r+1;
    end
    
    t(index_h) = h;
    index_h=index_h+1;
end

% Unconventional with rank 30:
        r = 30;
        h = 10^-3;

        UU0 = U0(:,1:r);
        VV0 = V0(:,1:r);
        SS0 = S0(1:r, 1:r);
        
        Y0 ={UU0, VV0, SS0};
        
        Y1_std = Y0; 
        
        Max = T/h;
        for i=1:Max
            clc;
            fprintf('r = %d, t = %f \n', rank(Y0{end}), i*h);
            
            Y1_std = Method_standard(Y1_std, (i-1)*h, i*h);

            
            % Unconventional:
            val_norm_std_hr(i,1) = abs(norm(Y1_std{end},'fro'));
            
            tmp1 = Build(Y1_std);
            tmp2 = 1i*fun(tmp1);
            val_energy_std_hr(i,1) = abs(tmp1(:)'*tmp2(:));
        end
            

%% Plotting:
    ymin = 1e-1*min( [min(min(Error)) min(min(Error_std)) ] );

% Standard method:
    subplot(1,3,1)
    loglog_conv(t, Error_std);

    ylim([ymin 1e-1]);

    title('Unconventional')
    xlabel('time step-size')
    ylabel('Error');

    % Slope-1
    hold on
    loglog(t,t)
    
% Error plot adaptive method:
    subplot(1,3,2)
    loglog_conv(t, Error);
    
    ylim([ymin 1e-1]);

    title('Adaptive');
    xlabel('time step-size')
    ylabel('Error');

    %Slope-1
    hold on
    loglog(t,t)
    loglog(t,t.^2, '--')
    
% Rank-plot
    subplot(1,3,3)
    for i=1:3:size(rr,2)
        
        y = rr(:,i);
        y = nonzeros(y);
        
        dt = T / length(y);
        x = (1:length(y))*dt;
        
        plot(x,y, 'DisplayName','h='+string(dt));
        
        hold on
    end
    xlabel('time');
    ylabel('Rank');
    ylim([min(min(rr)) max(rr(:))+2]);
    
    legend()
    
    saveas(gcf,'./png/ErrorPlot.fig')
    saveas(gcf,'./png/ErrorPlot.png')

    % Plot difference with Energy and Norm:
    index = find(tt == 10^-3);
    
    norm_adapt = val_norm(:,index);
    energy_adapt = val_energy(:,index);
    
    norm_uncov = val_norm_std(:,index);
    energy_unconv = val_energy_std(:,index);

    hold on
    figure
    subplot(1,2,1)
    %set(gca,'Ydir','reverse')
    semilogy(x, abs(norm_uncov - norm_ref), '--', 'DisplayName', 'Unconventional')
    
    hold on
    semilogy(x, abs(val_norm_std_hr - norm_ref), ':', 'DisplayName', 'Unconventional - rank 30')
    
    hold on
    semilogy(x, abs(norm_adapt - norm_ref), 'DisplayName', 'Adaptive')

    title('Norm')
    legend()
    
    hold on
    subplot(1,2,2)
    
    semilogy(x, abs(energy_unconv - ener_ref),'--', 'DisplayName', 'Unconventional')
    
    hold on
    semilogy(x, abs(val_energy_std_hr - ener_ref), ':', 'DisplayName', 'Unconventional - rank 30')
    
    hold on
    semilogy(x, abs(energy_adapt - ener_ref), 'DisplayName', 'Adaptive')
    
    title('Energy')
    legend()
    
    
    saveas(gcf,'./png/Conservation.png')
    saveas(gcf,'./png/Conservation.fig')
    
%% Methods and functions:    
function Y1 = Method_adaptive(Y0, t0, t1)
    
    %global D
    %global V_cos
    global fun
    
    U_0 = Y0{1};
    V_0 = Y0{2};
    S_0 = Y0{3};
    
    h = t1-t0;
    
    % K1-step:
    funK = @(K) fun(K*V_0')*V_0;
    
    K_0 = U_0*S_0;
    K_1 = rk4(funK, K_0, h);
    
    K_1 = [K_1, U_0];
    [U_1, ~] = qr(K_1,0);
    
    % K2-step:
    funK = @(K) fun(U_0*K')'*U_0;
    
    K_0 = V_0*S_0';
    K_1 = rk4(funK, K_0, h);
    
    K_1 = [K_1, V_0];
    [V_1, ~] = qr(K_1,0);
    
    
    % S-step:
    funS = @(S) U_1'*fun(U_1*S*V_1')*V_1;
    
    S_0 = (U_1'*U_0) * S_0 * (V_1'*V_0)';
    S_1 = rk4(funS, S_0, h);
    
    % Compute singular values of S1 and decide how to truncate:
    [U,S,V] = svd(S_1);
    
    tol = 1e-6;
    %tol = tol*norm(S);
    
    rmax = size(S_1,1)/2;
    sg = diag(S); %sum(sg)
    
    for j=1:2*rmax
        tmp = sqrt(sum(sg(j:2*rmax)).^2);
        if(tmp<tol)
            break;
        end
    end
    
    rmax = j;
    rmax = min(rmax,30);
    
    % To use in the rank-fixed way - just comment the previous part.
    
    %Truncation:
    U_1 = U_1*U;
    V_1 = V_1*V;
    
    S_1 = S(1:rmax,1:rmax);
    U_1 = U_1(:,1:rmax);
    V_1 = V_1(:,1:rmax);
    
    Y1 = {U_1, V_1, S_1};    

end

function Y1 = PSI(Y0, t0, t1)

 global fun

    U_0 = Y0{1};
    V_0 = Y0{2};
    S_0 = Y0{3};
    
    h = t1-t0;
    
    % K1-step:
    funK = @(K) fun(K*V_0')*V_0;
    
    K_0 = U_0*S_0;
    K_1 = rk4(funK, K_0, h);
    
    [U_1, S_0] = qr(K_1,0);
    
    % S-step:
    funS = @(S) -U_1'*fun(U_1*S*V_0')*V_0;
    S_1 = rk4(funS, S_0, h);

    
    % K2-step:
    funK = @(K) fun(U_1*K')'*U_1;
    
    K_0 = V_0*S_1';
    K_1 = rk4(funK, K_0, h);
    
    [V_1, S_1 ] = qr(K_1,0);    
    S_1 = S_1';
    
    Y1 = {U_1, V_1, S_1};    

end

function Y1 = Method_standard(Y0, t0, t1)
    
    global fun

    U_0 = Y0{1};
    V_0 = Y0{2};
    S_0 = Y0{3};
    
    h = t1-t0;
    
    % K1-step:
    funK = @(K) fun(K*V_0')*V_0;
    
    K_0 = U_0*S_0;
    K_1 = rk4(funK, K_0, h);
    
    [U_1, ~] = qr(K_1,0);
    
    % K2-step:
    funK = @(K) fun(U_0*K')'*U_0;
    
    K_0 = V_0*S_0';
    K_1 = rk4(funK, K_0, h);
    
    [V_1, ~ ] = qr(K_1,0);
    
    
    % S-step:
    funS = @(S) U_1'*fun(U_1*S*V_1')*V_1;
    
    S_0 = (U_1'*U_0) * S_0 * (V_1'*V_0)';
    S_1 = rk4(funS, S_0, h);
    
    Y1 = {U_1, V_1, S_1};    

end

function Y = rk2(f,Y, h)
     k1 = h*f(Y);
     k2 = h*f(Y+h*k1);

 
     Y = Y + 0.5*(k1+k2);   

end

function Y = rk4(f,Y, h)
     k1 = h*f(Y);
     k2 = h*f(Y + 0.5*k1);
     k3 = h*f(Y + 0.5*k2);
     k4 = h*f(Y + k3);
 
     Y = Y + (k1 + 2*k2 + 2*k3 + k4) / 6;   

end

function y = Build(X)
    U = X{1};
    V = X{2};
    S = X{3};
    
    y = U*S*V';
end

function loglog_conv(vect,err)

    [n1 n2]=size(err);

    symbols='sox.d+^*v><';
    ms=[6 6 6 8 6 6 6 6 6 6 6];
     gr=(linspace(.66,0,n1))';
     colors=[gr gr gr];

    for jj=1:n1
        loglog(vect,err(jj,:), ...
               'LineWidth',1,...
               'Marker',symbols(jj),...
               'MarkerSize',ms(jj),...
               'MarkerIndices', 1:1:length(vect), ...
               'MarkerSize',12,...
               'Color', colors(jj,:));
        if jj==1
            hold on;
        end
    end
    hold off;
end

