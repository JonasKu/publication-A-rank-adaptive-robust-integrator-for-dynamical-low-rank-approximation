clc; clear; close all

rng(0); %set seed.

N= 100;
T = 1; 

%% Matrices

global D

D = 0.5*full(gallery('tridiag',N,-1,2,-1));

global V_cos

    DD = -N/2 : N/2-1;
    dx = (2*pi*N^-1); 
    x = dx.*DD;

    V_cos = diag(1-cos(x));
    
global fun
    
    fun = @(X) -1i*(D*X+X*D'+V_cos*X*V_cos');
    

%% Initial Data:
U0 = orth(rand(N,N));
V0 = orth(rand(N,N));

S0 = zeros(N,N);
 
 for i=1:N
       S0(i,i) = 10^-i;
 end

S0 = S0 ./ norm(S0, 'fro');