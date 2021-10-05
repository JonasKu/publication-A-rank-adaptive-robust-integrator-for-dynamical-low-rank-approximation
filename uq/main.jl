include("settings.jl")
include("DLRSolver.jl")

using NPZ;
using PyPlot

s = Settings();

############################ run solver for fixed theta and maximal/minimal rank ############################

rMax = s.r;
s.epsAdapt = 1.5e-2
solver = DLRSolver(s);

@time tEnd,rankInTime, X,S,W = SolveForwardAdaptive(solver);

uDLRadapt = Array((X*S*W')');

############################
s.r = minimum(rankInTime[2,:]);
solver = DLRSolver(s);

@time tEnd, X,S,W = SolveForward(solver);

uDLRMin = Array((X*S*W')');

############################
s.r = maximum(rankInTime[2,2:end]);
solver = DLRSolver(s);

@time tEnd, X,S,W = SolveForward(solver);

uDLRMax = Array((X*S*W')');

############################

# run solver for different shoices of theta
# run solver for thetaBar = 1.2e-2
s.r = rMax;
s.epsAdapt = 1.2e-2;
solver = DLRSolver(s);

@time tEnd,rankInTimeF, X,S,W = SolveForwardAdaptive(solver);

uDLRadaptF = Array((X*S*W')');

# run solver for thetaBar = 1e-2
s.r = rMax;
s.epsAdapt = 1.0e-2;
solver = DLRSolver(s);

@time tEnd,rankInTimeFF, X,S,W = SolveForwardAdaptive(solver);

uDLRadaptFF = Array((X*S*W')');

##########################################################################################################################
######################################################## Plotting ########################################################
##########################################################################################################################

q = Quadrature(s.Nq,"Gauss");
basis = Basis(q,s);
qFine = Quadrature(100,"Gauss")

s.r = rankInTime[2,end]
b = Basis(q,s);
s.r = minimum(rankInTime[2,:]);
bMin = Basis(q,s);
s.r = maximum(rankInTime[2,:]);
bMax = Basis(q,s);

Nq = s.Nq;
Nx = s.Nx;
NxFine = 1000;
xFine = collect(range(s.a,s.b,length=NxFine))
uExact = zeros(NxFine);
varExact = zeros(NxFine);
uPlot = zeros(Nx);
varPlot = zeros(Nx);
uFPlot = zeros(Nx);
varFPlot = zeros(Nx);
uFFPlot = zeros(Nx);
varFFPlot = zeros(Nx);
vPlot = zeros(Nx);
varVPlot = zeros(Nx);
wPlot = zeros(Nx);
varWPlot = zeros(Nx);

# start plot
fig, ax = subplots(figsize=(15, 12), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung

# compute expected value and variance
for j = 1:Nx
    uVals = EvalAtQuad(b,uDLRadapt[:,j]);
    uPlot[j] = Integral(q,uVals*0.25);
    varPlot[j] = Integral(q,0.25*(uVals.-uPlot[j]).^2);
    vVals = EvalAtQuad(bMin,uDLRMin[:,j]);
    vPlot[j] = Integral(q,vVals*0.25);
    varVPlot[j] = Integral(q,0.25*(vVals.-vPlot[j]).^2);
    wVals = EvalAtQuad(bMax,uDLRMax[:,j]);
    wPlot[j] = Integral(q,wVals*0.25);
    varWPlot[j] = Integral(q,0.25*(wVals.-wPlot[j]).^2);
    uVals = EvalAtQuad(b,uDLRadaptF[:,j]);
    uFPlot[j] = Integral(q,uVals*0.25);
    varFPlot[j] = Integral(q,0.25*(uVals.-uPlot[j]).^2);
    uVals = EvalAtQuad(b,uDLRadaptFF[:,j]);
    uFFPlot[j] = Integral(q,uVals*0.25);
    varFFPlot[j] = Integral(q,0.25*(uVals.-uPlot[j]).^2);
end
varMax = maximum(varPlot);
expMax = maximum(uPlot);

exactState = zeros(NxFine,qFine.Nq,qFine.Nq);
for j = 1:NxFine
    for k = 1:qFine.Nq
        for l = 1:qFine.Nq
            exactState[j,k,l] = s.solutionExact(s.tEnd,xFine[j],qFine.xi[k],qFine.xi[l])[1];
        end
    end
end
for j = 1:NxFine
    for k = 1:qFine.Nq
        for l = 1:qFine.Nq
            uExact[j] += exactState[j,k,l]*0.25*qFine.w[k]*qFine.w[l];
        end
    end
    for k = 1:qFine.Nq
        for l = 1:qFine.Nq
            varExact[j] += (exactState[j,k,l]-uExact[j])^2 * 0.25*qFine.w[k]*qFine.w[l];
        end
    end
end
x = s.x

ax.plot(x,uPlot, "k--", linewidth=2, label="DLR adaptive", alpha=1.0)
ax.plot(x,vPlot, "g:", linewidth=2, label=L"DLR$_{9}$", alpha=1.0)
ax.plot(x,wPlot, "m-.", linewidth=2, label=L"DLR$_{25}$", alpha=1.0)
ylabel("Expectation", fontsize=30,color="red")
ax.plot(xFine,uExact, "r-", linewidth=2, alpha=0.5)
ax2 = ax[:twinx]() # Create another axis on top of the current axis
ylabel("Standard deviation", fontsize=30,color="blue")
ax2.plot(x,sqrt.(varPlot), "k--", linewidth=2, label="SG", alpha=1.0)
ax2.plot(x,sqrt.(varVPlot), "g:", linewidth=2, label="DLR", alpha=1.0)
ax2.plot(x,sqrt.(varWPlot), "m-.", linewidth=2, label="unconventional DLR", alpha=1.0)
#ax2[:set_position](new_position) # Position Method 2
setp(ax2[:get_yticklabels](),color="blue") # Y Axis font formatting
setp(ax[:get_yticklabels](),color="red")
ax2.plot(xFine,sqrt.(varExact), "b-", linewidth=2, alpha=0.5)
#ylimMinus = -0.5;
#ylimPlus = 16.0
#ax[:set_ylim]([ylimMinus,ylimPlus])
ax2.set_ylim([-0.2,6.3])
ax.set_xlim([s.a,s.b])
ax.set_xlabel("x", fontsize=30);
ax.legend(loc="upper right", fontsize=30)
ax.tick_params("both",labelsize=30) 
ax2.tick_params("both",labelsize=30)
fig.canvas.draw() # Update the figure

# start plot
fontsizeVal=20;
fig, ax = subplots(figsize=(15, 9), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax.plot(x,uPlot, "k--", linewidth=2, label=L"\bar{\vartheta} = 0.015", alpha=1.0)
ax.plot(x,uFPlot, "g-.", linewidth=2, label=L"\bar{\vartheta} = 0.012", alpha=1.0)
ax.plot(x,uFFPlot, "m:", linewidth=2, label=L"\bar{\vartheta} = 0.01", alpha=1.0)
ylabel("Expectation", fontsize=fontsizeVal,color="red")
ax.plot(xFine,uExact, "r-", linewidth=2, alpha=0.5)
ax2 = ax[:twinx]() # Create another axis on top of the current axis
ylabel("Standard deviation", fontsize=fontsizeVal,color="blue")
ax2.plot(x,sqrt.(varPlot), "k--", linewidth=2, label="SG", alpha=1.0)
ax2.plot(x,sqrt.(varFPlot), "g-.", linewidth=2, label="DLR", alpha=1.0)
ax2.plot(x,sqrt.(varFFPlot), "m:", linewidth=2, label="DLR", alpha=1.0)
setp(ax2[:get_yticklabels](),color="blue") # Y Axis font formatting
setp(ax[:get_yticklabels](),color="red")
ax2.plot(xFine,sqrt.(varExact), "b-", linewidth=2, alpha=0.5)
ax2.set_ylim([-0.2,6.3])
ax.set_xlim([s.a,s.b])
ax.set_xlabel("x", fontsize=fontsizeVal);
ax.legend(loc="upper right", fontsize=fontsizeVal)
ax.tick_params("both",labelsize=fontsizeVal) 
ax2.tick_params("both",labelsize=fontsizeVal)
fig.canvas.draw() # Update the figure

# plot rank over time
fig, ax = subplots(figsize=(15, 12), dpi=100)
ax.plot(rankInTime[1,:],rankInTime[2,:], "k--", linewidth=2, label=L"\bar{\vartheta} = 0.015", alpha=1.0)
ax.plot(rankInTimeF[1,:],rankInTimeF[2,:], "g-.", linewidth=2, label=L"\bar{\vartheta} = 0.012", alpha=1.0)
ax.plot(rankInTimeFF[1,:],rankInTimeFF[2,:], "m:", linewidth=2, label=L"\bar{\vartheta} = 0.01", alpha=1.0)
ax.set_xlim([rankInTime[1,1],rankInTime[1,end]])
#ax.set_ylim([8,30])
ax.legend(loc="upper left", fontsize=30)
ax.set_xlabel("time", fontsize=30);
ax.set_ylabel("rank", fontsize=30);
ax.tick_params("both",labelsize=30) 
fig.canvas.draw() # Update the figure
