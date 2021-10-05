using Base: Float64
include("settings.jl")
include("SolverDLRA.jl")

using PyPlot
using DelimitedFiles

s = Settings(251,251,1000); # create settings class with 251 x 251 spatial cells and a maximal rank of 1000

################################################################
######################### execute code #########################
################################################################

################### run explicit Euler ###################
s.epsAdapt = 5e-2;
solver = SolverDLRA(s);
@time u1, rankInTime1 = SolveUnconventionalAdaptive(solver);
u1 = Vec2Mat(s.NCellsX,s.NCellsY,u1)

s.epsAdapt = 2.5e-2;
solver = SolverDLRA(s);
@time u2, rankInTime2 = SolveUnconventionalAdaptive(solver);
u2 = Vec2Mat(s.NCellsX,s.NCellsY,u2)

s.epsAdapt = 1e-2;
solver = SolverDLRA(s);
@time u3, rankInTime3 = SolveUnconventionalAdaptive(solver);
u3 = Vec2Mat(s.NCellsX,s.NCellsY,u3)

################### run Heun ###################
s.epsAdapt = 5e-2;
solver = SolverDLRA(s);
@time u1H, rankInTime1H = SolveUnconventionalAdaptiveHeun(solver);
u1H = Vec2Mat(s.NCellsX,s.NCellsY,u1H)

s.epsAdapt = 2.5e-2;
solver = SolverDLRA(s);
@time u2H, rankInTime2H = SolveUnconventionalAdaptiveHeun(solver);
u2H = Vec2Mat(s.NCellsX,s.NCellsY,u2H)

s.epsAdapt = 1e-2;
solver = SolverDLRA(s);
@time u3H, rankInTime3H = SolveUnconventionalAdaptiveHeun(solver);
u3H = Vec2Mat(s.NCellsX,s.NCellsY,u3H)

############################################################
######################### plotting #########################
############################################################

################### read in reference solution ###################
lsRef = readdlm("exactLineSource.txt", ',', Float64);
xRef = lsRef[:,1];
phiRef = lsRef[:,2];
lsRefFull = readdlm("refPhiFull.txt", ',', Float64);

################### plot reference cut ###################
fig = figure("u cut ref",figsize=(10, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax.plot(xRef,phiRef, "k-", linewidth=2, label="exact", alpha=0.8)
ylabel("scalar flux", fontsize=20)
xlabel(L"$x$", fontsize=20)
ax.set_xlim([s.a,s.b])
ax.tick_params("both",labelsize=20) 
tight_layout()
show()
savefig("scalar_flux_reference_cut.pdf")

################### plot reference full ###################
fig = figure("scalar_flux_reference",figsize=(10,10),dpi=100)
ax = gca()
pcolormesh(s.xMid,s.yMid, lsRefFull,vmin=0.0,vmax=maximum(lsRefFull), cmap=:gray)
#pcolormesh(s.xMid,s.yMid, lsRefFull,vmin=0.0,vmax=maximum(lsRefFull))
ax.tick_params("both",labelsize=20) 
plt.xlabel("x", fontsize=20)
plt.ylabel("y", fontsize=20)
plt.title("scalar flux, reference", fontsize=25)
savefig("scalar_flux_reference.pdf")

################### plot scalar fluxes full ###################

## explicit Euler
fig = figure("scalar_flux_vartheta_0-05_exp_Euler",figsize=(10,10),dpi=100)
ax = gca()
pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u1[:,:,1],vmin=0.0,vmax=maximum(lsRefFull))
ax.tick_params("both",labelsize=20) 
plt.xlabel("x", fontsize=20)
plt.ylabel("y", fontsize=20)
plt.title(L"scalar flux, $\bar{\vartheta} = 0.05$, explicit Euler", fontsize=25)
savefig("scalar_flux_vartheta_0-05_exp_Euler_nx$(s.NCellsX).pdf")

fig = figure("scalar_flux_vartheta_0-025_exp_Euler",figsize=(10,10),dpi=100)
ax = gca()
pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u2[:,:,1],vmin=0.0,vmax=maximum(lsRefFull))
ax.tick_params("both",labelsize=20) 
plt.xlabel("x", fontsize=20)
plt.ylabel("y", fontsize=20)
plt.title(L"scalar flux, $\bar{\vartheta} = 0.025$, explicit Euler", fontsize=25)
savefig("scalar_flux_vartheta_0-025_exp_Euler_nx$(s.NCellsX).pdf")

fig = figure("scalar_flux_vartheta_0-01_exp_Euler",figsize=(10,10),dpi=100)
ax = gca()
pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u3[:,:,1],vmin=0.0,vmax=maximum(lsRefFull))
ax.tick_params("both",labelsize=20) 
plt.xlabel("x", fontsize=20)
plt.ylabel("y", fontsize=20)
plt.title(L"scalar flux, $\bar{\vartheta} = 0.01$, explicit Euler", fontsize=25)
savefig("scalar_flux_vartheta_0-01_exp_Euler_nx$(s.NCellsX).pdf")

## Heun's method
fig = figure("scalar_flux_vartheta_0-05_heun",figsize=(10,10),dpi=100)
ax = gca()
pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u1H[:,:,1],vmin=0.0,vmax=maximum(lsRefFull))
ax.tick_params("both",labelsize=20) 
plt.xlabel("x", fontsize=20)
plt.ylabel("y", fontsize=20)
plt.title(L"scalar flux, $\bar{\vartheta} = 0.05$, Heun's method", fontsize=25)
savefig("scalar_flux_vartheta_0-05_heun_nx$(s.NCellsX).pdf")

fig = figure("scalar_flux_vartheta_0-025_heun",figsize=(10,10),dpi=100)
ax = gca()
pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u2H[:,:,1],vmin=0.0,vmax=maximum(lsRefFull))
ax.tick_params("both",labelsize=20) 
plt.xlabel("x", fontsize=20)
plt.ylabel("y", fontsize=20)
plt.title(L"scalar flux, $\bar{\vartheta} = 0.025$, Heun's method", fontsize=25)
savefig("scalar_flux_vartheta_0-025_heun_nx$(s.NCellsX).pdf")

fig = figure("scalar_flux_vartheta_0-01_heun",figsize=(10,10),dpi=100)
ax = gca()
pcolormesh(s.xMid,s.yMid, 4.0*pi*sqrt(2)*u3H[:,:,1],vmin=0.0,vmax=maximum(lsRefFull))
ax.tick_params("both",labelsize=20) 
plt.xlabel("x", fontsize=20)
plt.ylabel("y", fontsize=20)
plt.title(L"scalar flux, $\bar{\vartheta} = 0.01$, Heun's method", fontsize=25)
savefig("scalar_flux_vartheta_0-01_heun_nx$(s.NCellsX).pdf")

################### plot scalar fluxes cut ###################

fig = figure("u cut Euler",figsize=(10, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax.plot(xRef,phiRef, "k-", linewidth=2, label="exact", alpha=0.8)
ax.plot(s.yMid,4.0*pi*sqrt(2)*u1[Int(floor(s.NCellsX/2+1)),:,1], "b--", linewidth=2, label=L"DLRA, $\bar{\vartheta} = 0.05$", alpha=0.8)
ax.plot(s.yMid,4.0*pi*sqrt(2)*u2[Int(floor(s.NCellsX/2+1)),:,1], "g-.", linewidth=2, label=L"DLRA, $\bar{\vartheta} = 0.025$", alpha=0.8)
ax.plot(s.yMid,4.0*pi*sqrt(2)*u3[Int(floor(s.NCellsX/2+1)),:,1], "m:", linewidth=2, label=L"DLRA, $\bar{\vartheta} = 0.01$", alpha=1.0)
ax.legend(loc="upper left", fontsize=20)
ylabel("scalar flux", fontsize=20)
xlabel(L"$x$", fontsize=20)
ax.set_xlim([s.a,s.b])
ax.tick_params("both",labelsize=20) 
tight_layout()
show()
savefig("u_cut_exp_Euler_nx$(s.NCellsX).pdf")

fig = figure("u cut heun",figsize=(10, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax.plot(xRef,phiRef, "k-", linewidth=2, label="exact", alpha=0.8)
ax.plot(s.yMid,4.0*pi*sqrt(2)*u1H[Int(floor(s.NCellsX/2+1)),:,1], "b--", linewidth=2, label=L"DLRA, $\bar{\vartheta} = 0.05$", alpha=0.8)
ax.plot(s.yMid,4.0*pi*sqrt(2)*u2H[Int(floor(s.NCellsX/2+1)),:,1], "g-.", linewidth=2, label=L"DLRA, $\bar{\vartheta} = 0.025$", alpha=0.8)
ax.plot(s.yMid,4.0*pi*sqrt(2)*u3H[Int(floor(s.NCellsX/2+1)),:,1], "m:", linewidth=2, label=L"DLRA, $\bar{\vartheta} = 0.01$", alpha=1.0)
ax.legend(loc="upper left", fontsize=20)
ylabel("scalar flux", fontsize=20)
xlabel(L"$x$", fontsize=20)
ax.set_xlim([s.a,s.b])
ax.tick_params("both",labelsize=20) 
tight_layout()
show()
savefig("u_cut_heun_nx$(s.NCellsX).pdf")

################### plot rank in time ###################
# explicit Euler
fig = figure("rank in time Euler",figsize=(10, 10), dpi=100)
ax = gca()
ax.plot(rankInTime1[1,:],rankInTime1[2,:], "b--", linewidth=2, label=L"$\bar{\vartheta} = 0.05$", alpha=1.0)
ax.plot(rankInTime2[1,:],rankInTime2[2,:], "g-.", linewidth=2, label=L"$\bar{\vartheta} = 0.025$", alpha=1.0)
ax.plot(rankInTime3[1,:],rankInTime3[2,:], "m:", linewidth=2, label=L"$\bar{\vartheta} = 0.01$", alpha=1.0)
ax.set_xlim([0.0,s.tEnd])
ax.set_xlabel("time", fontsize=20);
ax.set_ylabel("rank", fontsize=20);
ax.tick_params("both",labelsize=20) 
ax.legend(loc="upper left", fontsize=20)
tight_layout()
fig.canvas.draw() # Update the figure
savefig("rank_in_time_euler_nx$(s.NCellsX).pdf")

# Heun's method
fig = figure("rank in time Heun",figsize=(10, 10), dpi=100)
ax = gca()
ax.plot(rankInTime1H[1,:],rankInTime1H[2,:], "b--", linewidth=2, label=L"$\bar{\vartheta} = 0.05$", alpha=1.0)
ax.plot(rankInTime2H[1,:],rankInTime2H[2,:], "g-.", linewidth=2, label=L"$\bar{\vartheta} = 0.025$", alpha=1.0)
ax.plot(rankInTime3H[1,:],rankInTime3H[2,:], "m:", linewidth=2, label=L"$\bar{\vartheta} = 0.01$", alpha=1.0)
ax.set_xlim([0.0,s.tEnd])
ax.set_xlabel("time", fontsize=20);
ax.set_ylabel("rank", fontsize=20);
ax.tick_params("both",labelsize=20) 
ax.legend(loc="upper left", fontsize=20)
tight_layout()
fig.canvas.draw() # Update the figure
savefig("rank_in_time_heun_nx$(s.NCellsX).pdf")

println("main finished")
