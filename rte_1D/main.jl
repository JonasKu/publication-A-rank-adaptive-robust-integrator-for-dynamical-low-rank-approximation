include("settings.jl")
include("DLRSolver.jl")

using DelimitedFiles
using NPZ
using PyPlot

s = Settings();

# run solver for different end times and different epsilons

# end time tEnd = 2
s.epsAdapt = 1e-1;
s.tEnd = 2.0;
solver = DLRSolver(s);
@time tEnd, rankInTime, normInTime, u3 = SolveUnconventionalEfficient(solver);

s.epsAdapt = 5e-2;
solver = DLRSolver(s);
@time tEnd, rankInTime, normInTime, u3F = SolveUnconventionalEfficient(solver);

# end time tEnd = 2.75
s.epsAdapt = 1e-1;
s.tEnd = 2.75;
solver = DLRSolver(s);
@time tEnd, rankInTime, normInTime, u4 = SolveUnconventionalEfficient(solver);

s.epsAdapt = 5e-2;
solver = DLRSolver(s);
@time tEnd, rankInTime, normInTime, u4F = SolveUnconventionalEfficient(solver);

# end time tEnd = 5
s.epsAdapt = 1e-1;
s.tEnd = 5.0;
solver = DLRSolver(s);
@time tEnd, rankInTime, normInTime, u5 = SolveUnconventionalEfficient(solver);

s.epsAdapt = 5e-2;
solver = DLRSolver(s);
@time tEnd, rankInTimeF, normInTime, u5F = SolveUnconventionalEfficient(solver);

## plotting ##

# plot rank in time

fig, ax = subplots(figsize=(15, 12), dpi=100)
ax.plot(rankInTime[1,:],rankInTime[2,:], "k-", linewidth=2, alpha=1.0)
ax.plot([2,2],[0,100], "b--", linewidth=2, alpha=1.0)
ax.plot([2.75,2.75],[0,100], "r:", linewidth=2, alpha=1.0)
ax.plot([5,5],[0,100], "g-.", linewidth=2, alpha=1.0)
ax.set_xlim([rankInTime[1,1],rankInTime[1,end]+0.05])
ax.set_ylim([2.5,14.5])
ax.set_xlabel("time", fontsize=30);
ax.set_ylabel("rank", fontsize=30);
ax.tick_params("both",labelsize=30) 
fig.canvas.draw() # Update the figure

fig, ax = subplots(figsize=(15, 12), dpi=100)
ax.plot(rankInTimeF[1,:],rankInTimeF[2,:], "k-", linewidth=2, alpha=1.0)
ax.plot([2,2],[0,100], "b--", linewidth=2, alpha=1.0)
ax.plot([2.75,2.75],[0,100], "r:", linewidth=2, alpha=1.0)
ax.plot([5,5],[0,100], "g-.", linewidth=2, alpha=1.0)
ax.set_xlim([rankInTimeF[1,1],rankInTimeF[1,end]+0.05])
ax.set_ylim([2.5,maximum(rankInTimeF[2,:])+0.5])
ax.set_xlabel("time", fontsize=30);
ax.set_ylabel("rank", fontsize=30);
ax.tick_params("both",labelsize=30) 
fig.canvas.draw() # Update the figure

# plot solution

## read reference solution
# t = 5
v = readdlm("PlaneSourceRawT5", ',')
uEx5 = zeros(length(v));
for i = 1:length(v)
    if v[i] == ""
        uEx5[i] = 0.0;
    else
        uEx5[i] = Float64(v[i])
    end
end
x5 = collect(range(-5,5,length=(2*length(v)-1)));
uEx5 = [uEx5[end:-1:2];uEx5];

# t = 2.75
v = readdlm("PlaneSourceRawT2-75", ',')
uEx4 = zeros(length(v));
for i = 1:length(v)
    if v[i] == ""
        uEx4[i] = 0.0;
    else
        uEx4[i] = Float64(v[i])
    end
end
x4 = collect(range(-5,5,length=(2*length(v)-1)));
uEx4 = [uEx4[end:-1:2];uEx4];

# t = 2
v = readdlm("PlaneSourceRawT2Fine", ',')
uEx3 = zeros(length(v));
for i = 1:length(v)
    if v[i] == ""
        uEx3[i] = 0.0;
    else
        uEx3[i] = Float64(v[i])
    end
end
x3 = collect(range(-5,5,length=(2*length(v)-1)));
uEx3 = [uEx3[end:-1:2];uEx3];

# start plot
fig, ax = subplots(figsize=(15, 12), dpi=100)
ax.plot(s.xMid,u3[:,1], "b--", linewidth=2, label=L"t=2", alpha=1.0)
ax.plot(x3,uEx3, "b-", linewidth=2, alpha=0.5)
ax.plot(s.xMid,u4[:,1], "r:", linewidth=2, label=L"t=2.75", alpha=1.0)
ax.plot(x4,uEx4, "r-", linewidth=2, alpha=0.5)
ax.plot(s.xMid,u5[:,1], "g-.", linewidth=2, label=L"t=5", alpha=1.0)
ax.plot(x5,uEx5, "g-", linewidth=2, alpha=0.5)
ylabel("scalar flux", fontsize=30)
ax.set_xlim([-5,5])
ax.set_xlabel("x", fontsize=30);
ax.legend(loc="upper right", fontsize=30)
ax.tick_params("both",labelsize=30) 
fig.canvas.draw() # Update the figure
#PyPlot.savefig("results/res.png")

# start plot
fig, ax = subplots(figsize=(15, 12), dpi=100)
ax.plot(s.xMid,u3F[:,1], "b--", linewidth=2, label=L"t=2", alpha=1.0)
ax.plot(x3,uEx3, "b-", linewidth=2, alpha=0.5)
ax.plot(s.xMid,u4F[:,1], "r:", linewidth=2, label=L"t=2.75", alpha=1.0)
ax.plot(x4,uEx4, "r-", linewidth=2, alpha=0.5)
ax.plot(s.xMid,u5F[:,1], "g-.", linewidth=2, label=L"t=5", alpha=1.0)
ax.plot(x5,uEx5, "g-", linewidth=2, alpha=0.5)
ylabel("scalar flux", fontsize=30)
ax.set_xlim([-5,5])
ax.set_xlabel("x", fontsize=30);
ax.legend(loc="upper right", fontsize=30)
ax.tick_params("both",labelsize=30) 
fig.canvas.draw() # Update the figure

println("main finished")
