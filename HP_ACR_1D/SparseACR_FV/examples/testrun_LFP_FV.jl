using Pkg

working_dir = pwd()
Pkg.activate(working_dir)

# Pkg.activate("C:\\Users\\Sam\\RnD\\2DCH_SBM.jl\\2DCH_SBM\\")
# include("C:\\Users\\Sam\\RnD\\2DCH_SBM.jl\\2DCH_SBM\\src\\2DCH_SBM.jl")
srcfile = "\\src\\sparseACR_FV.jl"
include(working_dir*srcfile)

Nx =16384;

Lx = 1e-7 #m
x = LinRange(0.0, Lx, Nx)
dx = x[2] - x[1];


# Set up Material Parameters
c_max = 2.29e4 #mol/m^3
Ω_enthalpy = 0.115; #eV\
T = 298.15; #K
Vθ = 3.42; #V
#κ = 3.13e9; #eV/m
κ = 3.1314026858813996e9; #eV/m
# k0 = 1.6e-4 # A/m^2 zeng et al. 2014 SIAM
k0 =1.6e-4
# k0 =1.0;
# k0 = 0.1
α = 0.5;

D = 1e-14


# C-rate 
C_rate = 1e-6; #1/hr

mp = MaterialDomainParam(Lx,c_max,Ω_enthalpy,T,Vθ,κ,k0,α,C_rate,Nx,dx);
 
sp,sarg, sn, sc_rhs = set_up_sim(mp)



u0 = sarg.c0[:]
du0 = similar(u0);
p=[]
# ACR_RHS_vc=build_functions_vc(sp,sarg, sn, sc_rhs,Vfunc)

# println("RHS eval")
# @time ACR_RHS_vc(du0,u0,p,0.0);
# println("")

# @show sp
# CHR_RHS = build_functions_cc(sp,sarg, sn, sc_rhs)
# t0,tf_tilde = sarg.tspan
# ts = range(t0,tf_tilde;length=100)
# p = zeros(sp.Np)


# u0 = [sarg.c0[:] ; sp.Vθ];
# du0 = similar(u0);

# @warn "Model current is not parameterized through DifferentialEquations.jl parameter so p is blank"
# p =[]
# println("RHS eval")
# @time CHR_RHS(du0,u0,p,0.0);
# println("")

# dV = du0[end]

# prob = build_forward_solve_cc(u0,sarg,p,CHR_RHS)
T = 298.15; #K
kB = 8.617333262145e-5; #Boltzmann Constant [eV/K]
η=(-2.0*kB*T)
Vfunc(t)=(η+Vθ)/(kB*T)


sp,sarg,sn,sol =runSim_vc(mp,Vfunc);
# plot_time_step_hist(sol.t;Nbins=100)
# sol
xgrid = LinRange(0.0, 1.0, Nx);
# cflag=true
vflag=false
plot_anim_sol(xgrid, sol;cc_flag=vflag,title="Presentation_fast_lith")


# sp,sarg,sn,sol =runSim_cc(mp);

# plotres(xgrid,sol.u[end][1:end-1])

# tf = sol.t[end];
# ts = 0.0:tf/100:tf

