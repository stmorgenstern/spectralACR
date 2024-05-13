using Pkg
Pkg.activate(".")

working_dir = pwd()
srcfile = "\\src\\spectralCHR.jl"
include(working_dir*srcfile)

using ApproxFun,DifferentialEquations,LinearSolve,LinearAlgebra,SparseArrays,BandedMatrices
using CairoMakie

mutable struct sim_numer_S
    #Forward Transforms
    US0_f
    US1_f
    US2_f
    US3_f
    US4_f

    #Inverse Transform
    US0_i
    US1_i
    US2_i
    US3_i
    US4_i
    
    #Differential operators
    D1
    D2
    D3
    D4

    #Evaluation operators
    ChebEval_l
    ChebEval_r
    
    DChebEval_l
    DChebEval_r

    D2ChebEval_l
    D2ChebEval_r

    D3ChebEval_l
    D3ChebEval_r

end

mutable struct sim_cache_rhs_S
    dchat
    d2chat
    d3chat
    d4chat

    c
    dc
    d2c
    d3c
    d4c

    T1
    T2
    T3
    T4

    c_l
    c_r
    μ_l
    μ_r

    RHS
    

end

function init_cache_S(N;chunk_size=25)
    
    dchat = zeros(N); dchat_c = DiffCache(dchat,chunk_size)
    d2chat = zeros(N); d2chat_c = DiffCache(d2chat,chunk_size)
    d3chat= zeros(N); d3chat_c = DiffCache(d3chat,chunk_size)
    d4chat= zeros(N); d4chat_c = DiffCache(d4chat,chunk_size)

    c = zeros(N); c_c = DiffCache(c,chunk_size)
    dc = zeros(N); dc_c = DiffCache(dc,chunk_size)
    d2c = zeros(N); d2c_c = DiffCache(d2c,chunk_size)
    d3c = zeros(N); d3c_c = DiffCache(d3c,chunk_size)
    d4c = zeros(N); d4c_c = DiffCache(d4c,chunk_size)

    T1= zeros(N); T1_c = DiffCache(T1,chunk_size)
    T2= zeros(N); T2_c = DiffCache(T2,chunk_size)
    T3= zeros(N); T3_c = DiffCache(T3,chunk_size)
    T4= zeros(N); T4_c = DiffCache(T4,chunk_size)

    c_l = zeros(N); c_l_c = DiffCache(c_l,chunk_size)
    c_r = zeros(N); c_r_c = DiffCache(c_r,chunk_size)
    μ_l = zeros(N); μ_l_c = DiffCache(μ_l,chunk_size)
    μ_r = zeros(N); μ_r_c = DiffCache(μ_r,chunk_size)

    RHS = zeros(N); RHS_c = DiffCache(RHS,chunk_size)

    sim_c_rhs_S = sim_cache_rhs_S(dchat_c,d2chat_c,d3chat_c,d4chat_c,
    c_c,dc_c,d2c_c,d3c_c,d4c_c,
    T1_c,T2_c,T3_c,T4_c,
    c_l_c,c_r_c,μ_l_c,μ_r_c,RHS)
    
    return sim_c_rhs_S
end

function build_diff_mat(N,λ)
    off_diag2 = [λ+i for i=0:N-λ-1]
    return (2^(λ-1))*factorial(λ-1)*spdiagm(λ=>off_diag2)
end
function build_BC_eval(Nmodes,p)
    return [((-1)^(n+p))*BC_helper(n,p) for n=0:Nmodes-1]
end
function BC_helper(n,p)
    prod = 1;
    for k=0:(p-1)
        prod = prod*(n^2-k^2)/(2*k+1)
    end
    return prod
end

function build_sim_numerics_S(N)
    C= Chebyshev()
    US1 = Ultraspherical(1)
    US2 = Ultraspherical(2)
    US3 = Ultraspherical(3)
    US4 = Ultraspherical(4)

    US0_f = ApproxFun.plan_transform(C, N);
    US1_f = ApproxFun.plan_transform(US1, N);
    US2_f = ApproxFun.plan_transform(US2, N);
    US3_f = ApproxFun.plan_transform(US3, N);
    US4_f = ApproxFun.plan_transform(US4, N);

    US0_i = ApproxFun.plan_itransform(C, N);
    US1_i = ApproxFun.plan_itransform(US1, N);
    US2_i = ApproxFun.plan_itransform(US2, N);
    US3_i = ApproxFun.plan_itransform(US3, N);
    US4_i = ApproxFun.plan_itransform(US4, N);

    D1 = build_diff_mat(N,1)
    D2 = build_diff_mat(N,2)
    D3 = build_diff_mat(N,3)
    D4 = build_diff_mat(N,4)


    ChebEval_l = build_BC_eval(N,0)
    ChebEval_r = abs.(build_BC_eval(N,0))

    DChebEval_l = build_BC_eval(N,1)
    DChebEval_r = abs.(build_BC_eval(N,1))

    D2ChebEval_l = build_BC_eval(N,2)
    D2ChebEval_r = abs.(build_BC_eval(N,2))

    D3ChebEval_l = build_BC_eval(N,3)
    D3ChebEval_r = abs.(build_BC_eval(N,3))
    
    sn = sim_numer_S(US0_f,US1_f,US2_f,US3_f,US4_f,
                    US0_i,US1_i,US2_i,US3_i,US4_i,
                    D1,D2,D3,D4,
                    ChebEval_l,ChebEval_r,
                    DChebEval_l,DChebEval_r,
                    D2ChebEval_l,D2ChebEval_r,
                    D3ChebEval_l,D3ChebEval_r)

    return sn
end

@inline function M(c)
    return c*(1-c)
end
@inline function dM(c)
    return 1-2*c
end
@inline function d2M(c)
    return -2
end
@inline function μh(c,Ω)
    return log(c/(1-c)) + Ω*(1-2*c)
end
@inline function dμh(c,Ω)
    return 1/(c - c^2) -2*Ω
end
@inline function d2μh(c,Ω)
    return (2*c-1)/(c^2*(1-c)^2)
end

@inline function Term1(c,d4c,κ)
    return -κ*M(c)*d4c
end
@inline function Term2(c,dc,d3c,κ)
    return -κ*dM(c)*dc*d3c
end
@inline function Term3(c,d2c,Ω)
    return M(c)*dμh(c,Ω)*d2c
end
@inline function Term4(c,dc,Ω)
    return (dM(c)*dμh(c,Ω) + M(c)*d2μh(c,Ω))*dc^2
end

@inline function j0(c,k0)
    # return k0*c*sqrt((1.0-c))
    return k0*sqrt(c)*(1.0-c)
end

@inline function BVrxn(c,μ,V,k0,α,Vθ)
    return j0(c,μ,k0)*(exp(-α*(μ + V - Vθ)) - exp((1-α)*(μ + V - Vθ)))
end

function CHR_US_vc(du,u,p,t,sp,sn,sc,Vfunc)
    chat =view(u,1:sp.Nx)
    dchat_t = view(du,1:sp.Nx)
    # dchat_t=du
    V = Vfunc(t)

    dchat = get_tmp(sc.dchat,u)
    d2chat = get_tmp(sc.d2chat,u)
    d3chat = get_tmp(sc.d3chat,u)
    d4chat = get_tmp(sc.d4chat,u)

    c = get_tmp(sc.c,u)
    dc = get_tmp(sc.dc,u)
    d2c = get_tmp(sc.d2c,u)
    d3c = get_tmp(sc.d3c,u)
    d4c = get_tmp(sc.d4c,u)

    T1 = get_tmp(sc.T1,u)
    T2 = get_tmp(sc.T2,u)
    T3 = get_tmp(sc.T3,u)
    T4 = get_tmp(sc.T4,u)

    c_l = get_tmp(sc.c_l,u)
    c_r = get_tmp(sc.c_r,u)
    μ_l = get_tmp(sc.μ_l,u)
    μ_r = get_tmp(sc.μ_r,u)

    RHS = get_tmp(sc.RHS,u)

    #Compute derivatives in Chebyshev space -> O(N)
    mul!(dchat,sn.D1,chat);
    mul!(d2chat,sn.D2,chat);
    mul!(d3chat,sn.D3,chat);
    mul!(d4chat,sn.D4,chat);

    #convert back to physical space -> O(NlogN)
    mul!(c,sn.US0_i,chat)
    # mul!(dc,sn.US1_i,dchat)
    # mul!(d2c,sn.US2_i,d2chat)
    # mul!(d3c,sn.US3_i,d3chat)
    # mul!(d4c,sn.US4_i,d4chat)
    dc = 2*sn.US1_i*dchat
    d2c = 2^2*sn.US2_i*d2chat
    d3c = 2^3*sn.US3_i*d3chat
    d4c = 2^4*sn.US4_i*d4chat

    #compute nonlinearities -> O(N)
    vmapntt!((c,d4c)->Term1(c,d4c,sp.κ),T1,c,d4c)
    vmapntt!((c,dc,d3c)->Term2(c,dc,d3c,sp.κ),T2,c,dc,d3c)
    vmapntt!((c,d2c)->Term3(c,d2c,sp.Ω),T3,c,d2c)
    vmapntt!((c,dc)->Term4(c,dc,sp.κ),T4,c,dc)

    RHS= T1 + T2 + T3 + T4
    RHS = sn.US4_f*RHS # Transform back to ultraspherical space
    
    # mul!(dchat_t,sn.US4_f,dc) # Transform back to ultraspherical space
    # dchat_t = sn.US4_f*dc
    @. du = RHS

    # Set up boundary conditions

    # Neumann in C at x=+-1
    du[end-3] = dot(sn.DChebEval_r,chat)
    du[end-2] = dot(sn.DChebEval_l,chat)

    #Rxn flux at x=+-1
    c_l = dot(sn.ChebEval_l,chat)
    c_r = dot(sn.ChebEval_r,chat)
    μ_l = μh(c_l,sp.Ω) - sp.κ*dot(sn.D2ChebEval_l,chat)
    μ_r = μh(c_r,sp.Ω) - sp.κ*dot(sn.D2ChebEval_r,chat)
    
    #Left Flux BC
    # @show (-(M(c_l)*-sp.κ*dot(sn.D3ChebEval_l,chat)) -BVrxn(c_l,μ_l,V,sp.k0,sp.α,sp.Vθ))
    du[end-1] =(-(M(c_l)*-sp.κ*dot(sn.D3ChebEval_l,chat)) -BVrxn(c_l,μ_l,V,sp.k0,sp.α,sp.Vθ))

    #Right Flux BC
    du[end] = ((M(c_r)*-sp.κ*dot(sn.D3ChebEval_r,chat)) -BVrxn(c_r,μ_r,V,sp.k0,sp.α,sp.Vθ))
    # @show c_l,c_r
    return nothing
end

function Jvp_CHR_US_vc(du,u,p,t,sp,sn,sc,Vfunc)
    chat =view(u,1:sp.Nx)
    dchat_t = view(du,1:sp.Nx)
    # dchat_t=du
    V = Vfunc(t)

    dchat = get_tmp(sc.dchat,u)
    d2chat = get_tmp(sc.d2chat,u)
    d3chat = get_tmp(sc.d3chat,u)
    d4chat = get_tmp(sc.d4chat,u)

    c = get_tmp(sc.c,u)
    dc = get_tmp(sc.dc,u)
    d2c = get_tmp(sc.d2c,u)
    d3c = get_tmp(sc.d3c,u)
    d4c = get_tmp(sc.d4c,u)

    T1 = get_tmp(sc.T1,u)
    T2 = get_tmp(sc.T2,u)
    T3 = get_tmp(sc.T3,u)
    T4 = get_tmp(sc.T4,u)

    c_l = get_tmp(sc.c_l,u)
    c_r = get_tmp(sc.c_r,u)
    μ_l = get_tmp(sc.μ_l,u)
    μ_r = get_tmp(sc.μ_r,u)

    RHS = get_tmp(sc.RHS,u)

    #Compute derivatives in Chebyshev space -> O(N)
    mul!(dchat,sn.D1,chat);
    mul!(d2chat,sn.D2,chat);
    mul!(d3chat,sn.D3,chat);
    mul!(d4chat,sn.D4,chat);

    #convert back to physical space -> O(NlogN)
    mul!(c,sn.US0_i,chat)
    # mul!(dc,sn.US1_i,dchat)
    # mul!(d2c,sn.US2_i,d2chat)
    # mul!(d3c,sn.US3_i,d3chat)
    # mul!(d4c,sn.US4_i,d4chat)
    # dc = 2*sn.US1_i*dchat
    # d2c = 2^2*sn.US2_i*d2chat
    # d3c = 2^3*sn.US3_i*d3chat
    # d4c = 2^4*sn.US4_i*d4chat

    #compute nonlinearities -> O(N)
    vmapntt!((c,d4c)->Term1(c,d4c,sp.κ),T1,c,d4c)
    vmapntt!((c,dc,d3c)->Term2(c,dc,d3c,sp.κ),T2,c,dc,d3c)
    vmapntt!((c,d2c)->Term3(c,d2c,sp.Ω),T3,c,d2c)
    vmapntt!((c,dc)->Term4(c,dc,sp.κ),T4,c,dc)

    dc= T1 + T2 + T3 + T4
    # mul!(dchat_t,sn.US4_f,dc) # Transform back to ultraspherical space
    # dchat_t = sn.US4_f*dc
    @. du = RHS

    # Set up boundary conditions

    # Neumann in C at x=+-1
    du[end-3] = dot(sn.DChebEval_r,chat)
    du[end-2] = dot(sn.DChebEval_l,chat)

    #Rxn flux at x=+-1
    c_l = dot(sn.ChebEval_l,chat)
    c_r = dot(sn.ChebEval_r,chat)
    μ_l = μh(c_l,sp.Ω) - sp.κ*dot(sn.D2ChebEval_l,chat)
    μ_r = μh(c_r,sp.Ω) - sp.κ*dot(sn.D2ChebEval_r,chat)
    
    #Left Flux BC
    # @show (-(M(c_l)*-sp.κ*dot(sn.D3ChebEval_l,chat)) -BVrxn(c_l,μ_l,V,sp.k0,sp.α,sp.Vθ))
    du[end-1] =(-(M(c_l)*-sp.κ*dot(sn.D3ChebEval_l,chat)) -BVrxn(c_l,μ_l,V,sp.k0,sp.α,sp.Vθ))

    #Right Flux BC
    du[end] = ((M(c_r)*-sp.κ*dot(sn.D3ChebEval_r,chat)) -BVrxn(c_r,μ_r,V,sp.k0,sp.α,sp.Vθ))
    # @show c_l,c_r
    return nothing
end

N = 1000;
sn = build_sim_numerics_S(N);
sc = init_cache_S(N);

Nx=N;
Lx = 1e-6 #m
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
C_rate = 4e0; #1/hr

mp = MaterialDomainParam(Lx,c_max,Ω_enthalpy,T,Vθ,κ,k0,α,C_rate,D,Nx,dx);
 
sp,sarg, sn_FV, sc_rhs_FV = set_up_sim(mp);

Vfunc(t)=115.0

CHR_S_RHS(du,u,p,t) = CHR_US_vc(du,u,p,t,sp,sn,sc,Vfunc)

u0 = zeros(N);
u0[1] = 0.05;
du0 = zeros(N);
p=[]
t0=0.0;

CHR_S_RHS(du0,u0,p,t0)
# using BenchmarkTools
@benchmark CHR_S_RHS($du0,$u0,$p,$t0)
