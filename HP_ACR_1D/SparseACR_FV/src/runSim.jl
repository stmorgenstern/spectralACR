
function runSim_cc(mp::MaterialDomainParam)
    sp,sarg, sn, sc_rhs = set_up_sim(mp)
    @show sp
    ACR_RHS = build_functions_cc(sp,sarg, sn, sc_rhs)
    t0,tf_tilde = sarg.tspan
    ts = range(t0,tf_tilde;length=100)
    # p = zeros(sp.Np)

    cavg = sum(sarg.c0[:])/sp.Nx
    μavg = μ_regsoln(cavg,sp.Ω)
    Vguess = sp.Vθ - μavg -asinh(sp.I/(4*j0(cavg,μavg,sp.k0)))
    @show Vguess

    u0 = [sarg.c0[:] ; Vguess];
    du0 = similar(u0);

    @warn "Model current is not parameterized through DifferentialEquations.jl parameter so p is blank"
    p =[]
    println("RHS eval")
    @time CHR_RHS(du0,u0,p,0.0);
    println("")
    

    prob = build_forward_solve_cc(u0,sarg,p,CHR_RHS)

    condition(u, t, integrator) = mean(u[1:end-1]) > 0.99
    # condition(u, t, integrator) = condition_psi(u, t, integrator,sn.ψb)
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)

    @time sol = solve(prob,FBDF(),initializealg = BrownFullBasicInit(),callback=cb,progress=true)
    
    return sp,sarg,sn,sol
end
#Create Problem
function build_ic_delithiation(x)
    return 0.49 * tanh((x - 0.15) / 0.05) + 0.5
end
function build_ic_lithiation(x)
    return 0.49 * tanh((x - 0.85) / 0.05) + 0.5
end
function runSim_vc(mp::MaterialDomainParam,V_vc; nondimflag=false)
    sp,sarg, sn, sc_rhs = set_up_sim(mp)
    @show sp
    e = 1.0; #Elementary Charge [e]
    kB = 8.617333262145e-5; #Boltzmann Constant [eV/K]
    V(t) = nondimflag ? V_vc(t)/(kB*mp.T/e) : V_vc(t)

    CHR_RHS = build_functions_vc(sp,sarg, sn, sc_rhs,V)
    t0,tf_tilde = sarg.tspan
    sarg.tspan = (t0,tf_tilde)
    ts = range(t0,tf_tilde;length=100)
    # p = zeros(sp.Np)

    #u0 = sarg.c0[:]
    u0 = build_ic_lithiation.(LinRange(0.0, 1.0, sp.Nx))
    sarg.c0 = u0
    du0 = similar(u0);

    @warn "overriding tf"
    sarg.tspan = (0.0, 2.50)
    @warn "Model current is not parameterized through DifferentialEquations.jl parameter so p is blank"
    p =[]
    println("RHS eval")
    @time CHR_RHS(du0,u0,p,0.0);
    println("")
    @show du0

    prob = build_forward_solve_vc(u0,sarg,p,CHR_RHS)

    # condition_psi(u, t, integrator,ψb) = 
    condition(u, t, integrator) = mean(u) < .01
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition, affect!)

    println("Forward Solve")
    @time sol = solve(prob,FBDF(),callback=cb,progress=true)#,dtmax=sqrt(sp.κ/sp.Ω))
        # @btime solve($prob,FBDF(),initializealg = BrownFullBasicInit(),callback=$cb)
    return sp,sarg,sn,sol
end