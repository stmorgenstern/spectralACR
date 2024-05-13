using DiffEqDevTools
#Create Problem
function build_ic_delithiation(x)
    return 0.49 * tanh((x - 0.15) / 0.05) + 0.5
end
function build_ic_lithiation(x)
    return 0.49 * tanh((x - 0.85) / 0.05) + 0.5
end

function build_LFP_test_mp(N)
    
    Lx = 1e-7 #m
    x = LinRange(0.0, Lx, N)
    dx = x[2] - x[1];


    # Set up Material Parameters
    c_max = 2.29e4 #mol/m^3
    Ω_enthalpy = 0.115; #eV\
    T = 298.15; #K
    Vθ = 3.42; #V

    κ = 3.1314026858813996e9; #eV/m
    C_rate = 1e-12;

    k0 =1.6e-4

    α = 0.5;

    mp = MaterialDomainParam(Lx,c_max,Ω_enthalpy,T,Vθ,κ,k0,α,C_rate,N,dx);

    return mp
end
struct probtoken
    lithflag::Bool
    lowvoltage::Bool
end
function create_problem(N,probtoken)
    lithflag =probtoken.lithflag
    lowvoltage = probtoken.lowvoltage

    if lowvoltage && lithflag
        # Slow filling at lowvoltage
        eta_tilde = -0.1
        tf_tilde=758.44
    elseif lowvoltage && !lithflag
        # Slow emptying at lowvoltage
        eta_tilde = 0.1
        tf_tilde=780.23
    elseif !lowvoltage && lithflag
        # Fast filling at high voltage
        eta_tilde = -2.0
        tf_tilde=2.50
    else
        # Fast emptying at high voltage
        eta_tilde = 2.0
        tf_tilde=4.54
    end
    # @show N
    mp = build_LFP_test_mp(N)
    sp,sarg, sn, sc_rhs = set_up_sim_bench(mp,probtoken)


    Vfunc(t)=(eta_tilde+sp.Vθ)
    
    #@show sp.Nx
    ACR_RHS = build_functions_vc(sp,sarg, sn, sc_rhs,Vfunc)

    p=[]
    prob = build_forward_solve_vc(sarg.c0,sarg,p,ACR_RHS)
    return prob
end
function create_rhs(N,probtoken)
    lithflag =probtoken.lithflag
    lowvoltage = probtoken.lowvoltage

    if lowvoltage && lithflag
        # Slow filling at lowvoltage
        eta_tilde = -0.1
        tf_tilde=758.44
    elseif lowvoltage && !lithflag
        # Slow emptying at lowvoltage
        eta_tilde = 0.1
        tf_tilde=780.23
    elseif !lowvoltage && lithflag
        # Fast filling at high voltage
        eta_tilde = -2.0
        tf_tilde=2.50
    else
        # Fast emptying at high voltage
        eta_tilde = 2.0
        tf_tilde=4.54
    end
    # @show N
    mp = build_LFP_test_mp(N)
    sp,sarg, sn, sc_rhs = set_up_sim_bench(mp,probtoken)


    Vfunc(t)=(eta_tilde+sp.Vθ)
    
    #@show sp.Nx
    ACR_RHS = build_functions_vc(sp,sarg, sn, sc_rhs,Vfunc)

    p=[]

    return sarg,p,ACR_RHS
end
function create_work_precision(N,probtoken)
    prob = create_problem(N,probtoken)
    test_sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14)
    # High Precision Tols
    abstols = 1.0 ./ 10.0 .^ (7:12)
    reltols = 1.0 ./ 10.0 .^ (4:9)
    
    #Low Precision Tols
    # abstols = 1.0 ./ 10.0 .^ (5:8)
    # reltols = 1.0 ./ 10.0 .^ (1:4);

    setups=[
    #Stiff solvers -> Default Linsolve
    Dict(:alg=>CVODE_BDF()),
    Dict(:alg=>TRBDF2()),
    Dict(:alg=>QNDF()),
    Dict(:alg=>FBDF()),
    Dict(:alg=>KenCarp4()),
    Dict(:alg=>Rosenbrock23()),

    #Stiff solvers -> Sparse KLU Linsolve
    Dict(:alg=>TRBDF2(linsolve=KLUFactorization())),
    Dict(:alg=>QNDF(linsolve=KLUFactorization())),
    Dict(:alg=>FBDF(linsolve=KLUFactorization())),
    Dict(:alg=>KenCarp4(linsolve=KLUFactorization())),
    Dict(:alg=>Rosenbrock23(linsolve=KLUFactorization()))
    ]

    labels = [
            "CVODE_BDF", "TRBDF2", "QNDF", "FBDF", "KenCarp4", "Rosenbrock23",
            "TRBDF2_KLU", "QNDF_KLU", "FBDF_KLU", "KenCarp4_KLU", "Rosenbrock23_KLU"]
    wp = WorkPrecisionSet(prob,abstols,reltols,setups;names = labels,
                      save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),numruns=10)
    return wp
end
function getproblabel(probtoken::probtoken)
    if probtoken.lithflag
        if probtoken.lowvoltage
            return "Lithiation Low Voltage"
        else
            return "Lithiation High Voltage"
        end
    else
        if probtoken.lowvoltage
            return "Delithiation Low Voltage"
        else
            return "Delithiation High Voltage"
        end
    end
end
# function get_end_times(N)
#     lith_low = probtoken(true,true)
#     lith_high = probtoken(true,false)
#     delith_low = probtoken(false,true)
#     delith_high = probtoken(false,false)

#     probtokens = [lith_low,lith_high,delith_low,delith_high]    
#     sols = []
#     for i = 1:length(probtokens)
#         prob_token_i = probtokens[i]
#         prob = create_problem(N,prob_token_i)
#         prob_bool(u) = prob_token_i.lithflag ? mean(u) > .9875 : mean(u) < .0125

#         condition(u, t, integrator) = prob_bool(u)
#         affect!(integrator) = terminate!(integrator)
#         cb = DiscreteCallback(condition, affect!)
#         # test_sol = solve(prob,CVODE_BDF(),callback=cb)
#         test_sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14,callback=cb)
#         push!(sols,test_sol)
#         problabel = getproblabel(prob_token_i)
#         println("End time for "*problabel*": ",test_sol.t[end])
#     end
#     return sols
# end
function run_bench(N)
    lith_low = probtoken(true,true)
    lith_high = probtoken(true,false)
    delith_low = probtoken(false,true)
    delith_high = probtoken(false,false)
    probtokens = [lith_low,lith_high,delith_low,delith_high]
    wp_set = []
    for i=1:length(probtokens)
        prob_token_i = probtokens[i]
        problabel = getproblabel(prob_token_i)
        println("Benchmarking: ",problabel)
        wp = create_work_precision(N,prob_token_i)
        push!(wp_set,wp)
    end
    return wp_set
end
function generate_WP_diagrams(N)
    wp_set = run_bench(N)
    labels=["lith_low_wp","lith_high_wp","delith_low_wp","delith_high_wp"]
    for i=1:length(wp_set)
        wp_i = plot(wp_set[i])
        save("C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\figures\\work_precision\\"*labels[i]*".pdf",wp_i)
    end
end

function test_run(N)
    #Screen the solvers to see what to benchmark with w/ Work Precision Diagrams
    setups=[
    #Stabilized Explicit Solvers

    #Stiff solvers -> Default Linsolve
    Dict(:alg=>CVODE_BDF()),
    Dict(:alg=>TRBDF2()),
    Dict(:alg=>QNDF()),
    Dict(:alg=>FBDF()),
    Dict(:alg=>KenCarp4()),
    Dict(:alg=>Rosenbrock23()),
    Dict(:alg=>Rodas5()),
    #Stiff solvers -> Krylov Linsolve
    Dict(:alg=>CVODE_BDF(linear_solver=:GMRES)),
    Dict(:alg=>TRBDF2(linsolve=KrylovJL_GMRES())),
    Dict(:alg=>QNDF(linsolve=KrylovJL_GMRES())),
    Dict(:alg=>FBDF(linsolve=KrylovJL_GMRES())),
    Dict(:alg=>KenCarp4(linsolve=KrylovJL_GMRES())),
    Dict(:alg=>Rosenbrock23(linsolve=KrylovJL_GMRES())),
    Dict(:alg=>Rodas5(linsolve=KrylovJL_GMRES())),
    # Stiff solvers -> Sparse KLU Linsolve
    Dict(:alg=>TRBDF2(linsolve=KLUFactorization())),
    Dict(:alg=>QNDF(linsolve=KLUFactorization())),
    Dict(:alg=>FBDF(linsolve=KLUFactorization())),
    Dict(:alg=>KenCarp4(linsolve=KLUFactorization())),
    Dict(:alg=>Rosenbrock23(linsolve=KLUFactorization())),
    Dict(:alg=>Rodas5(linsolve=KLUFactorization()))
    ]

    lith_low = probtoken(true,true)
    lith_high = probtoken(true,false)
    delith_low = probtoken(false,true)
    delith_high = probtoken(false,false)
    probtokens = [lith_low,lith_high,delith_low,delith_high]

    labels = [
            "CVODE_BDF", "TRBDF2", "QNDF", "FBDF", "KenCarp4", "Rosenbrock23", "Rodas5",
            "CVODE_BDF_GMRES", "TRBDF2_GMRES", "QNDF_GMRES", "FBDF_GMRES", "KenCarp4_GMRES", "Rosenbrock23_GMRES", "Rodas5_GMRES",
            "TRBDF2_KLU", "QNDF_KLU", "FBDF_KLU", "KenCarp4_KLU", "Rosenbrock23_KLU", "Rodas5_KLU"]

    for i=1:length(probtokens)
        prob_token_i = probtokens[i]
        problabel = getproblabel(prob_token_i)
        prob = create_problem(N,prob_token_i)
        println("Benchmarking: ",problabel)
        for j=1:length(setups)
            try 
                sol_j = solve(prob,setups[j][:alg])
                #println("Passed on ",problabel," with ",setups[j])
                if sol_j.retcode != :Success
                    println("Failed on ",problabel," with ",labels[j])
                end
            catch 
                println("Failed on ",problabel," with ",labels[j])
            end
        end
    end

end