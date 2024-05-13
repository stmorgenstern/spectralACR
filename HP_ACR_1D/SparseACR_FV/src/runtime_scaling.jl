# mutable struct scale_tracker
#     label::String
#     N::Vector{Int}
#     sparse_time::Vector{Float64}
#     dense_time::Vector{Float64}
#     sparse_accuracy::Vector{Float64}
#     dense_accruacy::Vector{Float64}
# end


# function gen_scaling()
#     lith_low = probtoken(true,true)
#     lith_high = probtoken(true,false)
#     delith_low = probtoken(false,true)
#     delith_high = probtoken(false,false)
#     probtokens = [lith_low,lith_high,delith_low,delith_high]
#     N_prob = length(probtokens)
#     tracker=[scale_tracker("lith_low",[],[],[]),scale_tracker("lith_high",[],[],[]),scale_tracker("delith_low",[],[],[]),scale_tracker("delith_high",[],[],[])]
    
#     N_range = [2^i for i=6:12]
#     label_list = ["lith_low","lith_high","delith_low","delith_high"]

#     for i=1:N_prob
#         prob_token_i = probtokens[i]
#         problabel = getproblabel(prob_token_i)
#         for N in N_range
#             println("Benchmarking: ",problabel," with N=",N)
#             sparse_time,dense_time = calc_scaling(prob_token_i,N)
#             push!(tracker[i].N,N)
#             push!(tracker[i].sparse_time,sparse_time)
#             push!(tracker[i].dense_time,dense_time)
#         end
#     end
#     return tracker
# end
# function calc_scaling(probtoken,N)
#     sarg,p,ACR_RHS = create_rhs(N,probtoken)
#     sparse_prob = build_forward_solve_vc(sarg.c0,sarg,p,ACR_RHS)
#     dense_prob = ODEProblem(ACR_RHS,sarg.c0,sarg.tspan,p)
#     sparse_time = @belapsed solve($sparse_prob,QNDF(linsolve=KLUFactorization()),save_everystep=false)
#     dense_time = @belapsed solve($dense_prob,QNDF(),save_everystep=false)
#     return sparse_time,dense_time
# end
# function calc_accuracy(probtoken,N,ref_sol_func)
#     sarg,p,ACR_RHS = create_rhs(N,probtoken)
#     sparse_prob = build_forward_solve_vc(sarg.c0,sarg,p,ACR_RHS)
#     dense_prob = ODEProblem(ACR_RHS,sarg.c0,sarg.tspan,p)
#     t_interval = LinRange(sarg.tspan[1],sarg.tspan[2],100)
#     sparse_sol = solve(sparse_prob,QNDF(linsolve=KLUFactorization()),saveat=t_interval)
#     dense_sol = solve(dense_prob,QNDF(),saveat=t_interval)
#     sparse_accuracy = calc_accuracy_helper(sparse_sol,N,ref_sol_func)
#     dense_accuracy = calc_accuracy_helper(dense_sol,N,ref_sol_func)
#     return sparse_accuracy,dense_accuracy
# end
# function calc_accuracy_helper(sol,N,ref_sol_func)
#     xgrid = LinRange(0.0, 1.0, N)
#     ref_sol = ref_sol_func(xgrid)
#     sol_c = sol(1.0)[1:N]
#     return norm(sol_c-ref_sol)/norm(ref_sol)
    
# end