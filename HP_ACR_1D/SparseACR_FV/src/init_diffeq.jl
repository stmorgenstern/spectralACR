function build_functions_cc(sp,sarg, sn, scr)
    CHR_RHS(du,u,p,t) = ACR_FV_cc(du,u,p,t,sp, sn, scr)
    return  CHR_RHS
end
function build_functions_vc(sp,sarg, sn, scr,V)
    CHR_RHS(du,u,p,t) = ACR_FV_vc(du,u,p,t,sp, sn, scr,V)
    return CHR_RHS
end

function build_forward_solve_cc(u0,sarg,p,RHSfunc)
    Nvar = length(u0)
    M= Diagonal(ones(Nvar));M[end,end]=0.0;
    du0 = similar(u0)
    jac_sparsity_cache = Symbolics.jacobian_sparsity((du, u) -> RHSfunc(du,u,p,0.0),du0,u0);
    # @show jac_sparsity_cache
    colorvec_cache = matrix_colors(jac_sparsity_cache);
    f_cache = ODEFunction(RHSfunc;jac_prototype=jac_sparsity_cache,colorvec=colorvec_cache,mass_matrix = M);
    sparse_prob_cache = ODEProblem(f_cache,u0,sarg.tspan,p);
    return sparse_prob_cache
end

function build_forward_solve_vc(u0,sarg,p,RHSfunc)
    Nvar = length(u0)
    du0 = similar(u0)
    jac_sparsity_cache = Symbolics.jacobian_sparsity((du, u) -> RHSfunc(du,u,p,0.0),du0,u0);
    colorvec_cache = matrix_colors(jac_sparsity_cache);
    f_cache = ODEFunction(RHSfunc;jac_prototype=jac_sparsity_cache,colorvec=colorvec_cache);
    sparse_prob_cache = ODEProblem(f_cache,u0,sarg.tspan,p);
    return sparse_prob_cache
end