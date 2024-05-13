@inline function μ_regsoln(c,Ωmix)
    return log(c./(1.0 - c))+Ωmix*(1.0 - 2.0*c)
end


"""
    j0(c,μ,k0)

Calculate the exchange current density for a given concentration `c`, chemical potential `μ`, and rate constant `k0`.
"""
@inline function j0(c)
    # return k0*c*sqrt((1.0-c))
    return sqrt(max(1e-6,c))*(1.0-c)
end
"""
    BVrxn(c,μ,k0,α,Vθ)

Calculate the Butler-Volmer reaction rate for a given concentration `c`, chemical potential `μ`, rate constant `k0`, symmetry factor `α`, and standard potential `Vθ`.
"""

@inline function BVrxn(c,μ,V,α,Vθ)
    return j0(c)*(exp(-α*(μ + V - Vθ)) - exp((1-α)*(μ + V - Vθ)))
end

function ACR_FV_vc(du,u,p,t,sp,sn,sc,V)
    c =view(u,1:sp.Nx)
    dc = view(du,1:sp.Nx)

    CF_r = get_tmp(sc.CF_r_c, u)
    CF_l = get_tmp(sc.CF_l_c, u)
    μ = get_tmp(sc.μ_c, u)

    #Compute Concentration Fluxes
    mul!(CF_r,sn.∇x_f,c)
    mul!(CF_l,sn.∇x_b,c)

    vmapntt!(c->μ_regsoln(c,sp.Ω),μ,c)
    @. μ -= sp.κ*((CF_r-CF_l)/sp.dx)

    @. dc = BVrxn(c,μ,V(t),sp.α,sp.Vθ)

    nothing
end
function ACR_FV_cc(du,u,p,t,sp,sn,sc,V)
    c =view(u,1:sp.Nx)
    dc = view(du,1:sp.Nx)

    CF_r = get_tmp(sc.CF_r_c, u)
    CF_l = get_tmp(sc.CF_l_c, u)
    μ = get_tmp(sc.μ_c, u)

    #Compute Concentration Fluxes
    mul!(CF_r,sn.∇x_f,c)
    mul!(CF_l,sn.∇x_b,c)

    vmapntt!(c->μ_regsoln(c,sp.Ω),μ,c)
    @. μ -= sp.κ*((CF_r-CF_l)/sp.dx)

    @. dc = BVrxn(c,μ,V(t),sp.α,sp.Vθ)
    du[end] = sp.I - sum(dc)*sp.dx
    nothing
end


function CHR_FV_cc(du,u,p,t,sp,sn,sc)
    c =view(u,1:sp.Nx)
    dc = view(du,1:sp.Nx)


    CF_r = get_tmp(sc.CF_r_c, u)
    CF_l = get_tmp(sc.CF_l_c, u)

    μ = get_tmp(sc.μ_c, u)
    μF_r = get_tmp(sc.μF_r_c, u)
    μF_l = get_tmp(sc.μF_l_c, u)

    Cbar_r = get_tmp(sc.Cbar_r_c, u)
    Cbar_l = get_tmp(sc.Cbar_l_c, u)

    rxn = get_tmp(sc.rxn_c,u)
    
    #Compute Concentration Fluxes
    mul!(CF_r,sn.∇x_f,c)
    mul!(CF_l,sn.∇x_b,c)

    #Compute Chemical Potential
    vmapntt!(c->μ_regsoln(c,sp.Ω),μ,c)
    @. μ -= sp.κ*((CF_r-CF_l)/sp.dx)
    #Compute Chemical Potential Fluxes
    mul!(μF_r,sn.∇x_f, μ)
    mul!(μF_l,sn.∇x_b, μ)

   

    #Compute Interface Concentration Averages
    mul!(Cbar_r,sn.A_r,c)
    mul!(Cbar_l,sn.A_l,c)

    μF_r[end]=BVrxn(c[end],μ[end],u[end],sp.k0,sp.α,sp.Vθ)
    μF_l[1]=-BVrxn(c[1],μ[1],u[end],sp.k0,sp.α,sp.Vθ)

    @. dc = (getD(Cbar_r)*μF_r-getD(Cbar_l)*μF_l)/sp.dx

    du[end] = (μF_r[end] - μF_l[1]) - sp.I
    # du[end] = (-μF_l[1]) - sp.I
    # du[end] = (μF_r[end]) - sp.I
    return nothing
end

function CHR_FV_vc(du,u,p,t,sp,sn,sc,V)
    c =view(u,1:sp.Nx)
    dc = view(du,1:sp.Nx)


    CF_r = get_tmp(sc.CF_r_c, u)
    CF_l = get_tmp(sc.CF_l_c, u)

    μ = get_tmp(sc.μ_c, u)
    μF_r = get_tmp(sc.μF_r_c, u)
    μF_l = get_tmp(sc.μF_l_c, u)

    Cbar_r = get_tmp(sc.Cbar_r_c, u)
    Cbar_l = get_tmp(sc.Cbar_l_c, u)

    rxn = get_tmp(sc.rxn_c,u)
    
    #Compute Concentration Fluxes
    mul!(CF_r,sn.∇x_f,c)
    mul!(CF_l,sn.∇x_b,c)

    #Compute Chemical Potential
    vmapntt!(c->μ_regsoln(c,sp.Ω),μ,c)
    @. μ -= sp.κ*((CF_r-CF_l)/sp.dx)
    # @show μ
    #Compute Chemical Potential Fluxes
    mul!(μF_r,sn.∇x_f, μ)
    mul!(μF_l,sn.∇x_b, μ)

   

    #Compute Interface Concentration Averages
    mul!(Cbar_r,sn.A_r,c)
    mul!(Cbar_l,sn.A_l,c)

    μF_r[end]=BVrxn(c[end],μ[end],V(t),sp.k0,sp.α,sp.Vθ)
    μF_l[1]=-BVrxn(c[1],μ[1],V(t),sp.k0,sp.α,sp.Vθ)

    @. dc = (getD(Cbar_r)*μF_r-getD(Cbar_l)*μF_l)/sp.dx

    return nothing
end