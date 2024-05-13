mutable struct MaterialDomainParam
    #Physical Parameters
    Lx::Float64 #Reference length in x
    c_ref::Float64 #Reference concentration
    Ω::Float64 #Enthalpy of mixing
    T::Float64 #Temperature
    Vθ::Float64 #Standard Potential for intercalation
    κ::Float64 #Gradient Penalty
    k0_int::Float64 #Intrinsic rate constant
    α::Float64 #Symmetry factor
    C_rate::Float64 #C-rate

    #Numerical Parameters
    Nx::Int64 #Number of grid points in x direction
    dx::Float64 #Grid spacing in x direction

end

mutable struct sim_param
    Ω::Float64 #Enthalpy of mixing
    κ::Float64 #Gradient Penalty
    Vθ::Float64 #Standard Potential for intercalation
    #Lc_i::Float64 #Ratio of lengthscales

    α::Float64 #Symmetry factor
    I::Float64 #Dimensionless current

    # Numerical Parameters
    Nx::Int64 #Number of grid points in x direction

    dx::Float64 #Grid spacing in x direction

end
# Implement the show method for MaterialDomainParam
function Base.show(io::IO, mdp::MaterialDomainParam)
    print(io, "MaterialDomainParam(")
    println(io, "Lx = $(mdp.Lx),")
    println(io, "c_ref = $(mdp.c_ref),")
    println(io, "Ω = $(mdp.Ω),")
    println(io, "T = $(mdp.T),")
    println(io, "Vθ = $(mdp.Vθ),")
    println(io, "κ = $(mdp.κ),")
    println(io, "α = $(mdp.α),")
    println(io, "C_rate = $(mdp.C_rate),")
    println(io, "Nx = $(mdp.Nx),")
    println(io, "dx = $(mdp.dx),")

end

# Implement the show method for sim_param
function Base.show(io::IO, sp::sim_param)
    print(io, "sim_param(")
    println(io, "Ω = $(sp.Ω),")
    println(io, "κ = $(sp.κ),")
    println(io, "Vθ = $(sp.Vθ),")
    println(io, "α = $(sp.α),")
    println(io, "I = $(sp.I),")
    println(io, "Nx = $(sp.Nx),")
    println(io, "dx = $(sp.dx),")
end

mutable struct sim_numer_FV
    #Differential operators
    ∇x_b
    ∇x_f

end

mutable struct sim_cache_rhs_FV
    CF_r_c
    CF_l_c

    μ_c

    rxn_c

end


mutable struct sim_arg
    c0
    tspan
end

