@inline function ar_mean(x,y)
    return 0.5*(x+y)
end

"""
    nondimWrapper(c0,tf,param::simParamDim)

Normalizes the initial concentration `c0` and final time `tf`, and converts the dimensional parameters in `param` to dimensionless parameters. Returns the dimensionless initial concentration, final time, and a `simParamNonDim` object with the dimensionless parameters.
"""
function nondimWrapper(c0,param)
    param_nondim = nondimParam(param)
    c0_tilde,tf_tilde = nondimArg(c0,param,param_nondim)
    return c0_tilde,tf_tilde,param_nondim
end
"""
    nondimParam(param::simParamDim)

Converts the dimensional parameters in `param` to dimensionless parameters. Returns a `simParamNonDim` object with the dimensionless parameters.
"""
function nondimParam(param::MaterialDomainParam)
    #e = 1.60217662e-19; #Elementary Charge [C]
    #kB = 1.38064852e-23; #Boltzmann Constant [J/(molecules*K)]
    NA = 	6.02214076e23; #Avogadro's Number [molecules/mol]
    F = 9.64853321233100184e4 #C⋅mol−1.
    e = 1.0; #Elementary Charge [e]
    kB = 8.617333262145e-5; #Boltzmann Constant [eV/K]

    τ = 1/3600; # hour to s conversion [hr/s]

    L_ref = param.Lx

    Ω_tilde = param.Ω/(kB*param.T)
    κ_tilde = param.κ/(param.c_ref*((L_ref)^2)*NA*kB*param.T)
    Vθ_tilde = param.Vθ/(kB*param.T/e)

    # @warn "double check k0 nondim to make sure units are right"
    #k0_tilde = L_ref*param.k0_int/(param.D0*param.c_ref*F)

    # @warn "Sam: Double check the C-Rate non-dimensionalization there might need to be a k0 in there"
    # C_rate_tilde = param.C_rate*param.c_ref*param.H*F*τ/param.k0_int 
    # I_tilde = param.C_rate*τ*A/param.D0
    # I_tilde = param.C_rate*τ*L_ref^2/param.D0
    I_tilde = param.C_rate*τ*param.c_ref*F*L_ref/param.k0_int

    dx_tilde = param.dx/L_ref


    # @warn "Overriding dx with dx=.01"
    # dx_tilde =0.01;
    # @warn "Overriding κ with κ=.001"
    # κ_tilde =0.002;
    # @warn "Overriding Ω with Ω=3"
    # Ω_tilde =3.0;
    # @warn "Overriding k0 with k0=1e-1, previously was k0=$(k0_tilde)"
    # k0_tilde =1e-1;

    param_nondim = sim_param(Ω_tilde,κ_tilde,Vθ_tilde,param.α,I_tilde,param.Nx,dx_tilde)
    #param_nondim = sim_param(Ω_tilde,κ_tilde,Vθ_tilde,k0_tilde,param.α,I_tilde,param.Nx,param.Ny,dx_tilde,dy_tilde)
    return param_nondim
end
"""
    nondimArg(c0,tf,param::simParamDim)

Normalizes the initial concentration `c0` and final time `tf` by the reference concentration and diffusion time scale, respectively. Returns the dimensionless initial concentration and final time.
"""
function nondimArg(c0,mp,param_nondim)
    F = 9.64853321233100184e4
    trxn = (mp.Lx*F*mp.c_ref)/(mp.k0_int)
    tprocess = 3600*(1/mp.C_rate)
    tprocess_tilde = tprocess/trxn
    #tdiff= (mp.Lx^2)/(mp.D0)

    #trxn_tilde = trxn/tdiff
    # tprocess_tilde = tprocess/trxn_tilde
    #@show tprocess_tilde
    
    #tf_tilde = max(trxn_tilde,tprocess_tilde)
    tf_tilde = tprocess_tilde

    c0_tilde = c0/mp.c_ref;
    return c0_tilde,tf_tilde
end
"""
    createMaskWrapper(raw_domain,Δx, Δy;Δx_sg=1.0,Δy_sg=1.0, Δt=0.01, Nt=10000, υ=1e-7, ζ_fac=4.0)

Calculates the mask and its gradient for the given `raw_domain` and grid spacing `Δx` and `Δy`. Returns the mask, its gradient, the boundary mask, and the area of the mask.
"""


function filter(x)
    x[isnan.(x)] .= 0.0
    x[x.==Inf] .= 0.0
    return x
end

function build_sim_numerics(sim_params)
    dx= sim_params.dx
    Nx = sim_params.Nx

    # Forward finite difference operator in x direction
    ∇x_f = Bidiagonal([-1.0 for i in 1:Nx],[1.0 for i in 1:Nx-1],:U)
    ∇x_f[end,end] = 0.0
    ∇x_f ./= dx
    # Backwards finite difference operator in x direction
    ∇x_b =Bidiagonal([1.0 for i in 1:Nx],[-1.0 for i in 1:Nx-1], :L);
    ∇x_b[1,1] = 0.0
    ∇x_b ./= dx

    sn = sim_numer_FV(∇x_b,∇x_f)
    
    return sn
end

function init_cache(sim_params;chunk_size=25)
    Nx =sim_params.Nx

    #Set up Flux/Potential Caches
    CF_r = zeros(Nx); CF_r_c = DiffCache(CF_r, chunk_size);
    CF_l = zeros(Nx); CF_l_c = DiffCache(CF_l, chunk_size);
    
    μ = zeros(Nx); μ_c = DiffCache(μ, chunk_size);
    
    rxn = zeros(Nx); rxn_c = DiffCache(rxn, chunk_size);

    sim_c_rhs = sim_cache_rhs_FV(CF_r_c,CF_l_c,μ_c,rxn_c)
    
    return sim_c_rhs
end

function set_up_sim(mp)
    # c0 =mp.c_ref*(.5*ones(mp.Nx) + 0.01*rand(mp.Nx)) 
    c0 = mp.c_ref*0.95*ones(mp.Nx)
    t0 = 0.0;
    F = 9.64853321233100184e4
    trxn = (mp.Lx*F*mp.c_ref)/(mp.k0_int)
    tprocess = 3600*(1/mp.C_rate)
    #tdiff= (mp.Lx^2)/(mp.D0)

    #trxn_tilde = trxn/tdiff
    #tprocess_tilde = tprocess/tdiff

    # tf = (1/mp.C_rate)*3600#*A_frac; #s
    # tf = 1.0;
    @warn "Model is currently set up in charging mode only"


    c0_tilde,tf_tilde,param_nondim = nondimWrapper(c0,mp);
    #Wetting BCs
    c0_tilde[1] = 0.01;
    c0_tilde[end] = 0.01;

    @show tf_tilde
    # tf_tilde = 5000.0

    sarg = sim_arg(c0_tilde,(t0,tf_tilde))
    if param_nondim.dx > sqrt(param_nondim.κ/param_nondim.Ω)
        @warn "Warning: dx > sqrt(κ/Ω)"
    end

    sn = build_sim_numerics(param_nondim)

    sim_c_rhs = init_cache(param_nondim;chunk_size=25)

    return param_nondim,sarg, sn, sim_c_rhs
end

function set_up_sim_bench(mp,probtoken)
    lithflag =probtoken.lithflag
    lowvoltage = probtoken.lowvoltage

    xgrid = LinRange(0.0, 1.0, mp.Nx);
    c0_nondim = lithflag ? build_ic_lithiation.(xgrid) : build_ic_delithiation.(xgrid)

    c0 = mp.c_ref*c0_nondim
    t0 = 0.0;
    F = 9.64853321233100184e4
    trxn = (mp.Lx*F*mp.c_ref)/(mp.k0_int)
    
    
    c0_tilde,tf_tilde_dummy,param_nondim = nondimWrapper(c0,mp);

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
    
    sarg = sim_arg(c0_tilde,(t0,tf_tilde))
    if param_nondim.dx > sqrt(param_nondim.κ/param_nondim.Ω)
        @warn "Warning: dx > sqrt(κ/Ω)"
    end

    sn = build_sim_numerics(param_nondim)

    sim_c_rhs = init_cache(param_nondim;chunk_size=25)

    return param_nondim,sarg, sn, sim_c_rhs
end