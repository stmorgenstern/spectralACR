using CSV
using DataFrames
using ApproxFun
using BenchmarkTools
using Serialization

function generate_full_scaling()
    save_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\scaling_data\\"
    tracker = gen_scaling()
    for i=1:length(tracker)
        save_scale_tracker(tracker[i], save_path*"$(tracker[i].label)_scaling_data.jls")
    end
    return tracker
end
# Define the mutable struct with additional fields for accuracy
mutable struct scale_tracker
    label::String
    N::Vector{Int}
    sparse_time::Vector{Float64}
    dense_time::Vector{Float64}
    sparse_accuracy::Vector{Float64}
    dense_accuracy::Vector{Float64}
end

function grab_label(probtoken::probtoken)
    if probtoken.lithflag
        if probtoken.lowvoltage
            return "lith_low"
        else
            return "lith_high"
        end
    else
        if probtoken.lowvoltage
            return "delith_low"
        else
            return "delith_high"
        end
    end
end

# Function to load coefficients and time solutions
function load_solution(scenario::String, base_path::String)
    coeffs_path = joinpath(base_path, "$(scenario)_coeffs.csv")
    tsol_path = joinpath(base_path, "$(scenario)_tsol.csv")
    
    coeffs_df = CSV.read(coeffs_path, DataFrame,header=false,types=Float64)
    tsol_df = CSV.read(tsol_path, DataFrame,header=false,types=Float64)
    
    return coeffs_df, tsol_df
end

# Function to create and return a reference solution function
function create_ref_sol_func(scenario::String, base_path::String)
    coeffs_df, _ = load_solution(scenario, base_path)
    domain = Chebyshev(0..1)

    # Return a function that calculates the reference solution for given xgrid and time index
    return (xgrid, ti) -> begin
        output = zeros(length(xgrid),length(ti))
        for i=1:length(ti)
            coeffs_ti = Array(coeffs_df[:,i])
            chebyshev_func = Fun(domain, coeffs_ti)
            output[:,i] = chebyshev_func.(xgrid)
        end
        return output
    end
end

function gen_scaling()
    base_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SpectralACR\\ref_sols"
    label_list = ["lith_low", "lith_high", "delith_low", "delith_high"]
    probtokens = [probtoken(true,true), probtoken(true,false), probtoken(false,true), probtoken(false,false)]

    # label_list = ["lith_low"]
    # probtokens = [probtoken(true,false)]
    N_prob = length(probtokens)
    tracker = [scale_tracker(label,[],[],[],[],[]) for label in label_list]
    
    N_range = [2^i for i=6:12]
    # @warn "Currently using dev N_range for testing"
    # N_range = [2^i for i=6:7]

    for i=1:N_prob
        prob_token_i = probtokens[i]
        problabel = grab_label(prob_token_i)
        ref_sol_func = create_ref_sol_func(problabel, base_path)
        coeffs_df, tsol_df = load_solution(problabel, base_path)
        for N in N_range
            println("Benchmarking: ", problabel, " with N=", N)
            sparse_time, dense_time = calc_scaling(prob_token_i, N)
            sparse_accuracy, dense_accuracy = calc_accuracy(prob_token_i, N, ref_sol_func,tsol_df)
            
            push!(tracker[i].N, N)
            push!(tracker[i].sparse_time, sparse_time)
            push!(tracker[i].dense_time, dense_time)
            push!(tracker[i].sparse_accuracy, sparse_accuracy)
            push!(tracker[i].dense_accuracy, dense_accuracy)
        end
    end
    return tracker
end

function gen_scaling2()
    base_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SpectralACR\\ref_sols"
    label_list = [ "lith_high"]
    probtokens = [probtoken(true,false)]

    # label_list = ["lith_low"]
    # probtokens = [probtoken(true,false)]
    N_prob = length(probtokens)
    tracker = [scale_tracker(label,[],[],[],[],[]) for label in label_list]
    
    N_range = [2^i for i=6:15]
    # @warn "Currently using dev N_range for testing"
    # N_range = [2^i for i=6:7]

    for i=1:N_prob
        prob_token_i = probtokens[i]
        problabel = grab_label(prob_token_i)
        ref_sol_func = create_ref_sol_func(problabel, base_path)
        coeffs_df, tsol_df = load_solution(problabel, base_path)
        for N in N_range
            println("Benchmarking: ", problabel, " with N=", N)
            #sparse_time, dense_time = calc_scaling(prob_token_i, N)
            sparse_accuracy, dense_accuracy = calc_accuracy(prob_token_i, N, ref_sol_func,tsol_df)
            
            push!(tracker[i].N, N)
            #push!(tracker[i].sparse_time, sparse_time)
            #push!(tracker[i].dense_time, dense_time)
            push!(tracker[i].sparse_accuracy, sparse_accuracy)
            push!(tracker[i].dense_accuracy, dense_accuracy)
        end
    end
    return tracker
end


# Assume calc_scaling and calc_accuracy functions are defined and use `create_rhs`, `build_forward_solve_vc`, etc.
function calc_scaling(probtoken,N)
    sarg,p,ACR_RHS = create_rhs(N,probtoken)
    sparse_prob = build_forward_solve_vc(sarg.c0,sarg,p,ACR_RHS)
    dense_prob = ODEProblem(ACR_RHS,sarg.c0,sarg.tspan,p)
    sparse_time = @belapsed solve($sparse_prob,QNDF(linsolve=KLUFactorization()),save_everystep=false)
    dense_time = @belapsed solve($dense_prob,QNDF(),save_everystep=false)
    return sparse_time,dense_time
end
function calc_accuracy(probtoken,N,ref_sol_func,tsol_df)
    sarg,p,ACR_RHS = create_rhs(N,probtoken)
    sparse_prob = build_forward_solve_vc(sarg.c0,sarg,p,ACR_RHS)
    dense_prob = ODEProblem(ACR_RHS,sarg.c0,sarg.tspan,p)
    t_interval = LinRange(sarg.tspan[1],sarg.tspan[2],100)
    sparse_sol = solve(sparse_prob,QNDF(linsolve=KLUFactorization()),saveat=t_interval)
    dense_sol = solve(dense_prob,QNDF(),saveat=t_interval)
    # @show sparse_sol.t,Array(tsol_df)
    @assert sparse_sol.t ≈ Array(tsol_df)
    sparse_accuracy = calc_accuracy_helper(sparse_sol,N,ref_sol_func)
    dense_accuracy = calc_accuracy_helper(dense_sol,N,ref_sol_func)
    return sparse_accuracy,dense_accuracy
end
function calc_accuracy_helper(sol,N,ref_sol_func)
    xgrid = LinRange(0.0, 1.0, N)
    ref_sol = ref_sol_func(xgrid,sol.t)
    sol_c =Array(sol)
    return norm(sol_c-ref_sol,2)/norm(ref_sol,2)
end
function save_scale_tracker(tracker::scale_tracker, filename::String)
    io = open(filename, "w")
    serialize(io, tracker)
    close(io)
end

function load_scale_tracker(filename::String)::scale_tracker
    io = open(filename, "r")
    tracker = deserialize(io)
    close(io)
    return tracker
end

function load_scaling_data()
    save_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\scaling_data\\"
    label_list = ["lith_low", "lith_high", "delith_low", "delith_high"]
    tracker = []
    for label in label_list
        push!(tracker, load_scale_tracker(save_path*"$(label)_scaling_data.jls"))
    end
    return tracker
end

function gather_ref_sol_plotdat()
    base_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SpectralACR\\ref_sols"
    label_list = ["lith_low", "lith_high", "delith_low", "delith_high"]
    probtokens = [probtoken(true,true), probtoken(true,false), probtoken(false,true), probtoken(false,false)]

    N = 8192
    sol_collection = []
    xgrid = LinRange(0.0, 1.0,N)
    for i=1:length(label_list)
        prob_token_i = probtokens[i]
        sarg,p,ACR_RHS = create_rhs(N,prob_token_i)
        t_interval = LinRange(sarg.tspan[1],sarg.tspan[2],100)
        problabel = label_list[i]
        ref_sol_func = create_ref_sol_func(problabel, base_path)
        ref_sol = ref_sol_func(xgrid,t_interval)
        push!(sol_collection, ref_sol)
    end
    return sol_collection
end