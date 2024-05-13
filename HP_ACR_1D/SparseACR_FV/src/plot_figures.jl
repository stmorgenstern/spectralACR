using CSV, DataFrames,Colors

mutable struct method_plot_data
    N::Vector{Int}
    accuracy::Vector{Float64}
    time::Vector{Float64}
end
mutable struct scenario_plot_data
    label::String
    spectral::method_plot_data
    sparse::method_plot_data
    dense::method_plot_data
end

function load_spectral_data(scenario_label)


    # Define the path to the Spectral Collocation data CSV file
    spectral_data_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SpectralACR\\scaling_dat\\$(scenario_label)_scaledat.csv"
    
    # Load Spectral Collocation data from CSV
    spectral_df = CSV.read(spectral_data_path, DataFrame,header=false,types=Float64)
    
    # Prepare data from CSV
    spectral_N = collect(spectral_df[1, :])
    spectral_accuracy = collect(spectral_df[2, :])
    spectral_time = collect(spectral_df[3, :])
    
    return method_plot_data(spectral_N, spectral_accuracy, spectral_time)
end
function convert_tracker_plotdata(tracker::scale_tracker)
    N = tracker.N
    sparse_accuracy = tracker.sparse_accuracy
    dense_accuracy = tracker.dense_accuracy
    sparse_time = tracker.sparse_time
    dense_time = tracker.dense_time
    return tracker.label, method_plot_data(N, sparse_accuracy, sparse_time), method_plot_data(N, dense_accuracy, dense_time)
end
function gen_plot_data()
    tracker_list = load_scaling_data()
    plot_data = []
    for tracker in tracker_list
        scenario_label, sparse_plot_data, dense_plot_data = convert_tracker_plotdata(tracker)
        spectral_plot_data = load_spectral_data(scenario_label)
        push!(plot_data, scenario_plot_data(scenario_label, spectral_plot_data, sparse_plot_data, dense_plot_data))
    end
    return plot_data
end
using Makie

function plot_N_vs_Time(scenario_data::scenario_plot_data)
    CairoMakie.activate!() # Activate CairoMakie backend
    
    fig = Figure()
    ax = Axis(fig[1, 1], xscale=Makie.log10, yscale=Makie.log10,
              xlabel="N", ylabel="Simulation Time", title="$(scenario_data.label) - N vs Simulation Time")
    
    # Plotting for Spectral Collocation
    scatterlines!(ax, scenario_data.spectral.N, scenario_data.spectral.time, color=:blue, label="SC - PDE15s")
    
    # Plotting for Finite Volume Dense Solver
    scatterlines!(ax, scenario_data.dense.N, scenario_data.dense.time, color=:red, label="FV - QNDF() Dense")
    
    # Plotting for Finite Volume Sparse Solver
    scatterlines!(ax, scenario_data.sparse.N, scenario_data.sparse.time, color=:green, label="FV - QNDF() Sparse")
    
    Legend(fig[2,1], ax, "Methods",orientation=:horizontal)
    hidedecorations!(ax,label=false,ticks=false,ticklabels=false)
    file_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\figures\\scaling\\$(scenario_data.label)_N_vs_Time_plot.pdf"
    save(file_path, fig) # Save the figure to a file
    
    return fig
end
function plot_N_vs_Error(scenario_data::scenario_plot_data)
    CairoMakie.activate!() # Activate CairoMakie backend
    
    fig = Figure()
    ax = Axis(fig[1, 1], xscale=Makie.log10, yscale=Makie.log10,
              xlabel="N", ylabel="Relative Error", title="$(scenario_data.label) - N vs Relative Error")
    
    # Plotting for Spectral Collocation
    scatterlines!(ax, scenario_data.spectral.N, scenario_data.spectral.accuracy, color=:blue, label="SC - PDE15s")
    
    # Plotting for Finite Volume Sparse Solver
    scatterlines!(ax, scenario_data.sparse.N, scenario_data.sparse.accuracy, color=:green, label="FV - QNDF()")
    
    Legend(fig[2,1], ax, "Methods",orientation=:horizontal)
    hidedecorations!(ax,label=false,ticks=false,ticklabels=false)
    file_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\figures\\scaling\\$(scenario_data.label)_N_vs_Err_plot.pdf"
    save(file_path, fig) # Save the figure to a file
    
    return fig
end
function plot_Err_vs_Time(scenario_data::scenario_plot_data)
    CairoMakie.activate!() # Activate CairoMakie backend
    
    fig = Figure()
    ax = Axis(fig[1, 1], xscale=Makie.log10, yscale=Makie.log10,
              xlabel="Relative Error", ylabel="Time", title="$(scenario_data.label) - Relative Error vs Time")
    
    # Plotting for Spectral Collocation
    scatterlines!(ax, scenario_data.spectral.accuracy, scenario_data.spectral.time, color=:blue, label="SC - PDE15s")
    
    # Plotting for Finite Volume Dense Solver
    scatterlines!(ax, scenario_data.dense.accuracy, scenario_data.dense.time, color=:red, label="FV - QNDF() Dense")
    
    # Plotting for Finite Volume Sparse Solver
    scatterlines!(ax, scenario_data.sparse.accuracy, scenario_data.sparse.time, color=:green, label="FV - QNDF() Sparse")
    
    Legend(fig[2,1], ax, "Methods",orientation=:horizontal)
    hidedecorations!(ax,label=false,ticks=false,ticklabels=false)
    file_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\figures\\scaling\\$(scenario_data.label)_Err_vs_Time_plot.pdf"
    save(file_path, fig) # Save the figure to a file
    
    return fig
end


function plot_N_vs_Time_dual(scenario_data1::scenario_plot_data, scenario_data2::scenario_plot_data, label1::String, label2::String,filename)
    CairoMakie.activate!() # Activate CairoMakie backend
    
    # Define the size in inches
    width_in_inches = 6.25
    height_in_inches = 4.29

    # Convert inches to pixels (assuming 96 dpi)
    width_in_pixels = round(Int, width_in_inches * 96)
    height_in_pixels = round(Int, height_in_inches * 96)

    # Create a figure with the specified size
    fig = Figure(size = (width_in_pixels, height_in_pixels))
    ax1 = Axis(fig[1, 1], xscale=Makie.log10, yscale=Makie.log10,
              xlabel="N", ylabel="Simulation Time", title=label1)
    ax2 = Axis(fig[1, 2], xscale=Makie.log10, yscale=Makie.log10,
              xlabel="N", ylabel="Simulation Time", title=label2)
    linkyaxes!(ax1, ax2)
    # Plotting for Spectral Collocation
    scatterlines!(ax1, scenario_data1.spectral.N, scenario_data1.spectral.time, color=:blue, label="SC - PDE15s")
    
    # Plotting for Finite Volume Dense Solver
    scatterlines!(ax1, scenario_data1.dense.N, scenario_data1.dense.time, color=:red, label="FV - QNDF() Dense")
    
    # Plotting for Finite Volume Sparse Solver
    scatterlines!(ax1, scenario_data1.sparse.N, scenario_data1.sparse.time, color=:green, label="FV - QNDF() Sparse")
    
    # Plotting for Spectral Collocation
    scatterlines!(ax2, scenario_data2.spectral.N, scenario_data2.spectral.time, color=:blue, label="SC - PDE15s")
    
    # Plotting for Finite Volume Dense Solver
    scatterlines!(ax2, scenario_data2.dense.N, scenario_data2.dense.time, color=:red, label="FV - QNDF() Dense")
    
    # Plotting for Finite Volume Sparse Solver
    scatterlines!(ax2, scenario_data2.sparse.N, scenario_data2.sparse.time, color=:green, label="FV - QNDF() Sparse")

    yspace = maximum(tight_yticklabel_spacing!, [ax1, ax2])
    ax1.yticklabelspace = yspace
    ax2.yticklabelspace = yspace
    Legend(fig[2,1:2], ax1, "Methods",orientation=:horizontal)
    hidedecorations!(ax1,label=false,ticks=false,ticklabels=false)
    hidedecorations!(ax2,label=false,ticks=false,ticklabels=false)
    file_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\figures\\scaling\\prod\\$(filename)_N_vs_Time_plot.pdf"
    save(file_path, fig) # Save the figure to a file
    
    return fig
end
function plot_N_vs_Error_dual(scenario_data1::scenario_plot_data, scenario_data2::scenario_plot_data, label1::String, label2::String, filename)
    CairoMakie.activate!() # Activate CairoMakie backend

    # Define the size in inches
    width_in_inches = 6.25
    height_in_inches = 4.29

    # Convert inches to pixels (assuming 96 dpi)
    width_in_pixels = round(Int, width_in_inches * 96)
    height_in_pixels = round(Int, height_in_inches * 96)

    fig = Figure(size = (width_in_pixels, height_in_pixels))
    ax1 = Axis(fig[1, 1], xscale=Makie.log10, yscale=Makie.log10,
               xlabel="N", ylabel="Relative Error", title=label1)
    ax2 = Axis(fig[1, 2], xscale=Makie.log10, yscale=Makie.log10,
               xlabel="N", ylabel="Relative Error", title=label2)
    linkyaxes!(ax1, ax2)

    scatterlines!(ax1, scenario_data1.spectral.N, scenario_data1.spectral.accuracy, color=:blue, label="SC - PDE15s")
    scatterlines!(ax1, scenario_data1.sparse.N, scenario_data1.sparse.accuracy, color=:green, label="FV - QNDF()")
    scatterlines!(ax2, scenario_data2.spectral.N, scenario_data2.spectral.accuracy, color=:blue, label="SC - PDE15s")
    scatterlines!(ax2, scenario_data2.sparse.N, scenario_data2.sparse.accuracy, color=:green, label="FV - QNDF()")

    Legend(fig[2, 1:2], ax1, "Methods", orientation=:horizontal)
    hidedecorations!(ax1, label=false, ticks=false, ticklabels=false)
    hidedecorations!(ax2, label=false, ticks=false, ticklabels=false)

    file_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\figures\\scaling\\prod\\$(filename)_N_vs_Err_plot.pdf"
    save(file_path, fig) # Save the figure to a file
    
    return fig
end

function plot_Err_vs_Time_dual(scenario_data1::scenario_plot_data, scenario_data2::scenario_plot_data, label1::String, label2::String, filename)
    CairoMakie.activate!() # Activate CairoMakie backend

    # Define the size in inches
    width_in_inches = 6.25
    height_in_inches = 4.29

    # Convert inches to pixels (assuming 96 dpi)
    width_in_pixels = round(Int, width_in_inches * 96)
    height_in_pixels = round(Int, height_in_inches * 96)

    fig = Figure(size = (width_in_pixels, height_in_pixels))
    ax1 = Axis(fig[1, 1], xscale=Makie.log10, yscale=Makie.log10,
               xlabel="Relative Error", ylabel="Time", title=label1)
    ax2 = Axis(fig[1, 2], xscale=Makie.log10, yscale=Makie.log10,
               xlabel="Relative Error", ylabel="Time", title=label2)
    linkyaxes!(ax1, ax2)

    scatterlines!(ax1, scenario_data1.spectral.accuracy, scenario_data1.spectral.time, color=:blue, label="SC - PDE15s")
    scatterlines!(ax1, scenario_data1.dense.accuracy, scenario_data1.dense.time, color=:red, label="FV - QNDF() Dense")
    scatterlines!(ax1, scenario_data1.sparse.accuracy, scenario_data1.sparse.time, color=:green, label="FV - QNDF() Sparse")
    scatterlines!(ax2, scenario_data2.spectral.accuracy, scenario_data2.spectral.time, color=:blue, label="SC - PDE15s")
    scatterlines!(ax2, scenario_data2.dense.accuracy, scenario_data2.dense.time, color=:red, label="FV - QNDF() Dense")
    scatterlines!(ax2, scenario_data2.sparse.accuracy, scenario_data2.sparse.time, color=:green, label="FV - QNDF() Sparse")

    Legend(fig[2, 1:2], ax1, "Methods", orientation=:horizontal)
    hidedecorations!(ax1, label=false, ticks=false, ticklabels=false)
    hidedecorations!(ax2, label=false, ticks=false, ticklabels=false)

    file_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\figures\\scaling\\prod\\$(filename)_Err_vs_Time_plot.pdf"
    save(file_path, fig) # Save the figure to a file
    
    return fig
end


function generate_all_figs()
    plot_data = gen_plot_data()
    for scenario_data in plot_data
        plot_N_vs_Time(scenario_data)
        plot_N_vs_Error(scenario_data)
        plot_Err_vs_Time(scenario_data)
    end
end

function generate_prod_figs()

    plot_data = gen_plot_data()
    plot_N_vs_Time_dual(plot_data[1],plot_data[2],"Low Voltage","High Voltage","lithiation")
    plot_N_vs_Time_dual(plot_data[3],plot_data[4],"Low Voltage","High Voltage","delithiation")

    plot_N_vs_Error_dual(plot_data[1],plot_data[2],"Low Voltage","High Voltage","lithiation")
    plot_N_vs_Error_dual(plot_data[3],plot_data[4],"Low Voltage","High Voltage","delithiation")

    plot_Err_vs_Time_dual(plot_data[1],plot_data[2],"Low Voltage","High Voltage","lithiation")
    plot_Err_vs_Time_dual(plot_data[3],plot_data[4],"Low Voltage","High Voltage","delithiation")

end
function plot_sol_t(low_volt_dat, high_volt_dat, low_volt_label, high_volt_label, filename)
    CairoMakie.activate!() # Activate CairoMakie backend

    # Define the size in inches
    width_in_inches = 6.25
    height_in_inches = 4.29

    # Convert inches to pixels (assuming 96 dpi)
    width_in_pixels = round(Int, width_in_inches * 96)
    height_in_pixels = round(Int, height_in_inches * 96)

    fig = Figure(size = (width_in_pixels, height_in_pixels))
    ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = L"c(x)", title = low_volt_label)
    ax2 = Axis(fig[1, 2], xlabel = "x", ylabel = L"c(x)", title = high_volt_label)

    linkaxes!(ax1, ax2)

    xgrid = LinRange(0.0, 1.0, 8192)
    
    # Define color gradient
    color_gradient = range(colorant"lightblue", stop=colorant"purple", length=5)
    # Data selection
    low_ti = low_volt_dat[:,1]
    low_q = low_volt_dat[:,25]
    low_tint = low_volt_dat[:,50]
    low_3q = low_volt_dat[:,75]
    low_tf = low_volt_dat[:,end]
    
    # Plot lines with gradient
    lines!(ax1, xgrid, low_ti, color=color_gradient[1], label=L"t=t_0")
    lines!(ax1, xgrid, low_q, color=color_gradient[2])
    lines!(ax1, xgrid, low_tint, color=color_gradient[3])
    lines!(ax1, xgrid, low_3q, color=color_gradient[4])
    lines!(ax1, xgrid, low_tf, color=color_gradient[5], label=L"t_{f}")
    
    ylims!(ax1, 0.0, 1.0)

    high_ti = high_volt_dat[:,1]
    high_q = high_volt_dat[:,25]
    high_tint = high_volt_dat[:,50]
    high_3q = high_volt_dat[:,75]
    high_tf = high_volt_dat[:,end]
    
    # Plot lines with gradient
    lines!(ax2, xgrid, high_ti, color=color_gradient[1], label=L"t=t_0")
    lines!(ax2, xgrid, high_q, color=color_gradient[2])
    lines!(ax2, xgrid, high_tint, color=color_gradient[3])
    lines!(ax2, xgrid, high_3q, color=color_gradient[4])
    lines!(ax2, xgrid, high_tf, color=color_gradient[5], label=L"t_{f}")
    
    ylims!(ax2, 0.0, 1.0)

    Legend(fig[2, 1:2], ax1, "Time", orientation=:horizontal)

    hidedecorations!(ax1, label=false, ticks=false, ticklabels=false)
    hidedecorations!(ax2, label=false, ticks=false, ticklabels=false)

    file_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\figures\\sol\\$(filename)_sol_plot.pdf"
    save(file_path, fig) # Save the figure to a file
    
    return fig
end

function generate_sol_plot()
    sol_collection = gather_ref_sol_plotdat()
    plot_sol_t(sol_collection[1],sol_collection[2],"Lithiation Low Voltage","Lithiation High Voltage","lithiation")
    plot_sol_t(sol_collection[3],sol_collection[4],"Delithiation Low Voltage","Delithiation High Voltage","delithiation")
end