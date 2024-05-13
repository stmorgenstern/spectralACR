function plotres(grid,fvals)
    fig = Figure(size=(600, 500))
    ax = Axis(fig[1,1], xlabel = "x", ylabel = L"u(x)",title="Plot of Solution")#,xticks=0:100:N)
    lines!(ax,grid,fvals)
    ylims!(ax,0.0,1.0)
    return fig
end
function plot_anim_sol(xgrid, sol;cc_flag = false, title="sol")
    CairoMakie.activate!() # Activate CairoMakie backend
    
    # Ntime, Nspace = size(sol)
    Ntime = length(sol.t);
    Nspace = length(sol.u[1])

    soldat = Float32.(reduce(vcat,transpose.(sol.u)))
    if cc_flag
        soldat = soldat[:,1:end-1]
    end
    fig = Figure(size=(1200, 600))
    
    ax1 = Axis(fig[1, 1], xlabel="x", ylabel=L"c(x)", title="Concentration")
    # ax2 = Axis(fig[2, 1], xlabel="x", ylabel="Phase", title="Phase")

    # Line objects to update in loop
    idx = Observable(1);
    cdat = @lift(soldat[$idx, :])
    # ph = @lift(angle.(sol[$idx, :]))
    # ph = @lift(getangle.(sol[$idx, :]))

    amplitude_line = lines!(ax1, xgrid, cdat, color=:blue)
    ylims!(ax1,0.0,1.0)

    timestamps = 1:Ntime
    # # Loop over time steps and update plots
    record(fig, title*".mp4", timestamps; framerate=120) do frame
        
        idx[] = frame
        # autolimits!(ax1) # update limits
        # autolimits!(ax2) # update limits
        # Update plot titles with the current time step
        # ax1.title = @lift("Amplitude Squared at index k = $(10*($idx-1))")
        # ax2.title = @lift("Phase at index k = $(10*($idx-1))")
    end

    return fig
end
function plot_anim_sol_mat(xgrid, sol;cc_flag = false, title="sol")
    CairoMakie.activate!() # Activate CairoMakie backend
    
    Ntime, Nspace = size(sol)
    # Ntime = length(sol.t);
    # Nspace = length(sol.u[1])

    soldat = Float32.(sol)
    if cc_flag
        soldat = soldat[:,1:end-1]
    end
    fig = Figure(size=(1200, 600))
    
    ax1 = Axis(fig[1, 1], xlabel="x", ylabel=L"c(x)", title="Concentration")
    # ax2 = Axis(fig[2, 1], xlabel="x", ylabel="Phase", title="Phase")

    # Line objects to update in loop
    idx = Observable(1);
    cdat = @lift(soldat[$idx, :])
    # ph = @lift(angle.(sol[$idx, :]))
    # ph = @lift(getangle.(sol[$idx, :]))

    amplitude_line = lines!(ax1, xgrid, cdat, color=:blue)
    ylims!(ax1,0.0,1.0)

    timestamps = 1:Ntime
    # # Loop over time steps and update plots
    record(fig, title*".mp4", timestamps; framerate=120) do frame
        
        idx[] = frame
        # autolimits!(ax1) # update limits
        # autolimits!(ax2) # update limits
        # Update plot titles with the current time step
        # ax1.title = @lift("Amplitude Squared at index k = $(10*($idx-1))")
        # ax2.title = @lift("Phase at index k = $(10*($idx-1))")
    end

    return fig
end

function plot_time_step_hist(t;Nbins=100,title="test_sol")
    fig = Figure(size=(600, 500))
    ax = Axis(fig[1,1], xlabel = "t", ylabel = L"N",title="Timestep")#,xticks=0:100:N)
    hist!(ax,t,bins=Nbins)
    hidedecorations!(ax, label=false, ticks=false, ticklabels=false)
    file_path = "C:\\Users\\Sam\\School\\HP_ACR_1D\\SparseACR_FV\\figures\\"*title*"_t_hist.pdf"
    save(file_path, fig) # Save the figure to a file
    return fig
end