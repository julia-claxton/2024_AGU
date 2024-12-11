const TOP_LEVEL = dirname(@__DIR__)
@assert endswith(TOP_LEVEL, "2024_AGU_Multibounce") "Script in wrong directory! TOP_LEVEL = $TOP_LEVEL"

using Statistics
using LinearAlgebra
using Plots, Plots.PlotMeasures
using NPZ

include("./General_Functions.jl") # Provides general-purpose functions I find useful

# ---------------- Functions ----------------
function figure_event_detection()
    data = npzread("$TOP_LEVEL/data/event_detection_example.npz")
    precipitation = data["precipitation"]
    score = data["score"]
    event_mask = data["event_mask"]
    energy_bins_mean = data["energy_bins_mean"]

    # Heatmap
    yticks = Int.(round.(energy_bins_mean))
    yticks = yticks[1:3:16]
    heatmap(log10.(precipitation'),
        xlabel = "",
        xlims = (.5, size(precipitation)[1] + .5),

        ylabel = "Energy, keV",
        yticks = (1:3:16, yticks),
        ylims = (.5, 16.5),

        colorbar_title = "\nLog10 Jprec/Jtrap",
        colormap = cgrad(:RdPu, rev = true),
        colorbar = true,
        clims = (-1.25, .25),

        rightmargin = 10mm,
        leftmargin = 5mm,
        topmargin = -30mm,

        background_color_inside = RGB(.85, .8, .8),
        aspect_ratio = (size(precipitation)[1]/size(precipitation)[2]) * .15,
        framestyle = :box,
        tickdirection = :out
    )
    p1 = plot!()

    # Precipitation score
    plot(
        title = "Precipitation Score",

        xlabel = "Data Index",
        xlims = (.5, size(precipitation)[1] + .5),

        ylabel = "Score (arb. units)",
        ylims = (0, 4),

        leftmargin = 11mm,
        rightmargin = 6.5mm,

        aspect_ratio = (size(precipitation)[1]/4) * .15,
        framestyle = :box,
        tickdirection = :out
    )
    hline!([2],
        label = false,
        linestyle = :dash,
        linecolor = :grey
    )
    plot!(score,
        label = false,
        linewidth = 1.5,
        linecolor = RGB(0xc23aa7)
    )
    p2 = plot!()

    # Detection area
    plot(1 .* event_mask,
        fill = true,
        label = false,
        linetype = :steppre,
        fillcolor = :red,
        linecolor = :transparent,

        xlims = (.5, size(precipitation)[1] + .5),

        ylims = (0, 1),

        topmargin = 5mm,
        bottommargin = -10mm,

        framestyle = :none,
        grid = false,
    )
    p3 = plot!()

    spacer = plot(framestyle = :none,
                  axis = false,
    )

    layout = @layout [a{.25h} b{.01w}; c; d e{.01w}]
    plot(p3, spacer, p1, p2, spacer,
        suptitle = "2021-03-06 07:02:00 ELB",

        size = (2, 1) .* 400,
        layout = layout,
        background = :transparent,
        dpi = 350
    )
    png("$(TOP_LEVEL)/results/figures/event_detection.png")
    display(plot!())
end

function figure_beam_fitting()
    data = npzread("$(TOP_LEVEL)/data/beam_fitting.npz")
    e_edges = data["e_edges"]
    pa_edges = data["pa_edges"]
    data_nflux = data["data_nflux"]
    beam_coordinates_old = data["beam_coordinates"]
    beam_colors = data["beam_colors"]
    fit_energy_fraction = data["fit_energy_fraction"]
    fit_number_fraction = data["fit_number_fraction"]
    
    beam_coordinates = [(beam_coordinates_old[i,1], beam_coordinates_old[i,2]) for i in 1:size(beam_coordinates_old)[1]]

    scatter(beam_coordinates,
        title = "Event Fit\n",

        xlabel = "Pitch Angle, deg",
        xlims = (0, 80),

        ylabel = "Energy, keV",
        ylims = (40, 10e3),
        yscale = :log10,

        colorbar_title = "Beam Strength, # e-",
        zcolor = beam_colors,
        colormap = cgrad(:ice, rev = false),
        clims = (1, 7),


        label = false,
        markerstrokewidth = 0,
        markersize = 5,
        background_color_inside = :black,
        background_color_outside = :transparent,
        framestyle = :box,
        grid = false,

        leftmargin = 5mm,
        topmargin = 5mm,
        size = (1,1) .* 400,
        aspect_ratio = (80/(10e3-40)) * 1.2,
        dpi = 350
    )
    annotate!((40, 1.5e4, text("2021-03-06 07:02:00 ELB", :black, :center, 11)))
    p1 = plot!()

    heatmap(pa_edges, log10.(e_edges), log10.(data_nflux),
        title = "Event Data\n",

        xlabel = "Pitch Angle, deg",
        xlims = (0, 80),

        ylabel = "Energy, keV",
        ylims = log10.([50, 7e3]),

        colorbar_title = "ELFIN Fluence, #/(cm² str MeV)",
        colormap = cgrad(:RdPu, rev = true),
        clims = (2, 8),

        background_color = :transparent,
        framestyle = :box,

        leftmargin = 5mm,
        topmargin = 5mm,
        aspect_ratio = (80/log10.(7e3/50)) * 1.2,
        size = (1,1) .* 400,
        dpi = 350
    )
    annotate!((40, log10(1e4), text("2021-03-06 07:02:00 ELB", :black, :center, 11)))
    p2 = plot!()

    plot(p2, p1,
        layout = (1,2),
        size = (2, .92) .* 400,
        dpi = 300
    )
    png("$(TOP_LEVEL)/results/figures/beam_fitting.png")
    display(plot!())
end

function figure_number_residuals()
    # Number residuals
    results = npzread("$TOP_LEVEL/data/elfin_lifetime_simulation_statistics.npz")
    full_n_residuals = results["full_range_number_residuals"]
    elfin_n_residuals = results["elfin_range_number_residuals"]

    histogram_bin_edges = 10 .^ (-2:.05:2)
    full_frequencies = exact_1dhistogram(full_n_residuals, histogram_bin_edges)
    elfin_frequencies = exact_1dhistogram(elfin_n_residuals, histogram_bin_edges)
    
    plot(histogram_bin_edges, cat(full_frequencies, full_frequencies[end], dims = 1),
        title = "Total Backscattered Electrons",
        
        linetype = :steppost,
        label = "10 keV - 10 MeV",
        linecolor = :black,
        linealpha = .5,
        linewidth = 1.5,

        fill = true,
        fillcolor = RGB(0xfcc0ec),
        fillalpha = .5,

        xlabel = "Simulation/Data",
        xscale = :log10,
        xminorticks = true,
        xlims = (.01, 100),

        ylabel = "Frequency",
        ylims = (0, max(full_frequencies...)*1.1),

        background = :transparent,
        tickdirection = :out,
        framestyle = :box,
        grid = false,
        size = (1,1) .* 500,
        dpi = 300
    )
    plot!(histogram_bin_edges, cat(elfin_frequencies, elfin_frequencies[end], dims = 1),        
        linetype = :steppost,
        label = "50 keV - 7.2 MeV",
        linecolor = :black,
        linealpha = .5,
        linewidth = 1.5,

        fill = true,
        fillcolor = RGB(0xd1e9ff),
        fillalpha = .5,

        ylims = (0, max(elfin_frequencies...)*1.1),

        aspect_ratio = (100-.01)/(max(elfin_frequencies...)*1.1),
    )
    # Grid
    vline!(10. .^(-2:2),
        label = false,
        linecolor = :black,
        linealpha = .1,
        linewidth = .5
    )
    png("$(TOP_LEVEL)/results/figures/number_residuals.png")
    display(plot!())
end

function figure_energy_residuals()
    # Energy residuals
    results = npzread("$TOP_LEVEL/data/elfin_lifetime_simulation_statistics.npz")
    full_e_residuals = results["full_range_energy_residuals"]
    elfin_e_residuals = results["elfin_range_energy_residuals"]

    histogram_bin_edges = 10 .^ (-2:.05:2)
    full_frequencies = exact_1dhistogram(full_e_residuals, histogram_bin_edges)
    elfin_frequencies = exact_1dhistogram(elfin_e_residuals, histogram_bin_edges)
    
    plot(histogram_bin_edges, cat(full_frequencies, full_frequencies[end], dims = 1),
        title = "Total Backscattered Energy",
        
        linetype = :steppost,
        label = "10 keV - 10 MeV",
        linecolor = :black,
        linealpha = .5,
        linewidth = 1.5,

        fill = true,
        fillcolor = RGB(0xfcc0ec),
        fillalpha = .5,

        xlabel = "Simulation/Data",
        xscale = :log10,
        xminorticks = true,
        xlims = (.01, 100),

        ylabel = "Frequency",
        ylims = (0, max(full_frequencies...)*1.1),

        background = :transparent,
        aspect_ratio = (100-.01)/(max(full_frequencies...)*1.1),
        tickdirection = :out,
        framestyle = :box,
        grid = false,
        size = (1,1) .* 500,
        dpi = 300
    )
    plot!(histogram_bin_edges, cat(elfin_frequencies, elfin_frequencies[end], dims = 1),        
        linetype = :steppost,
        label = "50 keV - 7.2 MeV",
        linecolor = :black,
        linealpha = .5,
        linewidth = 1.5,

        fill = true,
        fillcolor = RGB(0xd1e9ff),
        fillalpha = .5,

        ylims = (0, max(elfin_frequencies...)*1.1),

        aspect_ratio = (100-.01)/(max(elfin_frequencies...)*1.1),
    )
    # Grid
    vline!(10. .^(-2:2),
        label = false,
        linecolor = :black,
        linealpha = .1,
        linewidth = .5
    )
    png("$(TOP_LEVEL)/results/figures/energy_residuals.png")
    display(plot!())
end

function figure_residual_distribution()
    results = npzread("$TOP_LEVEL/data/elfin_lifetime_simulation_statistics.npz")
    e_edges = results["energy_bin_edges"]
    pa_edges = results["pitch_angle_bin_edges"]
    residual_distribution = results["residual_distribution"]

    mean_residual(x) = 10^mean(log10.(x[isfinite.(x)]))
    residuals = [mean_residual(residual_distribution[e,pa,:]) for e in 1:length(e_edges)-1, pa in 1:length(pa_edges)-1]

    xlim = [110, 180]
    ylim = log10.([50, 1e4])
    heatmap(pa_edges, log10.(e_edges), log10.(residuals),
        title = "Simulation Flux Residuals",

        xlabel = "NH Pitch Angle, deg",
        xlims = xlim,
        xticks = 110:10:180,

        ylabel = "Energy, keV",
        yticks = ([2,3,4], ["10²", "10³", "10⁴"]),
        ylims = ylim,

        colorbar_title = "\nLog10 Mean(Simulated Flux/Recorded Flux)",
        colormap = cgrad(:RdBu_7, rev = true),
        clims = (-1, 1),

        leftmargin = -5mm,

        aspect_ratio = (diff(xlim)/diff(ylim))[1] * 1.2,
        size = (1.2, 1) .* 400,
        framestyle = :box,
        tickdirection = :out,
        grid = false,
        background = :transparent,
        dpi = 300
    )
    png("$(TOP_LEVEL)/results/figures/residual_distribution.png")
    display(plot!())
end

function figure_lc_mb_fraction()
    # Number loss cone multibounce fraction
    results = npzread("$TOP_LEVEL/data/elfin_lifetime_simulation_statistics.npz")
    n_lc_multibounce_fraction = results["number_lc_multibounce_fraction"]

    histogram_bin_edges = 0:.0025:.2
    to_plot = n_lc_multibounce_fraction[n_lc_multibounce_fraction .< .99] # Hopefully won't need this later
    frequencies = exact_1dhistogram(to_plot, histogram_bin_edges)
    plot(histogram_bin_edges, cat(frequencies, frequencies[end], dims = 1),
        title = "Loss Cone Multibounce Fraction",
        
        linetype = :steppost,
        label = false,
        linecolor = :black,
        linewidth = 1.5,

        fill = true,
        fillcolor = RGB(0xddd8ed),

        xlabel = "Multibounce Fraction",
        xlims = (0, .15),
        xticks = 0:.05:.15,
        xminorticks = true,

        ylabel = "Frequency",
        ylims = (0, max(frequencies...)*1.1),

        background = :transparent,
        aspect_ratio = (.15)/(max(frequencies...)*1.1),
        tickdirection = :out,
        framestyle = :box,
        grid = false,
        size = (1,1) .* 500,
        dpi = 300
    )
    # Grid
    vline!(0:.01:.15,
        label = false,
        linecolor = :black,
        linealpha = .1,
        linewidth = .5
    )
    # Points of interest
    μ = mean(to_plot)
    vline!([μ],
        linecolor = :grey,
        linestyle = :dash,
        linewidth = 2,
        label = "μ"
    )
    σ = std(to_plot)
    vline!([μ+3σ, μ-3σ],
        linecolor = :grey,
        linestyle = :dashdot,
        linealpha = .8,
        linewidth = 1.5,
        label = "±3σ"
    )
    png("$(TOP_LEVEL)/results/figures/number_lc_multibounce_fraction.png")
    display(plot!())
end

function figure_alc_mb_fraction()
    # Number anti loss cone multibounce fraction
    results = npzread("$TOP_LEVEL/data/elfin_lifetime_simulation_statistics.npz")
    n_lc_multibounce_fraction = results["number_lc_multibounce_fraction"]
    n_alc_multibounce_fraction = results["number_alc_multibounce_fraction"]

    histogram_bin_edges = LinRange(0, .1, 200)
    to_plot = n_alc_multibounce_fraction[n_lc_multibounce_fraction .< .99] # Hopefully won't need this later
    frequencies = exact_1dhistogram(to_plot, histogram_bin_edges)
    plot(histogram_bin_edges, cat(frequencies, frequencies[end], dims = 1),
        title = "Anti Loss Cone Multibounce Fraction",
        
        linetype = :steppost,
        label = false,
        linecolor = :black,
        linewidth = 1.5,

        fill = true,
        fillcolor = RGB(0xddd8ed),

        xlabel = "Multibounce Fraction",
        xlims = (0, .05),
        xticks = 0:.01:.05,
        xminorticks = true,

        ylabel = "Frequency",
        ylims = (0, max(frequencies...)*1.1),

        background = :transparent,
        aspect_ratio = (.05)/(max(frequencies...)*1.1),
        tickdirection = :out,
        framestyle = :box,
        grid = false,
        size = (1,1) .* 500,
        dpi = 300
    )
    # Grid
    vline!(0:.01:.15,
        label = false,
        linecolor = :black,
        linealpha = .1,
        linewidth = .5
    )
    # Points of interest
    μ = mean(to_plot)
    vline!([μ],
        linecolor = :grey,
        linestyle = :dash,
        linewidth = 2,
        label = "μ"
    )
    σ = std(to_plot)
    vline!([μ+3σ, μ-3σ],
        linecolor = :grey,
        linestyle = :dashdot,
        linealpha = .8,
        linewidth = 1.5,
        label = "±3σ"
    )
    png("$(TOP_LEVEL)/results/figures/number_alc_multibounce_fraction.png")
    display(plot!())
end

function figure_backscatter_fraction()
    # Loss rate figure
    results = npzread("$TOP_LEVEL/data/loss_cone_simulation_statistics.npz")
    backscatter_fraction = results["backscatter_fraction"]
    e_bin_edges = results["e_bin_edges"]
    pa_bin_edges = results["pa_bin_edges"]
    heatmap(pa_bin_edges, log10.(e_bin_edges), backscatter_fraction .* 100,
        title = "Backscatter Fraction",

        xlabel = "Pitch Angle, deg",
        xlims = (0, 70),

        ylabel = "Energy, keV",
        ylims = (1, 4),
        yticks = ([1, 2, 3, 4], ["10¹", "10²", "10³", "10⁴"]),

        colorbar_title = "\nBackscatter Fraction, %/bounce",
        colormap = cgrad(:RdPu, rev = true),
        clims = (0, 20),

        aspect_ratio = (180/3),
        background_color_inside = :grey,
        background_color_outside = :transparent,
        grid = false,
        framestyle = :box,
        size = (1.5, 2.5) .* 200,
        dpi = 300
    )
    png("$(TOP_LEVEL)/results/figures/backscatter_fraction.png")
    display(plot!())
end

function figure_multibounce_fraction()
    # Loss rate figure
    results = npzread("$TOP_LEVEL/data/loss_cone_simulation_statistics.npz")
    multibounce_fraction = results["multibounce_fraction"]
    e_bin_edges = results["e_bin_edges"]
    pa_bin_edges = results["pa_bin_edges"]
    heatmap(pa_bin_edges, log10.(e_bin_edges), multibounce_fraction .* 100,
        title = "Multibounce Fraction",

        xlabel = "Pitch Angle, deg",
        xlims = (0, 70),

        ylabel = "Energy, keV",
        ylims = (1, 4),
        yticks = ([1, 2, 3, 4], ["10¹", "10²", "10³", "10⁴"]),

        colorbar_title = "Multibounce Fraction, %",
        colormap = cgrad(:RdPu, rev = true),
        clims = (0, 6),

        aspect_ratio = (180/3),
        background_color_inside = :grey,
        background_color_outside = :transparent,
        grid = false,
        framestyle = :box,
        size = (1.5, 2.5) .* 200,
        dpi = 300
    )
    png("$(TOP_LEVEL)/results/figures/multibounce_fraction.png")
    display(plot!())
end

function figure_survival_fraction()
    # Loss rate figure
    results = npzread("$TOP_LEVEL/data/loss_cone_simulation_statistics.npz")
    survival_fraction = results["survival_fraction"]
    e_bin_edges = results["e_bin_edges"]
    pa_bin_edges = results["pa_bin_edges"]
    heatmap(pa_bin_edges, log10.(e_bin_edges), survival_fraction .* 100,
        title = "Survival Fraction",

        xlabel = "Pitch Angle, deg",
        xlims = (0, 70),

        ylabel = "Energy, keV",
        ylims = (1, 4),
        yticks = ([1, 2, 3, 4], ["10¹", "10²", "10³", "10⁴"]),

        colorbar_title = "\n10-bounce Survival Fraction, %",
        colormap = cgrad(:RdPu, rev = true),
        clims = (0, 4),

        aspect_ratio = (180/3),
        background_color_inside = :grey,
        background_color_outside = :transparent,
        grid = false,
        framestyle = :box,
        size = (1.5, 2.5) .* 200,
        dpi = 300
    )
    png("$(TOP_LEVEL)/results/figures/survival_fraction.png")
    display(plot!())
end


# ---------------- Script ----------------
println("========= Script Begin =========")

figure_event_detection()
figure_beam_fitting()
figure_number_residuals()
figure_energy_residuals()
figure_residual_distribution()
figure_lc_mb_fraction()
figure_alc_mb_fraction()
figure_backscatter_fraction()
figure_multibounce_fraction()
figure_survival_fraction()

println("========= Script End =========")