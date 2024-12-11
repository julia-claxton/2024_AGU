using Statistics
using LinearAlgebra
using Plots, Plots.PlotMeasures

# Import my ELFIN tools
include("./Julia_ELFIN_Tools/Events.jl")
include("./Julia_ELFIN_Tools/Visualization.jl")

function get_elfin_grid_bin_edges(event::Event)
    # Energy
    elfin_energy_bin_edges = cat(event.energy_bins_min[1], event.energy_bins_max, dims = 1)
        
    # Pitch angle
    distance_between_bins = diff(event.avg_pitch_angles)
    elfin_pa_bin_edges = event.avg_pitch_angles[1:end-1] .+ (distance_between_bins ./ 2)
    first_bin_edge = event.avg_pitch_angles[1] .- (distance_between_bins[1] ./ 2)
    last_bin_edge = event.avg_pitch_angles[end] .+ (distance_between_bins[end] ./ 2)
    elfin_pa_bin_edges = cat(first_bin_edge, elfin_pa_bin_edges, last_bin_edge, dims = 1)

    return elfin_energy_bin_edges, elfin_pa_bin_edges
end

function convert_elfin_grid_to_simulation_grid(event::Event; t = 1:event.n_datapoints, normalize_northern_hemisphere = true, verbose = false, cull_out_of_range = true)
    if verbose && (event.loss_cone_angles[t[1]] > 90); println("Southern Hemisphere"); end

    # Get flux and bin edges from data
    _, elfin_nflux = integrate_flux(event, time = true, time_idxs = t)
    elfin_energy_bin_edges, elfin_pa_bin_edges = get_elfin_grid_bin_edges(event)

    # Get simulation grid
    simulation_energy_nbins, simulation_energy_bin_edges, simulation_energy_bin_means, simulation_pa_nbins, simulation_pa_bin_edges, simulation_pa_bin_means, SIMULATION_α_MAX = get_simulation_bins()

    # Precalculate number of particles in and area of each ELFIN pixel to avoid redundant calculations
    ΔΩ_elfin = [2π * (cosd(elfin_pa_bin_edges[α]) - cosd(elfin_pa_bin_edges[α+1])) for α in 1:16]
    ΔE_elfin = [(elfin_energy_bin_edges[e+1] - elfin_energy_bin_edges[e]) / 1000 for e in 1:16] # !!! IN MeV NOT keV !!!
    elfin_pixel_area = [(elfin_energy_bin_edges[e+1] - elfin_energy_bin_edges[e]) * (elfin_pa_bin_edges[α+1] - elfin_pa_bin_edges[α]) for e in 1:16, α in 1:16]
    n_particles_elfin = [elfin_nflux[e, α] * ΔΩ_elfin[α] * ΔE_elfin[e] for e in 1:16, α in 1:16]

    # Go through every simulation pixel and sum contribution from every overlapping ELFIN pixel, weighted by overlap
    simulation_n_particles = zeros(simulation_energy_nbins, simulation_pa_nbins)
    for simulation_coords = CartesianIndices(simulation_n_particles) # For every simulation pixel
        simulation_e_idx = simulation_coords[1]
        simulation_α_idx = simulation_coords[2]

        # Add contributions of each intersecting ELFIN pixel to the current simulation pixel
        for elfin_coords = CartesianIndices(elfin_nflux) # For every ELFIN pixel
            elfin_e_idx = elfin_coords[1]
            elfin_α_idx = elfin_coords[2]

            overlap_area = rectangle_overlap_area([simulation_energy_bin_edges[simulation_e_idx], simulation_energy_bin_edges[simulation_e_idx+1]], [simulation_pa_bin_edges[simulation_α_idx], simulation_pa_bin_edges[simulation_α_idx+1]],
                                                [elfin_energy_bin_edges[elfin_e_idx], elfin_energy_bin_edges[elfin_e_idx+1]], [elfin_pa_bin_edges[elfin_α_idx], elfin_pa_bin_edges[elfin_α_idx+1]]
            )
            if overlap_area == 0; continue; end
            if overlap_area > elfin_pixel_area[elfin_e_idx, elfin_α_idx]; error("Simulation grid larger than ELFIN grid. This code doesn't currently support that."); end

            # Assuming uniform distribution of particles within ELFIN pixels, get number of particles in simulation pixel, then convert to flux units
            n_particles_in_sim_pixel = n_particles_elfin[elfin_e_idx, elfin_α_idx] * (overlap_area / elfin_pixel_area[elfin_e_idx, elfin_α_idx])
            
            if isnan(n_particles_in_sim_pixel)
                n_particles_in_sim_pixel = 0
            end

            # Add particles to this pixel
            simulation_n_particles[simulation_e_idx, simulation_α_idx] += n_particles_in_sim_pixel
        end
    end

    # Convert to flux units
    ΔE_simulation = [(simulation_energy_bin_edges[e+1] - simulation_energy_bin_edges[e]) / 1000 for e in 1:simulation_energy_nbins] # MeV
    ΔΩ_simulation = [2π * (cosd(simulation_pa_bin_edges[α]) - cosd(simulation_pa_bin_edges[α+1])) for α in 1:simulation_pa_nbins] # Str
    simulation_flux = [simulation_n_particles[e, α] / (ΔΩ_simulation[α] * ΔE_simulation[e]) for e in 1:simulation_energy_nbins, α in 1:simulation_pa_nbins]

    # Normalize to Northern hemisphere
    if normalize_northern_hemisphere && event.loss_cone_angles[t[1]] > 90
        if verbose; println("Normalizing input to Northern hemisphere..."); end
        reverse!(simulation_flux, dims = 2)
    end

    # Cull data outside of the loss cone if needed
    if cull_out_of_range == true
        simulation_flux[:, simulation_pa_bin_edges[begin:end-1] .> SIMULATION_α_MAX] .= 0
    end

    # Return
    return simulation_flux
end

function convert_simulation_grid_to_elfin_grid(simulation_distribution, event::Event)
    # Get bin edges from data
    elfin_energy_bin_edges, elfin_pa_bin_edges = get_elfin_grid_bin_edges(event)

    # Get simulation grid
    simulation_energy_nbins, simulation_energy_bin_edges, simulation_energy_bin_means, simulation_pa_nbins, simulation_pa_bin_edges, simulation_pa_bin_means, SIMULATION_α_MAX = get_simulation_bins()

    # Go through every simulation pixel and sum contribution to every overlapping ELFIN pixel, weighted by overlap
    elfin_distribution = zeros(16,16)
    for simulation_coords = CartesianIndices(simulation_distribution) # For every simulation pixel
        simulation_e_idx = simulation_coords[1]
        simulation_α_idx = simulation_coords[2]

        simulation_pixel_area = (simulation_energy_bin_edges[simulation_e_idx+1] - simulation_energy_bin_edges[simulation_e_idx]) * (simulation_pa_bin_edges[simulation_α_idx+1] - simulation_pa_bin_edges[simulation_α_idx]) # Calculate outside the loop to avoid wasting performance
        ΔΩ_simulation_pixel = 2π * (cosd(simulation_pa_bin_edges[simulation_α_idx]) - cosd(simulation_pa_bin_edges[simulation_α_idx+1]))
        ΔE_simulation_pixel = (simulation_energy_bin_edges[simulation_e_idx+1] - simulation_energy_bin_edges[simulation_e_idx]) / 1000 # !!! IN MeV NOT keV !!!

        # Add contributions of each intersecting ELFIN pixel to the current simulation pixel
        for elfin_coords = CartesianIndices(elfin_distribution) # For every ELFIN pixel
            elfin_e_idx = elfin_coords[1]
            elfin_α_idx = elfin_coords[2]

            overlap_area = rectangle_overlap_area([simulation_energy_bin_edges[simulation_e_idx], simulation_energy_bin_edges[simulation_e_idx+1]], [simulation_pa_bin_edges[simulation_α_idx], simulation_pa_bin_edges[simulation_α_idx+1]],
                                                [elfin_energy_bin_edges[elfin_e_idx], elfin_energy_bin_edges[elfin_e_idx+1]], [elfin_pa_bin_edges[elfin_α_idx], elfin_pa_bin_edges[elfin_α_idx+1]]
            )
            overlap_fraction = overlap_area / simulation_pixel_area
            if overlap_fraction == 0; continue; end

            # Add the number of particles in each overlapping simulation pixel to ELFIN pixel, weighted by overlap. We convert to flux units later
            n_particles_in_sim_pixel = simulation_distribution[simulation_e_idx, simulation_α_idx] * ΔΩ_simulation_pixel * ΔE_simulation_pixel
            elfin_distribution[elfin_e_idx, elfin_α_idx] += overlap_fraction * n_particles_in_sim_pixel
        end
    end
    # Convert ELFIN distribution from number of particles to flux units
    ΔΩ_elfin = [2π * (cosd(elfin_pa_bin_edges[α]) - cosd(elfin_pa_bin_edges[α+1])) for α in 1:16]
    ΔE_elfin = [(elfin_energy_bin_edges[e+1] - elfin_energy_bin_edges[e]) / 1000 for e in 1:16] # !!! IN MeV NOT keV !!!
    elfin_distribution = [elfin_distribution[e,α] / (ΔΩ_elfin[α] * ΔE_elfin[e]) for e in 1:16, α in 1:16]

    return elfin_distribution
end

function rectangle_overlap_area(rect1x, rect1y, rect2x, rect2y)
    @assert length(rect1x) == length(rect1y) == length(rect2x) == length(rect2y) == 2 "Rectangle dimensions incorrect"
    
    Δx = min(rect1x[2], rect2x[2]) - max(rect1x[1], rect2x[1])
    Δy = min(rect1y[2], rect2y[2]) - max(rect1y[1], rect2y[1])

    # Clamp to zero – negative values mean no overlap
    if (Δx < 0) || (Δy < 0)
        Δx = Δy = 0
    end

    return Δx * Δy
end

function compare_sim_to_data(event::Event, distributions; t = 1:event.n_datapoints, black_out_trapped = false)
    if length(size(distributions)) == 2
        input = copy(distributions)
        distributions = zeros(1, size(input)...)
        distributions[1,:,:] = input
    end

    sim_energy_nbins, sim_energy_bin_edges, sim_energy_bin_means, sim_pa_nbins, sim_pa_bin_edges, sim_pa_bin_means, SIMULATION_α_MAX = get_simulation_bins()

    # Data
    elfin_nflux = convert_elfin_grid_to_simulation_grid(event, t = t, cull_out_of_range = false)
    elfin_eflux = [elfin_nflux[e,α] * sim_energy_bin_means[e] for e in 1:sim_energy_nbins, α in 1:sim_pa_nbins]

    # Nothern hemisphere normalization
    α_lc = event.loss_cone_angles[t[1]]
    α_alc = event.anti_loss_cone_angles[t[1]]
    if α_lc > 90 # Southern hemisphere
        α_lc = 180 - α_lc
        α_alc = 180 - α_alc
    end

    to_plot = log10.(elfin_eflux)

    xlims = (0, 180)
    Δx = xlims[2] - xlims[1]

    ylims = log10.([event.energy_bins_min[1], 1e4])
    Δy = ylims[2] - ylims[1]
    
    clims = (max(to_plot...)-3, max(to_plot...))

    heatmap(sim_pa_bin_edges, log10.(sim_energy_bin_edges), to_plot,
        title = "Data",

        xlabel = "Pitch Angle, deg",
        xlims = xlims,

        ylabel = "Energy, keV",
        ylims = ylims,
        yticks = ([2, 3, 4], ["10²", "10³", "10⁴"]),

        colorbar_title = "Log10 Flux\n keV/(s cm2 str MeV)",
        clims = clims,
        colormap = :haline,

        aspect_ratio = Δx/Δy,

        topmargin = -5mm,
        bottommargin = -5mm,
        leftmargin = 5mm,
        rightmargin = 5mm,

        background_color_inside = :black,
        bordercolor = :transparent, # No axis border
        foreground_color_axis = :transparent, # No ticks
        framestyle = :box,
        grid = false
    )
    if black_out_trapped == true
        plot!(Shape([α_lc, α_alc, α_alc, α_lc], [ylims[1]+.015, ylims[1]+.015, ylims[2]-.015, ylims[2]-.015]), # Blocks trapped region for easier comparison
            fillcolor = :black,
            linecolor = :black,
            label = false
        )
    end
    vline!([α_lc],
        linecolor = :white,
        linewidth = 1.2,
        linestyle = :solid,
        label = false
    )
    vline!([α_alc],
        linecolor = :white,
        linewidth = 1.2,
        linestyle = :dash,
        label = false
    )
    data = plot!()

    # Simulation
    simulation_nflux = dropdims(sum(distributions, dims = 1), dims = 1)

    # Clip all simulation pixels below the minimum value from data. Like a noisegate to make comparison easier
    noisegate = min(elfin_nflux[elfin_nflux .≠ 0]...)
    simulation_nflux[simulation_nflux .< noisegate] .= 0

    # Convert to energy flux to see things happening at higher energy channels
    simulation_eflux = [simulation_nflux[e,α] * sim_energy_bin_means[e] for e in 1:sim_energy_nbins, α in 1:sim_pa_nbins]

    to_plot = log10.(simulation_eflux)
    heatmap(sim_pa_bin_edges, log10.(sim_energy_bin_edges), to_plot,
        title = "Simulation",

        xlabel = "Pitch Angle, deg",
        xlims = xlims,

        ylabel = "Energy, keV",
        ylims = ylims,
        yticks = ([2, 3, 4], ["10²", "10³", "10⁴"]),

        colorbar_title = "Log10 Flux\nkeV/(s cm2 str MeV)",
        clims = clims,
        colormap = :haline,

        aspect_ratio = Δx/Δy,

        topmargin = -5mm,
        bottommargin = -5mm,
        leftmargin = 5mm,
        rightmargin = 5mm,

        background_color_inside = :black,
        bordercolor = :transparent, # No axis border
        foreground_color_axis = :transparent, # No ticks
        framestyle = :box,
        grid = false
    )
    if black_out_trapped == true
        plot!(Shape([α_lc, α_alc, α_alc, α_lc], [ylims[1]+.015, ylims[1]+.015, ylims[2]-.015, ylims[2]-.015]), # Blocks trapped region for easier comparison
            fillcolor = :black,
            linecolor = :black,
            label = false
        )
    end
    vline!([α_lc],
        linecolor = :white,
        linewidth = 1.2,
        linestyle = :solid,
        label = false
    )
    vline!([α_alc],
        linecolor = :white,
        linewidth = 1.2,
        linestyle = :dash,
        label = false
    )
    simulation = plot!()

    # Get residuals (ratio)
    residuals = simulation_nflux ./ elfin_nflux
    to_plot = log10.(residuals)

    heatmap(sim_pa_bin_edges, log10.(sim_energy_bin_edges), to_plot,
        title = "Residuals\n",

        xlabel = "Pitch Angle, deg",
        xlims = xlims,

        ylabel = "Energy, eV",
        ylims = ylims,
        yticks = ([2, 3, 4], ["10²", "10³", "10⁴"]),
        
        colorbar_title = "\nLog10(Sim/Data)",
        colormap =  :diverging_bwr_55_98_c37_n256, #:diverging_cwm_80_100_c22_n256
        clims = (-1, 1),

        aspect_ratio = Δx/Δy,

        margin = -5mm,

        grid = false,
        framestyle = :box,
        background_color_inside = RGB(0xD0D0D0)
    )
    plot!(Shape([α_lc, α_alc, α_alc, α_lc], [ylims[1]+.015, ylims[1]+.015, ylims[2]-.015, ylims[2]-.015]), # Blocks trapped region for easier comparison
        fillcolor = RGB(0xD0D0D0),
        linecolor = RGB(0xD0D0D0),
        label = false
    )
    vline!([α_lc],
        linecolor = :black,
        linewidth = 1.2,
        linestyle = :solid,
        label = false
    )
    vline!([α_alc],
        linecolor = :black,
        linewidth = 1.2,
        linestyle = :dash,
        label = false
    )
    residual_plot = plot!()

    #=
    histogram_bin_edges = -1:.1:1
    alc_idxs = findall(event.avg_pitch_angles .≥ α_alc)
    alc_residual_frequencies = exact_1dhistogram(log10.(residuals)[:, alc_idxs][:], histogram_bin_edges)
    alc_residual_frequencies = cat(alc_residual_frequencies, alc_residual_frequencies[end], dims = 1)
    plot(histogram_bin_edges, alc_residual_frequencies,
        title = "ALC Residuals",
        linetype = :steppost,
        label = false,
        linecolor = :black,
        linewidth = 1.5,

        xlabel = "Log10(sim/data)",
        xlims = (-1, 1),

        ylabel = "Frequency",
        ylims = (0, 20),

        aspect_ratio = (2/20),
        framestyle = :box
    )
    hist = plot!()
    #println("Mean ALC sim/data = $(mean(filter(isfinite, residuals[:, alc_idxs])))")
    =#

    layout = @layout [a b;
                      c{.4h}]

    plot(data, simulation, residual_plot,
        background = :transparent,
        layout = layout,
        size = (1, .8) .* 700,
        dpi = 300,
    )
    display("image/png", plot!())
end