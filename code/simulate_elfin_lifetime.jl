const TOP_LEVEL = dirname(@__DIR__)
@assert endswith(TOP_LEVEL, "2024_AGU_Multibounce") "Script in wrong directory! TOP_LEVEL = $TOP_LEVEL"

# Library includes
using Statistics
using LinearAlgebra
using NPZ
using BenchmarkTools, Profile, TickTock

# Import my ELFIN tools
include("./Julia_ELFIN_Tools/Events.jl")
include("./Julia_ELFIN_Tools/Visualization.jl")

include("./EPPBackscatterSimulation/BackscatterSimulation.jl") # Provides backscatter simulation capability
include("./ELFIN_Backscatter_Tools.jl") # Provides interface between ELFIN data and backscatter simulation
include("./General_Functions.jl") # Provides general-purpose functions I find useful

# ---------------- Functions ----------------
function write_example_event_detection_data()
    event = create_event(DateTime("2021-03-06T07:02:00"), DateTime("2021-03-06T07:08:00"), "b")
    score = _precipitation_score(event, 1:event.n_datapoints)
    threshold = 2

    event_mask = score .> threshold
    event_mask[end] = 0 # If the last datapoint is high, call it the end of an event. If it's low, no impact.
    event_edges = diff(event_mask)
    if event_edges[1] == -1
        event_edges[1] = 0 # If first datapoint is high and second low, discard detection
    else
        event_edges[1] = event_mask[1] # If the first datapoint is high, call it the start of an event
    end

    npzwrite("$(TOP_LEVEL)/data/event_detection_example.npz",
        Dict("precipitation" => event.Jprec_over_Jtrap,
             "score" => score,
             "event_mask" => event_mask,
             "energy_bins_mean" => event.energy_bins_mean
        )
    )
    return
end

function calculate_elfin_simulation_stats()
    println("Detecting events in ELFIN lifetime and simulating...")

    # Get ELFIN data
    PROJECT_TOP_LEVEL = dirname(@__DIR__)
    @assert endswith(PROJECT_TOP_LEVEL, "Multibounce_Particle_Simulation") "File not in /Multibounce_Particle_Simulation. Current path: $PROJECT_TOP_LEVEL"

    elfin_data_directory = "$(PROJECT_TOP_LEVEL)/code/Julia_ELFIN_Tools/data/processed_scientific_data"
    datafiles = glob("*.npz", elfin_data_directory)
    datafiles = replace.(datafiles, "$(elfin_data_directory)/" => "")

    satellites = string.([datafiles[i][end-4] for i = eachindex(datafiles)])

    days_with_data_available = replace.(datafiles, 
        "_ela.npz" => "",
        "_elb.npz" => "",
    )
    days_with_data_available = Date.(days_with_data_available, "yyyymmdd")

    # Get simulation grid parameters
    simulation_energy_nbins, simulation_energy_bin_edges, simulation_energy_bin_means, simulation_pa_nbins, simulation_pa_bin_edges, simulation_pa_bin_means, SIMULATION_α_MAX = get_simulation_bins()

    # Detect and simulate event over ELFIN's liftime, store information for each.
    elfin_range_energy_residuals = []
    full_range_energy_residuals = []

    elfin_range_number_residuals = []
    full_range_number_residuals = []

    energy_lc_multibounce_fraction = []
    number_lc_multibounce_fraction = []

    energy_alc_multibounce_fraction = []
    number_alc_multibounce_fraction = []

    residual_distribution = zeros(simulation_energy_nbins, simulation_pa_nbins, 0)

    for datafile_idx = eachindex(days_with_data_available)
        print_progress_bar(datafile_idx/length(days_with_data_available))
        if days_with_data_available[datafile_idx] < Date("2020-09-01"); continue; end # for now. avoid bad data

        fullday_event = create_event(days_with_data_available[datafile_idx], satellites[datafile_idx])
        if fullday_event == nothing; continue; end
        if fullday_event.data_reliable == false; continue; end

        event_idxs = []
        for obs_number in 1:fullday_event.n_observations
            append!(event_idxs, _detect_events(fullday_event, obs_number, threshold = 2, show_plot = false))
        end

        # Simulate each event and evaluate
        for detection_idxs in event_idxs
            # Create event
            if detection_idxs == (); continue; end
            event = create_event(fullday_event.time_datetime[detection_idxs[1]], fullday_event.time_datetime[detection_idxs[2]], fullday_event.satellite)
            if event == nothing; continue; end
            if event.data_reliable == false; continue; end

            event_stats = Dict("elfin_range_e_residual" => NaN,
                               "full_range_e_residual" => NaN,
                               "elfin_range_n_residual" => NaN,
                               "full_range_n_residual" => NaN,
                               "residual_distribution" => fill(NaN, simulation_energy_nbins, simulation_pa_nbins),
                               "e_lc_multibounce_fraction" => NaN,
                               "n_lc_multibounce_fraction" => NaN,
                               "e_alc_multibounce_fraction" => NaN,
                               "n_alc_multibounce_fraction" => NaN
            )
            _evaluate_event_simulation(event, event_stats) # Mutates input dict
            if isnan(event_stats["elfin_range_e_residual"]) # If one key is NaN then all keys are NaN
                continue
            end

            push!(elfin_range_energy_residuals, event_stats["elfin_range_e_residual"])
            push!(full_range_energy_residuals, event_stats["full_range_e_residual"])
            push!(elfin_range_number_residuals, event_stats["elfin_range_n_residual"])
            push!(full_range_number_residuals, event_stats["full_range_n_residual"])
            push!(energy_lc_multibounce_fraction, event_stats["e_lc_multibounce_fraction"])
            push!(number_lc_multibounce_fraction, event_stats["n_lc_multibounce_fraction"])
            push!(energy_alc_multibounce_fraction, event_stats["e_alc_multibounce_fraction"])
            push!(number_alc_multibounce_fraction, event_stats["n_alc_multibounce_fraction"])
            residual_distribution = cat(residual_distribution, event_stats["residual_distribution"], dims = 3)
        end
        if datafile_idx % 20 == 0; GC.gc(); end
    end
    # Convert residual distributions into means for storage space
    mean_residual(x) = 10^mean(log10.(x[isfinite.(x)]))
    residuals = [mean_residual(residual_distribution[e,pa,:]) for e in 1:simulation_energy_nbins, pa in 1:simulation_pa_nbins]
    println()

    # save results
    results = Dict("elfin_range_energy_residuals" => float.(elfin_range_energy_residuals),
                   "full_range_energy_residuals" => float.(full_range_energy_residuals),
                   "elfin_range_number_residuals" => float.(elfin_range_number_residuals),
                   "full_range_number_residuals" => float.(full_range_number_residuals),
                   "energy_bin_edges" => simulation_energy_bin_edges,
                   "pitch_angle_bin_edges" => simulation_pa_bin_edges,
                   "residual_distribution" => float.(residuals),
                   "energy_lc_multibounce_fraction" => float.(energy_lc_multibounce_fraction),
                   "number_lc_multibounce_fraction" => float.(number_lc_multibounce_fraction),
                   "energy_alc_multibounce_fraction" => float.(energy_alc_multibounce_fraction),
                   "number_alc_multibounce_fraction" => float.(number_alc_multibounce_fraction),
    )
    npzwrite("$TOP_LEVEL/data/elfin_lifetime_simulation_statistics.npz", results)
    println("Results written to $TOP_LEVEL/data/elfin_lifetime_simulation_statistics.npz")
end

function _evaluate_event_simulation(event::Event, stats)
    # Guard block
    if event == nothing
        return
    end
    # Don't simulate events with limited pitch angle coverage
    if abs(event.avg_pitch_angles[end] - event.avg_pitch_angles[begin]) < 140
        return
    end
    # Don't simulate events that have a hemisphere flip
    if any(abs.(diff(event.loss_cone_angles)) .> 10);
        return
    end

    # Simulate event
    simulation_energy_nbins, simulation_energy_bin_edges, simulation_energy_bin_means, simulation_pa_nbins, simulation_pa_bin_edges, simulation_pa_bin_means, SIMULATION_α_MAX = get_simulation_bins()
    fitted_event = convert_elfin_grid_to_simulation_grid(event, cull_out_of_range = false)
    input = copy(fitted_event)
    input[:, simulation_pa_bin_edges[begin:end-1] .≥ SIMULATION_α_MAX] .= 0
    n_bounces = 10
    output = multibounce_simulation(input, n_bounces, show_progress = false)

    # Get total residuals for ALC
    residual_dict = _get_alc_residuals(event, fitted_event, output)

    # Get flux residual distribution
    sim_flux = dropdims(sum(output[1:2, :, :], dims = 1), dims = 1) # Only use the first 2 distros, as data only should only be bounced once for comparisons
    residual_distribution = sim_flux ./ fitted_event

    # Get multibounce fraction
    α_lc = event.avg_loss_cone_angle
    if α_lc > 90 # Southern hemisphere
        α_lc = 180 - α_lc
    end
    e_lc_multibounce_fraction, n_lc_multibounce_fraction, e_alc_multibounce_fraction, n_alc_multibounce_fraction = get_multibounce_fraction(output, α_lc = α_lc)

    # Mutate input dict to give values back to caller
    stats["elfin_range_e_residual"] = residual_dict["elfin_range_e_residual"]
    stats["full_range_e_residual"] = residual_dict["full_range_e_residual"]
    stats["elfin_range_n_residual"] = residual_dict["elfin_range_n_residual"]
    stats["full_range_n_residual"] = residual_dict["full_range_n_residual"]
    stats["residual_distribution"] = residual_distribution
    stats["e_lc_multibounce_fraction"] = e_lc_multibounce_fraction
    stats["n_lc_multibounce_fraction"] = n_lc_multibounce_fraction
    stats["e_alc_multibounce_fraction"] = e_alc_multibounce_fraction
    stats["n_alc_multibounce_fraction"] = n_alc_multibounce_fraction
    return
end

function _get_alc_residuals(event::Event, fitted_event, distributions)
    @assert length(size(distributions)) == 3 "Improper format"
    simulation_energy_nbins, simulation_energy_bin_edges, simulation_energy_bin_means, simulation_pa_nbins, simulation_pa_bin_edges, simulation_pa_bin_means, SIMULATION_α_MAX = get_simulation_bins()

    # Northern hemisphere normalization
    α_alc = event.avg_anti_loss_cone_angle
    if α_alc < 90 # Southern hemisphere
        α_alc = 180 - α_alc
    end

    # Calculate residuals for ELFIN's recordable range as well as all energies/pitch angles
    results = Dict()
    for integration_mode in ["elfin_range", "full_range"]
        # Get integration bounds
        if integration_mode == "elfin_range" # Only integrating where ELFIN can see
            e_means_idxs = event.energy_bins_min[begin] .≤ simulation_energy_bin_edges[begin:end-1] .≤ event.energy_bins_max[end]
            e_edges_idxs = event.energy_bins_min[begin] .≤ simulation_energy_bin_edges .≤ event.energy_bins_max[end]
   
            pa_means_idxs = α_alc .≤ simulation_pa_bin_edges[begin:end-1] .≤ event.avg_pitch_angles[end]
            pa_edges_idxs = α_alc .≤ simulation_pa_bin_edges .≤ event.avg_pitch_angles[end]
        elseif integration_mode == "full_range"
            e_means_idxs = eachindex(simulation_energy_bin_means)
            e_edges_idxs = eachindex(simulation_energy_bin_edges)

            pa_means_idxs = α_alc .≤ simulation_pa_bin_edges[begin:end-1]
            pa_edges_idxs = α_alc .≤ simulation_pa_bin_edges
        end

            
        e_edges = simulation_energy_bin_edges[e_edges_idxs]
        e_means = edges_to_means(e_edges)
        pa_edges = simulation_pa_bin_edges[pa_edges_idxs]

        ΔE = [(e_edges[e+1] - e_edges[e])/1000 for e in 1:length(e_edges)-1] # MeV
        ΔΩ = [2π * (cosd(pa_edges[α]) - cosd(pa_edges[α+1])) for α in 1:length(pa_edges)-1] # str


        # Get total energy and particle distributions for simulation
        sim_flux = dropdims(sum(distributions[1:2,:,:], dims = 1), dims = 1) # Only use the first 2 distros, as data only should only be bounced once for comparisons
        sim_flux = sim_flux[e_means_idxs, pa_means_idxs]

        sim_particles = [sim_flux[e,α] * ΔE[e] * ΔΩ[α] for e in 1:length(e_edges)-1, α in 1:length(pa_edges)-1]
        sim_energy = [e_means[e] * sim_particles[e,α] for e in 1:length(e_edges)-1, α in 1:length(pa_edges)-1]

        sim_particles = sum(sim_particles)
        sim_energy = sum(sim_energy)

        # Get total energy and particle distributions for data
        data_flux = copy(fitted_event)
        data_flux = data_flux[e_means_idxs, pa_means_idxs]

        data_particles = [data_flux[e,α] * ΔE[e] * ΔΩ[α] for e in 1:length(e_edges)-1, α in 1:length(pa_edges)-1]
        data_energy = [e_means[e] * data_particles[e,α] for e in 1:length(e_edges)-1, α in 1:length(pa_edges)-1]

        data_particles = sum(data_particles)
        data_energy = sum(data_energy)

        # Calculate overall ALC residuals
        energy_residual = sim_energy / data_energy
        particle_residual = sim_particles / data_particles

        results["$(integration_mode)_e_residual"] = energy_residual
        results["$(integration_mode)_n_residual"] = particle_residual
    end
    return results
end

function _detect_events(event::Event, observation_period; show_plot = false, threshold = 2)
    time_slice = event.observation_start_idxs[observation_period]:event.observation_stop_idxs[observation_period]
    if length(time_slice) < 2; return [()]; end

    score = _precipitation_score(event, time_slice)

    event_mask = score .> threshold
    event_mask[end] = 0 # If the last datapoint is high, call it the end of an event. If it's low, no impact.
    event_edges = diff(event_mask)
    if event_edges[1] == -1
        event_edges[1] = 0 # If first datapoint is high and second low, discard detection
    else
        event_edges[1] = event_mask[1] # If the first datapoint is high, call it the start of an event
    end

    if show_plot == true
        _plot_event_detections(event, time_slice, score, event_mask)
    end

    # Return detections
    event_start_idxs = time_slice[findall(event_edges .== 1)]
    event_stop_idxs = time_slice[findall(event_edges .== -1)]
    return collect(zip(event_start_idxs, event_stop_idxs)) # Return a vector of tuples with (start_idx, stop_idx) for each detected event
end

function _plot_event_detections(event::Event, time_slice, score, event_mask)
    # Heatmap
    yticks = Int.(round.(event.energy_bins_mean))
    yticks = yticks[1:3:16]
    heatmap(time_slice, 1:16, log10.(event.Jprec_over_Jtrap[time_slice,:]'),
        background_color_inside = RGB(.75,.75,.75),
        ylabel = "Energy Channel",
        colorbar = false,
        clims = (-1.25, .25),
        xlims = (time_slice[1], time_slice[end]),
        ylims = (1,16),
        yticks = (1:3:16, yticks),
        aspect_ratio = ((time_slice[end] - time_slice[1])/16) * .15,
        c = :ice
    )
    p1 = plot!()

    # Precipitation score
    plot(
        title = "Precipitation Score",
        xlabel = "Data Index",
        ylabel = "Score (arb. units)",
        xlims = (time_slice[1], time_slice[end]),
        ylims = (0, 4),
        aspect_ratio = ((time_slice[end] - time_slice[1])/4) * .15,
    )
    plot!(time_slice, score, label = false)
    p2 = plot!()

    # Detection area
    plot(time_slice, 1 .* event_mask,
        fill = true,
        linetype = :steppre,
        fillcolor = :red,
        linecolor = :transparent,
        xlims = (time_slice[1], time_slice[end]),
        ylims = (0, 1),
        framestyle = :none,
        grid = false,
        label = false,
        topmargin = 5mm,
        bottommargin = -10mm
    )
    p3 = plot!()

    layout = @layout [a{.1h}; b; c]
    plot(p3, p1, p2,
        layout = layout,
        suptitle = "$(event.name)",
        dpi = 350
    )
    display("image/png",plot!())
end    

function _precipitation_score(event::Event, time_slice; smoothing_radius = 4)
    if _check_for_bad_data(event, time_slice) == true; return zeros(length(time_slice)); end

    Jprec_over_Jtrap = copy(event.Jprec_over_Jtrap[time_slice, :])

    # Kill pixels with low amounts of precipitation
    Jprec_over_Jtrap[Jprec_over_Jtrap .< 10^(-.75)] .= NaN

    # Find pixels with few neighbors and kill them
    old_J_Jprec_over_Jtrap = copy(Jprec_over_Jtrap) # To avoid feedback loop of pixel killing
    for t = 1:length(time_slice)
        t_to_check = clamp.(t-1:t+1, 1, length(time_slice))
        t_to_check = unique(t_to_check) # Don't double count pixels on the edges

        for e = 1:16
            e_to_check = clamp.(e-1:e+1, 1, 16)
            e_to_check = unique(e_to_check) # Don't double count pixels on the edges

            n_neighbors = sum(isfinite.(old_J_Jprec_over_Jtrap[t_to_check, e_to_check])) - 1 # Minus one because this includes current pixel. It can go below zero if current pixel is non-finite but that's fine.
            if n_neighbors < 3
                Jprec_over_Jtrap[t,e] = NaN
            end
        end
    end

    # Kill nonfinite pixels
    nonfinite_idxs = findall(.!isfinite.(Jprec_over_Jtrap))
    Jprec_over_Jtrap[nonfinite_idxs] .= 0

    # Sum vertically for score
    score = dropdims(sum(Jprec_over_Jtrap, dims = 2), dims = 2)

    # Apply smoothing
    smoothed_score = copy(score) .* 0
    for i = 1:length(smoothed_score)
        lower_idx = clamp(i - smoothing_radius, 1, length(time_slice))
        upper_idx = clamp(i + smoothing_radius, 1, length(time_slice))
        smoothed_score[i] = max(score[lower_idx:upper_idx]...)
    end

    # Return
    return smoothed_score
end

function _check_for_bad_data(event::Event, slice)
    # Low energy channels hot
    e_flux, _ = integrate_flux(event, energy = true)
    slice_avg = mean(e_flux[slice,:])
    slice_std = std(e_flux[slice,:])

    if log10(slice_std) <= log10(slice_avg) - .5 # These are events with hot low channels and a very hot and uniform PAD.
        return true
    end
end



# ---------------- Script ----------------
println("========= Script Begin =========")

write_example_event_detection_data() # Save example of event detection for poster figure
calculate_elfin_simulation_stats() # Save analysis statistics

println("========== Script End ==========")