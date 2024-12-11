const TOP_LEVEL = dirname(@__DIR__)
@assert endswith(TOP_LEVEL, "2024_AGU_Multibounce") "Script in wrong directory! TOP_LEVEL = $TOP_LEVEL"

# Library includes
using Statistics
using LinearAlgebra
using NPZ

# Import my ELFIN tools
include("./Julia_ELFIN_Tools/Events.jl")
include("./Julia_ELFIN_Tools/Visualization.jl")

include("./EPPBackscatterSimulation/BackscatterSimulation.jl") # Provides backscatter simulation capability
include("./ELFIN_Backscatter_Tools.jl") # Provides interface between ELFIN data and backscatter simulation
include("./General_Functions.jl") # Provides general-purpose functions I find useful

# ---------------- Functions ----------------
function write_fitting_data()
    event = create_event(DateTime("2021-03-06T07:02:00"), DateTime("2021-03-06T07:08:00"), "b")
    data_nflux = convert_elfin_grid_to_simulation_grid(event, cull_out_of_range = false)

    _, e_edges, e_means, _, pa_edges, pa_means, SIMULATION_α_MAX = get_simulation_bins()
    alc_idxs = pa_means .≥ event.avg_anti_loss_cone_angle

    sim_backscatter, beam_coordinates, beam_strengths, fit_energy_fraction, fit_number_fraction = simulate_NH_backscatter(e_edges, pa_edges, data_nflux, return_fit_info = true)
    beam_coordinates = [[beam_coordinates[i][2], beam_coordinates[i][1]] for i in eachindex(beam_coordinates)]
    beam_coordinates = cat(beam_coordinates..., dims = 2)'
    beam_colors = replace(log10.(beam_strengths), -Inf => -100) # -Inf breaks zcolor arg in scatter!()

    npzwrite("$(TOP_LEVEL)/data/beam_fitting.npz",
        Dict("e_edges" => e_edges,
             "pa_edges" => pa_edges,
             "data_nflux" => data_nflux,
             "beam_coordinates" => beam_coordinates,
             "beam_colors" => beam_colors,
             "fit_energy_fraction" => fit_energy_fraction,
             "fit_number_fraction" => fit_number_fraction
        )
    )
end


# ---------------- Script ----------------
println("========= Script Begin =========")
write_fitting_data()
println("========== Script End ==========")