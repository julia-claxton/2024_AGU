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
function calculate_simulation_statistics(; loss_e_nbins = 33, loss_pa_nbins = 36, e_min = .1, n_bounces = 10)
    println("Calculating backscatter and multibounce statistics as a function of energy & pitch angle...")
    
    # Rebin backscatter data with lower minimum energy bound to get accurate loss rates for low energy input beams
    println("Rebinning backscatter data with expanded energy bounds...")
    set_simulation_bins(e_min = e_min) # keV
    simulation_energy_nbins, simulation_energy_bin_edges, simulation_energy_bin_means, simulation_pa_nbins, simulation_pa_bin_edges, simulation_pa_bin_means, SIMULATION_α_MAX = get_simulation_bins()
    
    loss_e_bin_edges = 10 .^ LinRange(log10.(e_min), 4, loss_e_nbins+1)
    loss_e_bin_means = edges_to_means(loss_e_bin_edges)
    loss_pa_bin_edges = LinRange(0, 180, loss_pa_nbins+1)
    loss_pa_bin_means = edges_to_means(loss_pa_bin_edges)

    println("Calculating loss rate, multibounce fraction, and survival rate distributions...")
    # Preallocate results array
    backscatter_fraction = zeros(loss_e_nbins, loss_pa_nbins)
    multibounce_fraction = zeros(loss_e_nbins, loss_pa_nbins)
    survival_fraction = zeros(loss_e_nbins, loss_pa_nbins)

    coordinates = CartesianIndices(backscatter_fraction)
    for i in 1:length(coordinates)
        print_progress_bar(i/length(coordinates))

        e_idx = coordinates[i][1]
        α_idx = coordinates[i][2]

        # Kick out beams outside of simulation range
        if loss_pa_bin_edges[α_idx] > SIMULATION_α_MAX
            continue
        end

        # Create input beam
        input_beam = zeros(simulation_energy_nbins, simulation_pa_nbins)
        
        _, sim_e_idx = findmin(abs.(loss_e_bin_means[e_idx] .- simulation_energy_bin_means))
        sim_e_idx = sim_e_idx[1]
        if sim_e_idx == nothing; continue; end

        sim_pa_idx = findfirst(loss_pa_bin_edges[α_idx] .≤ simulation_pa_bin_means .< loss_pa_bin_edges[α_idx+1])
    
        input_beam[sim_e_idx, sim_pa_idx] = 1e6

        # Get backscatter
        output_distributions = multibounce_simulation(input_beam, n_bounces)

        # Get loss rate
        backscatter_fraction[e_idx, α_idx] = 1 - get_atmosphere_loss_rate(output_distributions)

        # Get multibounce fraction
        _, lc_multibounce_fraction_n, _, _ = get_multibounce_fraction(output_distributions)
        multibounce_fraction[e_idx, α_idx] = lc_multibounce_fraction_n

        # Get survival fraction
        survival_fraction[e_idx, α_idx] = get_steady_state_survival_fraction(output_distributions)
    end
    println()

    print("Saving results to $(TOP_LEVEL)/data/loss_cone_simulation_statistics.npz... ")
    npzwrite("$TOP_LEVEL/data/loss_cone_simulation_statistics.npz",
        Dict("backscatter_fraction" => backscatter_fraction,
             "multibounce_fraction" => multibounce_fraction,
             "survival_fraction" => survival_fraction,
             "e_bin_edges" => loss_e_bin_edges,
             "pa_bin_edges" => loss_pa_bin_edges
        )
    )
     # Reset bins to default before exiting
     println("Rebinning backscatter data to default...")
    set_simulation_bins()

    println("Done\n")
end

# ---------------- Script ----------------
println("========= Script Begin =========")
calculate_simulation_statistics()

println("========== Script End ==========")