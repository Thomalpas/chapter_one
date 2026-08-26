
# For calculating all extinction events

# Open timeseries

using DataFrames, Sobol
using stoch_befw
using Arrow, LinearAlgebra
using CSV, StatsBase

include("../src/utils.jl")

script_location = try
    ARGS[3]
catch
    "local"
end

function find_extinction_points(mat::Matrix{Float64})
    extinction_points = []
    for i in 1:size(mat,2) # each column
        nzeros = length(findall(x -> x == 0, mat[:,i])) # number of zeros in the column
        if nzeros > 0 # if there are any zeros
            extinction_point = size(mat, 1) - nzeros
            push!(extinction_points, extinction_point)
        end
    end
    return extinction_points
end

function matrix_nonzero_length(mat::Matrix{Float64})
    nzeros = Int64.(zeros(size(mat,2)))
    for i in 1:size(mat,2) # each column
        nzeros[i] = length(findall(x -> x != 0, mat[:,i])) # number of zeros in the column
    end
    nzeros
end

println("packages loaded!")

# Set up

first_sim, last_sim = try parse(Int, ARGS[1]), parse(Int, ARGS[2])
catch
    1, 20
end

# timeseries vec contains timeseries from 1 to either 
# a) termination (ususally if everything dead)
# b) end of simulation. 16400 time steps

hill = 2
step_length = 150

for win_size in [100, 200, 400]
    for sample_resolution in [1, 2, 4]

        all_foodwebs = []
        original_species_id = []
        realised_species_id = []
        foodweb_id = []
        extinction_id = []

        overall_df = DataFrame()
        no_extinctions = []

        if script_location == "HPC"
            file = string("../../../../mnt/parscratch/users/bi1tma/hill_$(hill)_stress_timeseries_", first_sim, "_", last_sim, ".arrow")
        elseif script_location == "local"
            file = string("./out/hill_$(hill)_stress_timeseries_", first_sim, "_", last_sim, ".arrow")
        end
        
        all_stress_timeseries = DataFrame(Arrow.Table(file))

        all_stress_timeseries.M = Vector{Vector{Float64}}(all_stress_timeseries.M)

        # the timeseries may not all be the same length - could terminate early!

        timeseries_vec = []
        for i in [1:1:size(all_stress_timeseries)[1];]

            push!(timeseries_vec, reshape(all_stress_timeseries.biomasses[i], 
                                    Int64(length(all_stress_timeseries.biomasses[i])/length(all_stress_timeseries.M[i])),
                                    length(all_stress_timeseries.M[i])))
        end
        timeseries_vec

        food_webs = []
        for i in [1:1:size(all_stress_timeseries)[1];]
        push!(food_webs, FoodWeb(reshape(all_stress_timeseries.A[i], length(all_stress_timeseries.M[i]), length(all_stress_timeseries.M[i])), M = all_stress_timeseries.M[i]))
        end

        for timeseries in eachindex(timeseries_vec)

            extinction_timepoints = find_extinction_points(Matrix{Float64}(timeseries_vec[timeseries]))

            if length(extinction_timepoints) > 0 && minimum(extinction_timepoints) > 6401

                # of the original timeseries, the first window ends 1900 time points in.
                # n_steps is how many more windows can we fit in before the final extinction point?

                n_steps = Int64(floor((maximum(extinction_timepoints) - 1900)/step_length))

                if n_steps >= 0

                    # First 6401 points are pre-stress - so first stress window ends at 6400

                    # we want 30 windows before? so 30 * 150 = 4500; 1900
                    # but with a maximum window size of 400; we need 4900 pre-stress timesteps [1500:6400] + all stress timesteps until last extinction

                    mat = Matrix{Float64}(timeseries_vec[timeseries][1500:maximum(extinction_timepoints),:])


                    # First window is 400 long, sampled every 2 timesteps (sampling_resolution) - then lag-1 autocorrelation
                    # Leaves a matrix 100 points x species richness

                    # We can then move this window by 150 steps (step_length) and calculate again - until we reach the last consumer extinction

                    biomass_change_vec = []
                    mean_biomass_vec = []
                    rel_biomass_change_vec = []

                    autocor_vec = []
                    SD_vec = []
                    SK_vec = []
                    KURT_vec = []
                    DEVCM_vec = []

                    for i in 0:n_steps - 1 #(to get from the 1st to the nth window needs n-1 steps)
                        
                        window = mat[1 + (i*step_length): win_size + (i*step_length), :][1:sample_resolution:end,:]

                        # Let's put in change in biomass here: log10.(window)
                        push!(biomass_change_vec, (log10.(window[end,:]) .- log10.(window[1,:]))')
                        push!(mean_biomass_vec, mapslices(mean, log10.(window); dims = 1))
                        push!(rel_biomass_change_vec, ((log10.(window[end,:]) .- log10.(window[1,:]))') ./ mapslices(mean, log10.(window); dims = 1)    )

                        # Autocorrelation
                        push!(autocor_vec, autocor(window, [1]))
                        # SD of log biomass for a block
                        push!(SD_vec, mapslices(std, log10.(window); dims=1))
                        # Skewness of a block
                        push!(SK_vec, mapslices(skewness, log10.(window); dims=1))
                        # Kurtosis of a block
                        push!(KURT_vec, mapslices(kurtosis, log10.(window); dims=1))
                        # Dominant Eigenvalue of the Variance Covariance Matrix
                        push!(DEVCM_vec, repeat([DEVCM(window)], size(window, 2)))
                    end

                    biomass_change_vec = reduce(vcat, biomass_change_vec)
                    mean_biomass_vec = reduce(vcat, mean_biomass_vec)
                    rel_biomass_change_vec = reduce(vcat, rel_biomass_change_vec)
                    autocor_vec = reduce(vcat, autocor_vec)
                    SD_vec = reduce(vcat, SD_vec)
                    SK_vec = reduce(vcat, SK_vec)
                    KURT_vec = reduce(vcat, KURT_vec)
                    DEVCM_vec = reduce(vcat, DEVCM_vec)
                
                    # time to put it in a DataFrame
                    foodweb = repeat([all_stress_timeseries.rep[timeseries]], inner = richness(food_webs[timeseries]) * (n_steps))
                    species_id = repeat(1:richness(food_webs[timeseries]), inner = n_steps)
                    timepoint_vec = repeat([1:1:n_steps;], outer = richness(food_webs[timeseries]))

                    tidy_biomass_vec = reduce(vcat, biomass_change_vec)
                    tidy_mean_biomass_vec = reduce(vcat, mean_biomass_vec)
                    tidy_rel_biomass_change_vec = reduce(vcat, rel_biomass_change_vec)
                    tidy_autocor_vec = reduce(vcat, autocor_vec)
                    tidy_SD_vec = reduce(vcat, SD_vec)
                    tidy_SK_vec = reduce(vcat, SK_vec)
                    tidy_KURT_vec = reduce(vcat, KURT_vec)
                    tidy_DEVCM_vec = reduce(vcat, DEVCM_vec)
                    
                    tidy_biomass_vec = round.(tidy_biomass_vec, sigdigits = 5)
                    tidy_mean_biomass_vec = round.(tidy_mean_biomass_vec, sigdigits = 5)
                    tidy_rel_biomass_change_vec = round.(tidy_rel_biomass_change_vec, sigdigits = 5)
                    tidy_autocor_vec = round.(tidy_autocor_vec, sigdigits = 5)
                    tidy_SD_vec = round.(tidy_SD_vec, sigdigits = 5)
                    tidy_SK_vec = round.(tidy_SK_vec, sigdigits = 5)
                    tidy_KURT_vec = round.(tidy_KURT_vec, sigdigits = 5)
                    tidy_DEVCM_vec = round.(tidy_DEVCM_vec, sigdigits = 5)


                    # Sort out when extinctions occur

                    # will all be values between 0 & n_steps
                    all_extinction_points = Int64.(floor.((extinction_timepoints .- 1900)./step_length))

                    # starts with just 0s - will fill in 1s where extinction occurs.
                    # essentially. each species in order, n_steps 0s
                    extinction = zeros(richness(food_webs[timeseries]) * (n_steps))

                    # Need to calculate indices of extinct species

                    # matrix_nonzero_length finds how many nonzero entries there are in each column of the matrix timeseries. 
                    # minusing the timeseries length means if nonzeros == timeseries length (doesn't go extinct), then doesn't flag here
                    # which columns contain zeros - asking whether the last entry is a zero might do the same job?

                    extinct_species_idx = findall(x -> x != 0, matrix_nonzero_length(Matrix{Float64}(timeseries_vec[timeseries])) .- size(timeseries_vec[timeseries], 1))

                    points_of_extinction = ((extinct_species_idx .- 1) * (n_steps)) .+ all_extinction_points

                    # Put the extinctions in the right place
                    extinction[Int64.(points_of_extinction)] .= 1.0

                    # Now to remove the timesteps after the extinction has occurred

                    post_extinct_timesteps = []
                    for i in eachindex(extinct_species_idx)
                        push!(post_extinct_timesteps, [points_of_extinction[i] + 1 : 1 : extinct_species_idx[i] * (n_steps);])
                    end
                    post_extinct_timesteps = reduce(vcat, post_extinct_timesteps)

                    occurred_extinctions = zeros(n_steps)
                    for i in eachindex(occurred_extinctions)
                        occurred_extinctions[i] = length(findall(x -> x < i, unique(sort(all_extinction_points)))) .+ 1
                    end

                    approaching_extinction = Int64.(repeat(occurred_extinctions, outer = richness(food_webs[timeseries])))

                    # Time to calculate all foodwebs as extinctions occur

                    endpoint_of_windows_with_extinction = sort(unique(((all_extinction_points) .* step_length) .+ 1900))

                    for i in eachindex(endpoint_of_windows_with_extinction)

                        extant_species = setdiff([1:1:richness(food_webs[timeseries]);], findall(x -> x == 0, timeseries_vec[timeseries][endpoint_of_windows_with_extinction[i], :]))
                        
                        push!(original_species_id, extant_species)
                        push!(realised_species_id, [1:1:length(extant_species);])
                        push!(foodweb_id, repeat([all_stress_timeseries.rep[timeseries]], inner = length(extant_species)))
                        push!(extinction_id, repeat([i], inner = length(extant_species)))

                        push!(all_foodwebs, FoodWeb(food_webs[timeseries].A[extant_species, extant_species]))

                    end

                    df = DataFrame()

                    df.network = foodweb
                    df.species = species_id
                    df.time = timepoint_vec
                    df.approaching_extinction = approaching_extinction
                    df.extinction = extinction
                    df.biomass_change = tidy_biomass_vec
                    df.mean_biomass_vec = tidy_mean_biomass_vec
                    df.rel_biomass_change_vec = tidy_rel_biomass_change_vec

                    df.autocor = tidy_autocor_vec
                    df.SD = tidy_SD_vec
                    df.SK = tidy_SK_vec
                    df.KURT = tidy_KURT_vec
                    df.DEVCM = tidy_DEVCM_vec

                    deleteat!(df, post_extinct_timesteps)

                    append!(overall_df, df)
                else
                    push!(no_extinctions, timeseries)
                end

            else
                push!(no_extinctions, timeseries)
            end

            if rem(timeseries, 100) == 0
                println(timeseries)
            end

        end


        # This is the overall object we produce:

        overall_df

        # Just confirm that every time we have an extinction event we have a new foodweb
        all_foodwebs

        @assert nrow(unique(overall_df[!, [:network, :approaching_extinction]])) == length(all_foodwebs)

        # When we calculate network properties they will be for the realised species index i.e. within the food web after extinction has occurred
        # We need to be able to match up these new indices to the same original species
        # i.e. these are the network properties associated with species i from food web j after k extinctions...

        foodweb_id = reduce(vcat, foodweb_id)
        extinction_id = reduce(vcat, extinction_id)
        realised_species_id = reduce(vcat, realised_species_id)
        original_species_id = reduce(vcat, original_species_id)

        species_index_conversion = DataFrame()

        species_index_conversion.foodweb_id = foodweb_id
        species_index_conversion.extinction_id = extinction_id
        species_index_conversion.realised_species_id = realised_species_id
        species_index_conversion.original_species_id = original_species_id

        # Let's take a look

        species_index_conversion

        # We need something for the realised foodweb id

        # There is a realised food web for each combination of food web and extinction
        species_index_conversion.realised_foodweb = string.(species_index_conversion.foodweb_id).*"_".*string.(species_index_conversion.extinction_id)

        temp_df = DataFrame()
        temp_df.realised_foodweb = unique(string.(species_index_conversion.foodweb_id).*"_".*string.(species_index_conversion.extinction_id))
        temp_df.realised_foodweb_id = [1:1:length(unique(string.(species_index_conversion.foodweb_id).*"_".*string.(species_index_conversion.extinction_id)));]
        temp_df

        combined_df = innerjoin(species_index_conversion, temp_df, on = :realised_foodweb)

        # Drop the temporary column string I used to join

        select!(combined_df, Not(:realised_foodweb))

        all_foodwebs

        CSV.write("./out/timeseries_properties/hill_$(hill)_winsize_$(win_size)_resolution_$(sample_resolution)_timeseries_properties_$(first_sim)_$(last_sim).csv", overall_df)
        CSV.write("./out/timeseries_properties/hill_$(hill)_winsize_$(win_size)_resolution_$(sample_resolution)_index_conversion_$(first_sim)_$(last_sim).csv", combined_df)

        # Can't save food webs like FoodWeb anymore...

        realised_foodweb_id = []
        A = []
        M = []

        for i in eachindex(all_foodwebs)
            push!(realised_foodweb_id, i)
            push!(A, all_foodwebs[i].A)
            push!(M, all_foodwebs[i].M)

        end

        df = DataFrame(realised_foodweb_id = realised_foodweb_id, A = A, M = M)

        Arrow.write("./out/timeseries_properties/hill_$(hill)_winsize_$(win_size)_resolution_$(sample_resolution)_foodwebs_$(first_sim)_$(last_sim).arrow", df)
    end

end
