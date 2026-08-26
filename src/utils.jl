"""
Script to generate sobol sequences for parameter combinations

Order of parameters in sobol sequence are as follows:

1. Number of species: 9.5 - 50.5, rounded to nearest integer
2. Consumer-resource body-mass ratio (Z): 10 - 100
3. Connectance (C): 0.1 - 0.3
REMOVED: 4. Hill exponent (h): 1 - 2
REMOVED: 5. Noise decay rate (θ): 0.01 - 1.0
REMOVED: 6. Noise variance (σ): 0.02 - 0.5
"""

function generate_parameters(df::Union{DataFrame, Nothing} = nothing; n::Int = 2000)

    # Generate an empty dataframe if needed
    if isnothing(df)
        df = DataFrame(n_species = [], Z = [], C = [])
    end

    # Define the sequence
    seq = SobolSeq([9.5, 10, 0.1],[50.5, 100, 0.3])

    # Skip the already used values
    Sobol.next!(skip(seq, nrow(df), exact = true))

    # Generate the next n values in the sequence
    params_ = reduce(hcat, Sobol.next!(seq) for i = 1:n)'

    # Convert to dataframe 
    params_df = DataFrame(params_, [:n_species, :Z, :C])
    df = vcat(df, params_df)

    # Round & convert number of species to an integer
    df.:n_species = Int.(round.(df.:n_species))

    return df
end

# Quick test
begin
    df1 = generate_parameters(n = 5)
    df2 = generate_parameters(df1, n = 5)
    df3 = generate_parameters(n = 10)
    @assert df2 == df3
end

"""
Internal function:
Removes unconnected producers from the burn-in outputs during prune_foodweb

 • Find out who is extinct
 • For each producer, count number of consumers ( > 0 in that column of the interaction matrix)
 • Set biomass to 0 for all producers that have no extant consumers
 • Returns vector of final biomasses corrected for unconnected producers
"""
function remove_unconnected_producers(foodweb::FoodWeb, final_biomasses::Vector{Float64})

    # Step 1 - who is extinct?
    extinct_species = findall(x -> x == 0, final_biomasses)

    # For each producer...
    for i in producers(foodweb)

        # Step 2 - who consumes each producer?
        consumers = findall(x -> x > 0, foodweb.A[:,i])

        # Step 3 - which consumers are extant?
        extant_cons = setdiff(consumers, extinct_species)

        # Step 4 - set producer biomass to 0 if all consumers extinct (or never existed)
        if length(extant_cons) == 0
            final_biomasses[i] = 0
        end
    end
    return final_biomasses
end

"""
Internal function:
Removes unconnected producers from the burn-in outputs during prune_foodweb
    
 • An unconnected consumer has no prey
 • A producer has no prey
 • producers returns a vector of "producers" - species with no prey
 • The challenge is to differentiate between these 

Potential methods:

 ◎ By mass? (1.0 by default in producers, greater in consumers if Z > 1) - but not in all cases
 ◎ By functional response parameter (ω, eᵢⱼ, Fᵢⱼ); > 0 in consumers
    ∴ if sum(ω) for "producers" is > 0, some are unconnected consumers
"""
function remove_unconnected_consumers(p::ModelParameters, biomass::Vector{Float64})
    
    # Step 1 - Who has no prey?
    prods = producers(p.network)

    # While some "producers" have a resource preference (i.e. are really consumers)...
    while sum(p.functional_response.ω[prods,:]) > 0
        # Step 2 - What is the index within the "producer" vector of the consumers
        unconnected_consumer_index = findall(x -> x > 0, sum(p.functional_response.ω, dims = 2)[prods])
        
        # Step 3 - Which species are not unconnected
        connected_species = setdiff(1:richness(p.network), prods[unconnected_consumer_index])

        # Step 4 - Retain species interactions, masses, consumer preferences, & biomasses
        cropped_fw = p.network.A[connected_species,connected_species]
        cropped_mass = p.network.M[connected_species]
        cropped_preference = p.functional_response.ω[connected_species, connected_species]
        biomass = biomass[connected_species]
        
        # Step 5 - Rebuild the ModelParameters
        fw = FoodWeb(cropped_fw, M = cropped_mass)
        br = BioRates(fw, d = 0)
        fr = BioenergeticResponse(fw, ω = cropped_preference, h = p.h)

        if p.stochasticity.addstochasticity
            as = AddStochasticity(fw, addstochasticity = true, target = "producers", θ = p.stochasticity.θ[1], σe = p.stochasticity.σe[1])
        else
            as = AddStochasticity(fw, addstochasticity = false)
        end

        if p.stressor.addstressor
            s = Stressor(fw, addstressor = true, slope = -0.0001, start = 1000.0)
        else
            s = Stressor(fw, addstressor = false)
        end

        p = ModelParameters(fw,
            biorates = br,
            functional_response = fr,
            stochasticity = as,
            stressor = s)
        
        prods = producers(p.network)
    end
    p, biomass
end

"""
Prunes the food webs in burn-in outputs by removing species that are extinct or unconnected.
With my current use of defaults a species is defined by:
 ◎ Interactions
 ◎ Body mass
 ◎ Resource preference (ω)

Even if all species are extinct, needs to return a vector of same length as goes in...
Otherwise would mismatch with params_set 
"""
function prune_foodweb(params::ModelParameters, final_biomasses::Vector{Float64}, i)
    
    # Only interested in species biomasses
    final_biomasses = final_biomasses[1:richness(params.network)]
    
    # Set biomasses of unconnected producers to 0
    extinct_or_unconnected = remove_unconnected_producers(params.network, final_biomasses)

    # Now extant species is only connected species with biomass
    extant_species = findall(x -> x != 0, extinct_or_unconnected)

    # Retained Characteristics
    retained_interactions = params.network.A[extant_species, extant_species]
    retained_bodymass = params.network.M[extant_species]
    retained_preference = params.functional_response.ω[extant_species,extant_species]

    retained_decay = params.stochasticity.θ[1]
    retained_variance = params.stochasticity.σe[1]

    # If some species remain, rebuild ModelParameters
    if size(retained_interactions, 1) != 0
        pruned_web = FoodWeb(retained_interactions, M = retained_bodymass)
        pruned_fr = BioenergeticResponse(pruned_web, ω = retained_preference, h = params.functional_response.h)
        pruned_rates = BioRates(pruned_web, d = 0)

        if params.stochasticity.addstochasticity
            pruned_stoch = AddStochasticity(pruned_web, addstochasticity = true, target = "producers", θ = retained_decay, σe = retained_variance)
        else
            pruned_stoch = AddStochasticity(pruned_web, addstochasticity = false)
        end

        if params.stressor.addstressor
            pruned_stress = Stressor(pruned_web, addstressor = true, slope = -0.0001, start = 1000.0)
        else
            pruned_stress = Stressor(pruned_web, addstressor = false)
        end

        pruned_params = ModelParameters(pruned_web,
            biorates = pruned_rates,
            functional_response = pruned_fr,
            stochasticity = pruned_stoch,
            stressor = pruned_stress)

        pruned_biomasses = final_biomasses[extant_species]

        if iszero(length(producers(pruned_web)))
            pruned_params = "Nothing left!"
        elseif unique(pruned_params.network.M[producers(pruned_params.network)]) != [1.0]
            # Can remove_unconnected_consumers here because it won't unconnect any producers
            pruned_params, pruned_biomasses = remove_unconnected_consumers(pruned_params, pruned_biomasses)
        end
        if typeof(pruned_params) == ModelParameters && richness(pruned_params.network) > 1
            return (rep = i, params = pruned_params, pruned_biomasses = pruned_biomasses)
        end # Here remove_unconnected_consumers has removed all remaining species
    end # No species remain
end

"""
Calculates the Dominant Eigenvalue of the Variance-Covariance Matrix (D.E.V.C.M)

A variance-covariance matrix for the timeseries of a 3 species food web is constructed as follows:

S1 = Species 1's timeseries
S2 = Species 2's timeseries
S3 = Species 3's timeseries

[ var(S1)      cov(S1, S2)  cov(S1, S3) 
  cov(S1, S2)  var(S2)      cov(S2, S3) 
  cov(S1, S3)  cov(S2, S3)  var(S3)     ]

Note: cov(S1, S2) == cov(S2, S1)
"""
function DEVCM(timeseries::Matrix{Float64})

    S = size(timeseries, 2)

    mat = zeros(S, S)

    for i in axes(timeseries, 2)
        for j in axes(timeseries, 2)
            if i == j
                mat[i, j] = var(timeseries[:,i])
            else
                mat[i, j] = cov(timeseries[:,i], timeseries[:,j])
            end
        end
    end

    maximum(eigvals(mat))
end
