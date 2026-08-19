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
