println("Rundir is $(pwd())")

using DataFrames, Arrow

using stoch_befw, Distributions, Sobol

include("../src/utils.jl")

params_vec = []

for θ in [0.02, 0.1, 0.5]
    for cov in [0.05, 0.1, 0.2]

        restored_df = DataFrame(Arrow.Table("hill_2_deterministic_burnin_combined.arrow"))

          biomasses = Vector{Vector{Float64}}(restored_df.biomasses)

          reps = Vector{Int64}(restored_df.rep)
          interaction_matrices = Matrix.(reshape.(restored_df.A, length.(biomasses), length.(biomasses)))
          body_masses = Vector{Float64}.(restored_df.M)
          h = Float64.(restored_df.h)
          ω_matrices = Matrix{Float64}.(reshape.(restored_df.ω, length.(biomasses), length.(biomasses)))

          model_params = []

          for i in eachindex(restored_df.A)
          
            fw = FoodWeb(interaction_matrices[i], M = body_masses[i])
            br = BioRates(fw, d = 0)
            fr = BioenergeticResponse(fw, h = h[i], ω = ω_matrices[i])
            as = AddStochasticity(fw, addstochasticity = true, target = "producers", θ = θ, σe = sqrt(2 * θ * cov))
            MP = ModelParameters(fw, biorates = br, functional_response = fr, stochasticity = as)
            push!(model_params, MP)
          
          end
          
          pruned_foodwebs = []
          for i in eachindex(model_params)
              pruned = prune_foodweb(model_params[i], biomasses[i], reps[i])
              if !isnothing(pruned) # Here we only want food webs 
                push!(pruned_foodwebs, pruned)
              end
              if rem(i, 1000) == 0
                println(i)
              end
          end

          new_pruned_foodwebs = []
          for i in eachindex(pruned_foodwebs)
              fw = pruned_foodwebs[i].params.network
              br = BioRates(fw, d = 0)
              fr = BioenergeticResponse(fw, h = h[i], ω = pruned_foodwebs[i].params.functional_response.ω)
              as = AddStochasticity(fw, addstochasticity = true, target = "producers", θ = θ, σe = sqrt(2 * θ * cov))
              MP = ModelParameters(fw, biorates = br, functional_response = fr, stochasticity = as)
              rep = pruned_foodwebs[i].rep
              biomasses = pruned_foodwebs[i].pruned_biomasses
              push!(new_pruned_foodwebs, (rep = rep, params = MP, pruned_biomasses = biomasses))
          end

          push!(params_vec, new_pruned_foodwebs)

    end
end

params_vec = reduce(vcat, params_vec)

rep = []
A = []
M = []
h = []
ω = []
σ = []
θ = []
biomasses = []
last_extinction_point = []

for i in eachindex(params_vec)
  push!(rep, params_vec[i].rep)
  push!(A, params_vec[i].params.network.A)
  push!(M, params_vec[i].params.network.M)
  push!(h, params_vec[i].params.functional_response.h)
  push!(ω, params_vec[i].params.functional_response.ω)
  push!(σ, params_vec[i].params.stochasticity.σe)
  push!(θ, params_vec[i].params.stochasticity.θ)
  push!(biomasses, params_vec[i].pruned_biomasses)

end

df = DataFrame(rep = rep, A = A, M = M, h = h, ω = ω, σ = σ, θ = θ, biomasses = biomasses)

Arrow.write("hill_2_params_set.arrow", df)

