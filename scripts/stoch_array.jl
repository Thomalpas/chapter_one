println("Rundir is $(pwd())")

using Distributed, Serialization

first_sim, last_sim = try parse(Int, ARGS[1]), parse(Int, ARGS[2])
catch
    1, 90
end

try 
  println("Running parameters combination from $(ARGS[1]) to $(ARGS[2])")
catch
  println("Running parameters combination from $first_sim to $last_sim")
end
  
println(length(Sys.cpu_info()))

#ncpu = maximum([length(Sys.cpu_info()), 15])
ncpu = 10

#Flag enables all the workers to start on the project of the current dir
#flag = "--project=~/julia_coding/chapter_one/" # HPC flag
flag = "--project=~/Documents/julia_coding/chapter_one/" # local machine flag
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(ncpu - 1, exeflags=flag)
println("Using $(ncpu -1) cores")

@everywhere import Pkg, Random.seed!
@everywhere Pkg.instantiate()
@everywhere using DataFrames, Sobol
@everywhere using stoch_befw
@everywhere using StatsBase, Distributions
@everywhere using ProgressMeter

@everywhere include("../src/stoch_sim.jl")

@everywhere include("../src/utils.jl")

# Step one, bring in data

# params_set = try DataFrame(Arrow.Table("../h_1_params_set.arrow"))
# catch
#   DataFrame(Arrow.Table("h_1_params_set.arrow"))
# end

using Arrow

params_set = try DataFrame(Arrow.Table("../hill_2_params_set.arrow"))
catch
  DataFrame(Arrow.Table("hill_2_params_set.arrow"))
end

biomasses = Vector{Vector{Float64}}(params_set.biomasses)

reps = Vector{Int64}(params_set.rep)
interaction_matrices = Matrix.(reshape.(params_set.A, length.(biomasses), length.(biomasses)))
body_masses = Vector{Float64}.(params_set.M)
h = Float64.(params_set.h)
ω_matrices = Matrix{Float64}.(reshape.(params_set.ω, length.(biomasses), length.(biomasses)))
θ = first.(params_set.θ)
σ = first.(params_set.σ)

model_params = []

for i in eachindex(params_set.A)

  fw = FoodWeb(interaction_matrices[i], M = body_masses[i])
  br = BioRates(fw, d = 0)
  fr = BioenergeticResponse(fw, h = h[i], ω = ω_matrices[i])
  as = AddStochasticity(fw, addstochasticity = true, target = "producers", θ = θ[i], σe = σ[i])
  MP = ModelParameters(fw, biorates = br, functional_response = fr, stochasticity = as)
  push!(model_params, MP)

end

temp_df = (rep = reps, params = model_params, pruned_biomasses = biomasses)

pruned_foodwebs = []
for i in eachindex(temp_df.rep)
    push!(pruned_foodwebs, (rep = temp_df.rep[i], params = temp_df.params[i], pruned_biomasses = temp_df.pruned_biomasses[i]))
end

pm = sample(pruned_foodwebs)
println("Running warmup")

warmup = mega_stoch_sim(pm.params, pm.pruned_biomasses;
      tmax = 50
     )
println("$(warmup)")


if last_sim > length(pruned_foodwebs)
      last_sim = length(pruned_foodwebs)
end
println("Running param sim from lines $first_sim to $last_sim")

if last_sim - first_sim > 5000
  throw(DomainError(last_sim - first_sim, "Why are you trying to do more than 5000 simulations? Bad SHARC 🦈"))
end 

timing = @elapsed sim = @showprogress pmap(p ->
                         merge(
                               (rep = p.rep, params = p.params),
                               mega_stoch_sim(p.params, p.pruned_biomasses;
                               tmax = 15000,
                               gc_thre = .02
                               )
                              ),
                              pruned_foodwebs[first_sim:last_sim],
                         batch_size = 50
                        )

println("$(length(sim)) simulations took $(round(timing /60, digits = 2)) minutes to run")

rep = []
A = []
M = []
h = []
ω = []
σ = []
θ = []
biomasses = []
last_extinction_point = []

for i in eachindex(sim)
  push!(rep, sim[i].rep)
  push!(A, sim[i].params.network.A)
  push!(M, sim[i].params.network.M)
  push!(h, sim[i].params.functional_response.h)
  push!(ω, sim[i].params.functional_response.ω)
  push!(σ, sim[i].params.stochasticity.σe)
  push!(θ, sim[i].params.stochasticity.θ)
  push!(biomasses, sim[i].final_biomasses)
  push!(last_extinction_point, sim[i].last_extinction_point)

end

df = DataFrame(rep = rep, A = A, M = M, h = h, ω = ω, σ = σ, θ = θ, biomasses = biomasses, last_extinction_point = last_extinction_point)

file = string("./out/hill_2_stochastic_burnin_", first_sim, "_", last_sim, ".arrow")

Arrow.write(file, df)
