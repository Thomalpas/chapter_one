println("Rundir is $(pwd())")

import Pkg

try
  Pkg.develop(path = "/users/bi1tma/julia_coding/stoch_befw")
catch
  Pkg.develop(path = "/users/bop21tdm/Documents/julia_coding/stoch_befw")
end

Pkg.instantiate()

using Distributed, Serialization

first_sim, last_sim = try parse(Int, ARGS[1]), parse(Int, ARGS[2])
catch
    1, 100
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
flag = "--project=~/chapter_one/"
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

using Arrow

#@everywhere using DifferentialEquations, SparseArrays
#@everywhere using LinearAlgebra, DataFrames
#@everywhere using Distributions, ProgressMeter
#@everywhere using StatsBase, Arrow
#@everywhere using JSON3, FiniteDiff

@everywhere include("../src/determ_sim.jl")

@everywhere include("../src/utils.jl")

# Step one, bring in data
hill_exponent = 2

p_vec = generate_parameters(n = 10000)

model_params = []
for i in eachindex(p_vec.n_species)
    fw = FoodWeb(nichemodel, p_vec.:n_species[i], C = p_vec.:C[i], Z = p_vec.:Z[i])
    br = BioRates(fw, d = 0)
    fr = BioenergeticResponse(fw, h = hill_exponent)
    MP = ModelParameters(fw, biorates = br, functional_response = fr)
    push!(model_params, MP)
end

reps = [1:1:length(model_params);]

data_vec = []

for i in eachindex(reps)
    push!(data_vec, (rep = reps[i], params = model_params[i]))
end


pm = sample(data_vec)
println("Running warmup")

warmup = determ_sim(pm.params;
      tmax = 50
     )
println("$(warmup)")


if last_sim > length(model_params)
      last_sim = length(model_params)
end
println("Running param sim from lines $first_sim to $last_sim")

if last_sim - first_sim > 1000
  throw(DomainError(last_sim - first_sim, "Why are you trying to do more than 1000 simulations? Bad SHARC 🦈"))
end 

timing = @elapsed sim = @showprogress pmap(p ->
                         merge(
                               (rep = p.rep, params = p.params),
                               determ_sim(p.params;
                               tmax = 100000,
                               gc_thre = .02
                               )
                              ),
                              data_vec[first_sim:last_sim],
                         batch_size = 10
                        )

println("$(length(sim)) simulations took $(round(timing /60, digits = 2)) minutes to run")

rep = []
A = []
M = []
h = []
ω = []
biomasses = []
last_extinction_point = []

for i in eachindex(sim)
  push!(rep, sim[i].rep)
  push!(A, sim[i].params.network.A)
  push!(M, sim[i].params.network.M)
  push!(h, sim[i].params.functional_response.h)
  push!(ω, sim[i].params.functional_response.ω)
  push!(biomasses, sim[i].final_biomasses)
  push!(last_extinction_point, sim[i].last_extinction_point)
end

df = DataFrame(rep = rep, A = A, M = M, h = h, ω = ω, biomasses = biomasses, last_extinction_point = last_extinction_point)

file = string("hill_", hill_exponent, "_deterministic_burnin_", first_sim, "_", last_sim, ".arrow")

Arrow.write(file, df)
