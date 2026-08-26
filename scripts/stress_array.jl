println("Rundir is $(pwd())")

using Distributed, Serialization

first_sim, last_sim = try parse(Int, ARGS[1]), parse(Int, ARGS[2])
catch
    1, 20
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
script_location = try
    ARGS[3]
catch
    "local"
end

if script_location == "HPC"
    flag = "--project=~/julia_coding/chapter_one/" # HPC flag
elseif script_location == "local"
flag = "--project=~/Documents/julia_coding/chapter_one/" # local machine flag
end
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(ncpu - 1, exeflags=flag)
println("Using $(ncpu -1) cores")

@everywhere import Pkg, Random.seed!
@everywhere Pkg.instantiate()
@everywhere using DataFrames, Sobol
@everywhere using stoch_befw
@everywhere using StatsBase, Distributions
@everywhere using ProgressMeter, Arrow

@everywhere include("../src/stress_sim.jl")

@everywhere include("../src/utils.jl")

# Step one, bring in data

hill = 2

#for hill in hill_exponent

  burnin_timeseries = try DataFrame(Arrow.Table("./out/hill_$(hill)_stochastic_burnin_combined.arrow"))
  catch
    DataFrame(Arrow.Table("../out/hill_$(hill)_stochastic_burnin_combined.arrow"))
  end

  biomasses = Vector{Vector{Float64}}(burnin_timeseries.biomasses)

  reps = Vector{Int64}(burnin_timeseries.rep)
  interaction_matrices = Matrix.(reshape.(burnin_timeseries.A, length.(biomasses), length.(biomasses)))
  body_masses = Vector{Float64}.(burnin_timeseries.M)
  h = Float64.(burnin_timeseries.h)
  ω_matrices = Matrix{Float64}.(reshape.(burnin_timeseries.ω, length.(biomasses), length.(biomasses)))
  θ = first.(burnin_timeseries.θ)
  σ = first.(burnin_timeseries.σ)

  model_params = []

  for i in eachindex(burnin_timeseries.A)

    fw = FoodWeb(interaction_matrices[i], M = body_masses[i])
    br = BioRates(fw, d = 0)
    fr = BioenergeticResponse(fw, h = h[i], ω = ω_matrices[i])
    as = AddStochasticity(fw, addstochasticity = true, target = "producers", θ = θ[i], σe = σ[i])
    s = Stressor(fw, addstressor = true, slope = -0.0001, start = 6400.0)
    MP = ModelParameters(fw, biorates = br, functional_response = fr, stochasticity = as, stressor = s)
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

  pm = sample(pruned_foodwebs)
  println("Running warmup")

  warmup = stress_sim(pm.params, pm.pruned_biomasses;
        tmax = 50
      )
  
  println("$(warmup)")

  if last_sim > length(pruned_foodwebs)
        global last_sim = length(pruned_foodwebs)
  end
  println("Running param sim from lines $first_sim to $last_sim")

  if last_sim - first_sim > 5000
    throw(DomainError(last_sim - first_sim, "Why are you trying to do more than 5000 simulations? Bad SHARC 🦈"))
  end 

  timing = @elapsed sim = @showprogress pmap(p ->
                          merge(
                                (rep = p.rep, params = p.params),
                                stress_sim(p.params, p.pruned_biomasses;
                                tmax = 16400,
                                gc_thre = .02
                                )
                                ),
                                pruned_foodwebs[first_sim:last_sim],
                          batch_size = 1
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
  first_extinction_point = []

  for i in eachindex(sim)
    push!(rep, sim[i].rep)
    push!(A, sim[i].params.network.A)
    push!(M, sim[i].params.network.M)
    push!(h, sim[i].params.functional_response.h)
    push!(ω, sim[i].params.functional_response.ω)
    push!(σ, sim[i].params.stochasticity.σe)
    push!(θ, sim[i].params.stochasticity.θ)
    push!(biomasses, sim[i].final_biomasses)
    push!(first_extinction_point, sim[i].first_extinction_point)

  end

  df = DataFrame(rep = rep, A = A, M = M, h = h, ω = ω, σ = σ, θ = θ, biomasses = biomasses, first_extinction_point = first_extinction_point)

  if script_location == "HPC"
    file = string("../../../../mnt/parscratch/users/bi1tma/hill_$(hill)_stress_timeseries_", first_sim, "_", last_sim, ".arrow")
  elseif script_location == "local"
    file = string("./out/hill_$(hill)_stress_timeseries_", first_sim, "_", last_sim, ".arrow")
  end

  Arrow.write(file, df)

#end
