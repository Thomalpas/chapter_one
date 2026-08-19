import Pkg

Pkg.develop(path = "/users/bop21tdm/Documents/julia_coding/stoch_befw/")

using stoch_befw

fw = FoodWeb(nichemodel, 10, C = 0.2)

trophic_levels(fw)

params = ModelParameters(fw)

B0 = repeat([0.5], richness(params.network))

simulate(params, B0, tmax = 100)

θ = 0.1
σ = 0.1

stochasticity_object = AddStochasticity(params.network, addstochasticity = false, target = "producers", θ = θ, σe = σ)
no_mortality = BioRates(params.network, d = 0)
params.stochasticity = stochasticity_object
params.biorates = no_mortality

sol = simulate(params, B0, tmax = 100, dt = 0.1)

tidy = tidy_output(params, sol)
