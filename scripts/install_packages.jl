import Pkg

try
  Pkg.develop(path = "/users/bi1tma/julia_coding/stoch_befw")
catch
  Pkg.develop(path = "/users/bop21tdm/Documents/julia_coding/stoch_befw")
end

Pkg.add(
	["Distributed","Serialization", "DataFrames", "Sobol",
	 "StatsBase", "Distributions", "ProgressMeter", "Arrow"]
	)
