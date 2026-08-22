
using Arrow, DataFrames

hill = 2

dat = DataFrame()

for start in [1:1000:89001;]
    if start == 89001
        temp = DataFrame(Arrow.Table("./out/hill_$(hill)_stochastic_burnin_$(start)_89568.arrow"))
    else
        temp = DataFrame(Arrow.Table("./out/hill_$(hill)_stochastic_burnin_$(start)_$(start+999).arrow"))
    end
    dat = vcat(dat, temp)
end

Arrow.write("./out/hill_$(hill)_stochastic_burnin_combined.arrow", dat)
