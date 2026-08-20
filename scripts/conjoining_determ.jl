using Arrow, DataFrames

dat1 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_1_1000.arrow"))
dat2 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_1001_2000.arrow"))
dat3 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_2001_3000.arrow"))
dat4 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_3001_4000.arrow"))
dat5 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_4001_5000.arrow"))
dat6 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_5001_6000.arrow"))
dat7 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_6001_7000.arrow"))
dat8 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_7001_8000.arrow"))
dat9 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_8001_9000.arrow"))
dat10 = DataFrame(Arrow.Table("hill_2_deterministic_burnin_9001_10000.arrow"))

dat = DataFrame()

dat = vcat(dat, dat1)
dat = vcat(dat, dat2)
dat = vcat(dat, dat3)
dat = vcat(dat, dat4)
dat = vcat(dat, dat5)
dat = vcat(dat, dat6)
dat = vcat(dat, dat7)
dat = vcat(dat, dat8)
dat = vcat(dat, dat9)
dat = vcat(dat, dat10)

Arrow.write("hill_2_deterministic_burnin_combined.arrow", dat)
