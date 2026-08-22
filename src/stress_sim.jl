function stress_sim(params,  B0;
    tmax = 11000,
    gc_thre = .02)

    if rand(Distributions.Uniform(0, 1)) < gc_thre
        println("")
        GC.gc()
        GC.safepoint()
    end

    stress_timeseries = try simulate(params, B0, tmax = tmax, verbose = false, dt = 0.1)
    catch
        stress_timeseries = missing
    end

    if !ismissing(stress_timeseries)
        
        timeseries_matrix = reduce(hcat, stress_timeseries.u)'

        first_extinction_point = 15000
        for i in 1:richness(params.network)
            if any(timeseries_matrix[:,i] .== 0)
                extinction_point = length(findall(x -> x != 0, timeseries_matrix[:,i])) + 1
                if extinction_point < first_extinction_point
                    first_extinction_point = extinction_point
                end
            end
        end

        tidy = tidy_output(params, stress_timeseries)

        biomasses = tidy.B[findall(x -> x == 1, rem.(round.(tidy.t, sigdigits = 6), 1) .== 0), :]
        
        biomasses = round.(biomasses, sigdigits = 6)

    else
        biomasses = missing
        first_extinction_point = missing
    end

    out =  merge(
        (final_biomasses = biomasses, first_extinction_point = first_extinction_point)
        )
    out
end