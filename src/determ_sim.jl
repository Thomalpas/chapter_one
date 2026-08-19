function determ_sim(params;
    tmax = 100000,
    gc_thre = .02)

    if rand(Distributions.Uniform(0, 1)) < gc_thre
        println("")
        GC.gc()
        try ccall(:malloc_trim, Cvoid, (Cint,), 0)
        catch
        end
        GC.safepoint()
    end

    B0 = repeat([0.5], richness(params.network))

    burnin_timeseries = try simulate(params, B0, tmax = tmax, verbose = false)
    catch
        burnin_timeseries = missing
    end

    if !ismissing(burnin_timeseries)

        timeseries_matrix = reduce(hcat, burnin_timeseries.u)'

        last_extinction_point = 0
        for i in 1:richness(params.network)
            if any(timeseries_matrix[:,i] .== 0)
                extinction_point = length(findall(x -> x != 0, timeseries_matrix[:,i])) + 1
                if extinction_point > last_extinction_point
                    last_extinction_point = extinction_point
                end
            end
        end

        biomasses = tidy_output(params, burnin_timeseries).B[end, :]
    else
        final_biomasses = missing
        last_extinction_point = missing
    end

    out =  merge(
        (final_biomasses = biomasses, last_extinction_point = last_extinction_point)
        )
    out
end
