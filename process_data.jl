function fetch_data(f)
    sh1 = f["Sheet1"]

    pop = []
    t   = Vector{Float64}(sh1["E"][27:285])

    for col in eachcol(sh1["F:CU"])
        dat = Vector{Float64}(collect(skipmissing(col[27:285])))
        batch_id = nothing
        if contains(col[1], "#1021")
            batch_id = 1
        elseif contains(col[1], "#0611")
            batch_id = 2
        elseif contains(col[1], "#1107")
            batch_id = 3
        end

        p = Subject(;
            id = col[1],
            covariates = (cov_f_frac = col[3], 
                            cov_τ_f=col[4],
                            cov_τ_s=col[5],
                            cov_kon_DMn=col[8],
                            cov_koff_DMn=col[9],
                            cov_Kd_DMn=col[7],
                            cov_koff_PP=col[12],
                            cov_DMn_t=col[6],
                            cov_CaM_t=col[21],
                            cov_OGB5_t=col[14],
                            cov_Fluo4FF_t=0,
                            cov_Ca_free=col[20],
                            PCD=col[23],
                            batch_id=batch_id,
                            cov_τ_decay = 0,
                            cov_Ca_i    = 0,
                            CaM_equil   = true
                        ),
            observations = (F_F0 = dat,),
            time = t[1:length(dat)]
        )
        push!(pop, p)
    end

    return [pop...]
end


function subsample_start(pop, slice)

    new_pop = []

    for p in pop
        new_p = Subject(;
            id = p.id,
            covariates = p.covariates.u,
            observations = (F_F0 = [p.observations.F_F0[slice]; p.observations.F_F0[slice[end]+1:end]],),
            time = [p.time[slice]; p.time[slice[end]+1:end]]
        )
        push!(new_pop, new_p)
    end
    return [new_pop...]
end
