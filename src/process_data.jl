function fetch_faas(fname)

    f = XLSX.readxlsx(fname)

    sh1 = f["Sheet1"]

    pop = []
    t   = Vector{Float64}(sh1["E"][27:285])

    for col in eachcol(sh1["F:CU"])
        dat = Vector{Float64}(collect(skipmissing(col[27:285])))

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
                            cov_τ_decay = 0,
                            cov_Ca_i    = 0,
                            CaM_equil   = true,
                            cov_kon_D  = 8.77e5,         ### in M^-1 ms^-1
                            cov_koff_D =33.2,
                            cov_Kd_Fluo4FF = 23e-6, ## https://pubmed.ncbi.nlm.nih.gov/29604967/
                            cov_kon_Fluo4FF  = 1
                        ),
            observations = (F_F0 = dat, Ca_CaM = repeat([missing], length(dat))),
            time = t[1:length(dat)]
        )
        push!(pop, p)
    end

    return [pop...]
end


function fetch_shifman(fname)

    f = XLSX.readxlsx(fname)

    pop = []
    t   = [0.0; 1.0]

    X = f["Shifman_2006.csv"]["A"][2:end] .* 1e-6;
    Y = [i for i in f["Shifman_2006.csv"]["B"][2:end]];

    for (Ca_i, CaM_i) in zip(X, Y)

        p = Subject(;
            id = Random.randstring(),
            covariates = (cov_f_frac = 1.0, 
                            cov_τ_f= 1.0,
                            cov_τ_s= 1.0,
                            cov_kon_DMn= 1.0,
                            cov_koff_DMn=1.0,
                            cov_Kd_DMn=1.0,
                            cov_koff_PP=1.0,
                            cov_DMn_t=0.0,
                            cov_CaM_t=5e-6,
                            cov_OGB5_t=0.0,
                            cov_Fluo4FF_t=5e-6,
                            cov_Ca_free=Ca_i,
                            PCD=1.0,
                            cov_τ_decay = 0,
                            cov_Ca_i    = 0,
                            CaM_equil   = true,
                            cov_kon_D   = 8.77e5,         ### in M^-1 ms^-1
                            cov_koff_D  = 33.2,
                            cov_Kd_Fluo4FF   = 23e-6, ## https://pubmed.ncbi.nlm.nih.gov/29604967/
                            cov_kon_Fluo4FF  = 1
                        ),
            observations = (F_F0 = [missing; missing], Ca_CaM = [CaM_i; CaM_i]),
            time = t
        )
        push!(pop, p)
    end

    return [pop...]
end


function fetch_crouch(fname)

    f = XLSX.readxlsx(fname)

    pop = []
    t   = [0.0; 1.0]

    X = f["Shifman_2006.csv"]["A"][2:end] .* 1e-6;
    Y = [i for i in f["Shifman_2006.csv"]["B"][2:end]];

    for (Ca_i, CaM_i) in zip(X, Y)

        p = Subject(;
            id = Random.randstring(),
            covariates = (cov_f_frac = 1.0, 
                            cov_τ_f= 1.0,
                            cov_τ_s= 1.0,
                            cov_kon_DMn= 1.0,
                            cov_koff_DMn=1.0,
                            cov_Kd_DMn=1.0,
                            cov_koff_PP=1.0,
                            cov_DMn_t=0.0,
                            cov_CaM_t=5e-6,
                            cov_OGB5_t=0.0,
                            cov_Fluo4FF_t=5e-6,
                            cov_Ca_free=Ca_i,
                            PCD=1.0,
                            cov_τ_decay = 0,
                            cov_Ca_i    = 0,
                            CaM_equil   = true,
                            cov_kon_D   = 8.77e5,         ### in M^-1 ms^-1
                            cov_koff_D  = 33.2,
                            cov_Kd_Fluo4FF   = 23e-6, ## https://pubmed.ncbi.nlm.nih.gov/29604967/
                            cov_kon_Fluo4FF  = 1
                        ),
            observations = (F_F0 = [missing; missing], Ca_CaM = [CaM_i; CaM_i]),
            time = t
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

