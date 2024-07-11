#############################
### SOME HELPER FUNCTIONS ###
#############################

###function idx_to_groups(full_pop, idxs)
###
###    A = full_pop[idxs[idxs .> 0  .&& idxs .<= 13]]
###    B = full_pop[idxs[idxs .> 13 .&& idxs .<= 24]]
###    C = full_pop[idxs[idxs .> 24 .&& idxs .<= 35]]
###    D = full_pop[idxs[idxs .> 35 .&& idxs .<= 51]]
###    E = full_pop[idxs[idxs .> 51 .&& idxs .<= 65]]
###    F = full_pop[idxs[idxs .> 65 .&& idxs .<= 80]]
###    G = full_pop[idxs[idxs .> 81 .&& idxs .<= 94]]
###
###    return [A, B, C, D, E, F, G]
###end
###
###function process_fpms(fpms, pop)
###    all_pre = []
###    all_sims = []
###   
###    for j in 1:length(fpms)
###    
###        fpm = fpms[j]
###   
###        ebes = empirical_bayes(fpm.model, pop, coef(fpm), LaplaceI(), diffeq_options=(;alg=Rodas5P(), abstol=1e-16))
###        pre_i = [fpm.model.pre(coef(fpm), ebes[i], pop[i]).ret for i in 1:length(pop)]
###        sim_i = simobs(fpm.model, 
###            pop,
###            coef(fpm),
###            ebes,
###            diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
###        )
###    
###        push!(all_pre, pre_i)
###        push!(all_sims, sim_i)
###    end
###    return all_pre, all_sims
###end
###
###
###function gen_dummy_data(Ca_range)
###
###    pop = []
###
###    for val in Ca_range
###        sub = Subject(;
###            id = string(val),
###            covariates = (cov_f_frac    = 0, 
###                            cov_τ_f     = 1,
###                            cov_τ_s     = 1,
###                            cov_kon_DMn = 0,
###                            cov_koff_DMn= 0,
###                            cov_Kd_DMn  = 0,
###                            cov_koff_PP = 0,
###                            cov_DMn_t   = 0,
###                            cov_CaM_t   = 5e-6,
###                            cov_OGB5_t  = 0,
###                            cov_Fluo4FF_t = 5e-6,
###                            cov_Ca_free = val,
###                            PCD         = 0,
###                            batch_id    = 1,
###                            cov_τ_decay = 0,
###                            cov_Ca_i    = 0,
###                            cov_kon_D   = 1.0,
###                            cov_koff_D  = 1.0,
###                            cov_Kd_Fluo4FF  = 23e-6,
###                            cov_kon_Fluo4FF = 1.0,
###                            CaM_equil   = true
###                            
###            ), 
###            observations = (F_F0 = [1.0;],),
###            time = [1.0;]
###        )
###        push!(pop, sub)
###    end
###
###    return [pop...]
###end
###
###
###function get_bwell_eqs(fpm, pop)
###    no_bound    = []
###    part_bound  = []
###    full_bound  = []
###
###    for p in pop
###        pre = fpm.model.pre(coef(fpm), (η = [0.0; 0.0; 0.0],), p).ret
###
###        push!(no_bound,     pre.CaM₀_all[1] / pre.CaM_t * 0)
###        push!(part_bound,   pre.CaM₀_all[2] / pre.CaM_t * 2)
###        push!(full_bound,   pre.CaM₀_all[3] / pre.CaM_t * 4)
###    end
###
###    return no_bound, part_bound, full_bound
###end
###
###
###function get_bhalla_eqs(fpm, pop)
###    no_bound    = []
###    part_bound  = []
###    full_bound  = []
###
###    for p in pop
###        pre = fpm.model.pre(coef(fpm), (η = [0.0; 0.0; 0.0],), p).ret
###
###        push!(no_bound,     pre.CaM₀_all[1] / pre.CaM_t * 0)
###        push!(part_bound,   2 * pre.CaM₀_all[2] / pre.CaM_t 
###                          + 3 * pre.CaM₀_all[3] / pre.CaM_t)
###        push!(full_bound,   4 * pre.CaM₀_all[4] / pre.CaM_t)
###    end
###
###    return no_bound, part_bound, full_bound
###end
###
###
###function get_shifman_eqs(fpm, pop)
###    no_bound    = []
###    part_bound  = []
###    full_bound  = []
###
###    for p in pop
###        pre = fpm.model.pre(coef(fpm), (η = [0.0; 0.0; 0.0],), p).ret
###
###        push!(no_bound,     0 * pre.CaM₀_all[1] / pre.CaM_t)
###        push!(part_bound,   1 * pre.CaM₀_all[2] / pre.CaM_t 
###                          + 2 * pre.CaM₀_all[3] / pre.CaM_t
###                          + 3 * pre.CaM₀_all[4] / pre.CaM_t)
###        push!(full_bound,   4 * pre.CaM₀_all[5] / pre.CaM_t)
###    end
###
###    return no_bound, part_bound, full_bound
###end
###
###
###function get_pepke_m2_eqs(fpm, pop)
###    no_bound    = []
###    part_bound  = []
###    full_bound  = []
###
###    for p in pop
###        pre = fpm.model.pre(coef(fpm), (η = [0.0; 0.0; 0.0],), p).ret
###
###        push!(no_bound,     pre.c0[1] / pre.CaM_t * 0)
###        push!(part_bound,   (pre.c0[2] / pre.CaM_t + pre.c0[2] / pre.CaM_t))
###        push!(full_bound,   pre.c0[4] / pre.CaM_t * 4)
###    end
###
###    return no_bound, part_bound, full_bound
###end
###
###
###function get_byrne_eqs(fpm, pop)
###    no_bound    = []
###    part_bound  = []
###    full_bound  = []
###
###    for p in pop
###        pre = fpm.model.pre(coef(fpm), (η = [0.0; 0.0; 0.0],), p).ret
###
###        push!(no_bound,     0)
###        b1 = (pre.N₀[2] + pre.C₀[2]) / 2pre.CaM_t * 1
###        b2 = (pre.N₀[3] + pre.C₀[3]) / 2pre.CaM_t * 2
###        push!(part_bound,   b1+b2)
###        push!(full_bound,   (pre.N₀[4] + pre.C₀[4]) / 2pre.CaM_t * 4)
###    end
###
###    return no_bound, part_bound, full_bound
###end
###
###
###function get_faas_eqs(fpm, pop)
###    no_bound    = []
###    part_bound  = []
###    full_bound  = []
###
###    for p in pop
###
###        pre = fpm.model.pre(coef(fpm), (η = [0.0],), p).ret
###
###        CaM₀  = 0 
###        CaMp₀ = pre.N₀_all[2] / pre.CaM_t * 1 + 
###                pre.N₀_all[3] / pre.CaM_t * 2 + 
###                pre.C₀_all[2] / pre.CaM_t * 1 + 
###                pre.C₀_all[3] / pre.CaM_t * 2 
###        CaMf₀ = 0
###        push!(no_bound,     CaM₀)
###        push!(part_bound,   CaMp₀)
###        push!(full_bound,   CaMf₀)
###    end
###
###    return no_bound, part_bound, full_bound
###end
###
###
###function get_multiple_eqs(fpms, pop, single_func)
###
###    no_bound    = []
###    part_bound  = []
###    full_bound  = []
###
###    for fpm in fpms
###        n, p, f = single_func(fpm, pop)
###        push!(no_bound, n)
###        push!(part_bound, p)
###        push!(full_bound, f)
###    end
###
###    return no_bound, part_bound, full_bound
###end
###
###
###function aic_boxplots(triples)
###
###    all_x = []
###    all_y = []
###
###    for trip in triples
###        fpm = trip[1]
###        pop = vcat(trip[2]...)
###        name = trip[3]
###        
###        y_i = my_aic(fpm; pop=pop)
###
###        push!(all_x, name)
###        push!(all_y, y_i)
###    end
###
###    return all_x, all_y
###
###end
###
###
###function gen_f_data(f, amp)
###
###    sub = Subject(;
###        id = "1",
###        ###events = DosageRegimen([DosageRegimen(35e-6/f, time=i) for i in LinRange(20, 1020, f)]...),
###        events = DosageRegimen([DosageRegimen(amp, time=i) for i in LinRange(20, 1020, f)]...),
###        covariates = (cov_f_frac    = 0, 
###                        cov_τ_f     = 1,
###                        cov_τ_s     = 1,
###                        cov_kon_DMn = 0,
###                        cov_koff_DMn= 0,
###                        cov_Kd_DMn  = 0,
###                        cov_koff_PP = 0,
###                        cov_DMn_t   = 0,
###                        cov_CaM_t   = 20e-6,
###                        cov_OGB5_t  = 0,
###                        cov_Fluo4FF_t  = 0,
###                        cov_Ca_free = 100e-9,
###                        PCD         = 0,
###                        batch_id    = 1,
###                        cov_τ_decay = 0.0833,
###                        cov_Ca_i    = 100e-9,
###                        cov_kon_D   = 1.0,
###                        cov_koff_D  = 1.0,
###                        cov_Kd_Fluo4FF  = 23e-6,
###                        cov_kon_Fluo4FF = 1.0,
###                        CaM_equil   = true
###
###        ), 
###        observations = (F_F0 = repeat([missing], 10001),),
###        time = collect(LinRange(0, 1500, 10001))
###    )
###
###    return [sub;]
###end
###
###
###function calc_integration(fpms, f_range, amp, getter)
###
###    lines = map(fpms) do j
###       
###        partial_j   = Vector{Float64}()
###        full_j      = Vector{Float64}()
###
###        for ff in f_range
###            d_ff = gen_f_data(ff, amp)
###
###            _, p_ff, f_ff = getter(j, d_ff)
###            i_p = integrate(d_ff[1].time, p_ff)
###            i_f = integrate(d_ff[1].time, f_ff)
###
###            push!(partial_j, i_p)
###            push!(full_j, i_f)
###        end
###
###        (partial_j, full_j)
###    end
###
###    return lines 
###end
###
###
###function get_bwell_dyns(fpm, pop)
###
###    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]
###
###    return abs.(obs.dynamics.CaM), abs.(obs.dynamics.CaM2Ca), abs.(obs.dynamics.CaM4Ca)
###end
###
###
###function get_bhalla_dyns(fpm, pop)
###
###    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]
###
###    CaM0 = abs.(obs.dynamics.CaM)
###    CaMp = abs.(obs.dynamics.CaM2Ca) + abs.(obs.dynamics.CaM3Ca)
###    CaMf = abs.(obs.dynamics.CaM4Ca)
###
###    return CaM0, CaMp, CaMf
###end
###
###
###function get_shifman_dyns(fpm, pop)
###
###    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]
###
###    CaM0 = abs.(obs.dynamics.CaM)
###    CaMp = abs.(obs.dynamics.CaM1Ca) + abs.(obs.dynamics.CaM2Ca) + abs.(obs.dynamics.CaM3Ca)
###    CaMf = abs.(obs.dynamics.CaM4Ca)
###
###    return CaM0, CaMp, CaMf
###end
###
###
###function get_bwell_CN_dyns(fpm, pop)
###
###    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]
###
###    CaM0 = abs.(obs.dynamics.CaM0)
###    CaMp = abs.(obs.dynamics.CaM2C) + abs.(obs.dynamics.CaM2N)
###    CaMf = abs.(obs.dynamics.CaM4)
###
###    return CaM0, CaMp, CaMf
###end
###
###
###function get_faas_dyns(fpm, pop)
###
###    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]
###
###    CaM0 = (abs.(obs.dynamics.NtNt) + abs.(obs.dynamics.CtCt)) / 2
###    CaMp = (abs.(obs.dynamics.NtNr) + abs.(obs.dynamics.CtCr)) / 2
###    CaMf = (abs.(obs.dynamics.NrNr) + abs.(obs.dynamics.CrCr)) / 2
###
###    return CaM0, CaMp, CaMf
###end
###
###
###function get_pepke_m1_dyns(fpm, pop)
###
###    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]
###
###    CaM0 = (abs.(obs.dynamics.N0) + abs.(obs.dynamics.C0)) / 2
###    CaMp = (abs.(obs.dynamics.N1) + abs.(obs.dynamics.C1)) / 2
###    CaMf = (abs.(obs.dynamics.N2) + abs.(obs.dynamics.C2)) / 2
###
###    return CaM0, CaMp, CaMf
###end
###
###
###function get_byrne_dyns(fpm, pop)
###
###    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]
###    d = obs.dynamics
###
###    CaM0 = (abs.(d.N0) + abs.(d.C0)) / 2
###    CaMp = (abs.(d.N1) + abs.(d.C1) + abs.(d.N2) + abs.(d.C2)) / 4
###    CaMf = (abs.(d.N3) + abs.(d.C3)) / 2
###
###    return CaM0, CaMp, CaMf
###end
###
###
###function get_coef_lists(fpms, symbols; 
###    sort=false, 
###    kd_C_name=nothing,
###    kd_N_name=nothing,
###    counternames=nothing)
###
###    params = []
###
###    for s in symbols
###
###        p_s = map(fpms) do f
###           coef(f)[s] 
###        end
###
###        if length(unique(p_s)) == 1
###            p_s = p_s[1:1]
###        end
###
###        push!(params, p_s)
###    end
###
###    nt = NamedTuple{symbols}(params)
###
###    if sort
###        kd_c = nt[kd_C_name]
###        kd_n = nt[kd_N_name]
###        
###        repl_idx = kd_c .> kd_n
###
###        d = Dict()
###
###        for (n_i, n_j) in counternames
###            temp_n_i = nt[n_i]
###            temp_n_j = nt[n_j]
###
###            i_to_j = temp_n_i[repl_idx]
###            j_to_i = temp_n_j[repl_idx]
###
###            temp_n_i[repl_idx] .= j_to_i
###            temp_n_j[repl_idx] .= i_to_j
###
###            d[n_i] = temp_n_i
###            d[n_j] = temp_n_j
###        end
###
###        ret_nt = NamedTuple(d)
###    else
###        ret_nt = nt
###    end
###
###    return ret_nt
###end
###
###
###function get_ppcor_matrix(nt)
###    
###    names    = keys(nt)
###    name_mat = Iterators.product(names, names)
###    corr_mat = zeros(size(name_mat))
###    pval_mat = zeros(size(name_mat))
###
###    map(name_mat) do (i,j)
###
###        x = nt[i]
###        y = nt[j]
###        Z = hcat(nt[setdiff(names, (i,j))]...)
###        res = CorrelationTest(x, y, Z)
###
###        i_idx = findfirst(names .== i)
###        j_idx = findfirst(names .== j)
###
###        corr_mat[i_idx, j_idx] = res.r
###        pval_mat[i_idx, j_idx] = HypothesisTests.pvalue(res) 
###    end
###
###    return name_mat, corr_mat, pval_mat
###end
###
###
###function byrne_N_C_lobe_binding_comp(fpms)
###
###    all_C0_c = []
###    all_C0_t = []
###    all_C1_c = []
###    all_C1_t = []
###    all_C2_c = []
###    all_C2_t = []   
###    all_C3_c = []
###    all_C3_t = []
###
###    all_N0_c = []
###    all_N0_t = []
###    all_N1_c = []
###    all_N1_t = []
###    all_N2_c = []
###    all_N2_t = []
###    all_N3_c = []
###    all_N3_t = []
###
###    all_catot = Vector{Float64}()
###
###    idxs = occursin.("#", [i.id for i in fpms[1].data])
###
###    for f in fpms
###        sims = simobs(f)[1:sum(idxs)]
###        for s in sims
###
###            CaM_t = s.icoefs.CaM_t[1]
###            ss = s.sol
###            tt = s.time
###
###            Ca_tot = ss[1,1] + ss[3,  1] + ss[5,  1] + ss[7,  1] + ss[9,  1] + 
###                   ss[11, 1] + ss[13, 1] + ss[14, 1] +2ss[15, 1] + ss[17, 1] + ss[18, 1] + 2ss[19, 1]
###
###            push!(all_catot, Ca_tot)
###
###            push!(all_N0_c, ss[end-7,:] ./ CaM_t)
###            push!(all_N1_c, ss[end-6,:] ./ CaM_t)
###            push!(all_N2_c, ss[end-5,:] ./ CaM_t)
###            push!(all_N3_c, ss[end-4,:] ./ CaM_t)
###
###            push!(all_C0_c, ss[end-3,:] ./ CaM_t)
###            push!(all_C1_c, ss[end-2,:] ./ CaM_t)
###            push!(all_C2_c, ss[end-1,:] ./ CaM_t)
###            push!(all_C3_c, ss[end-0,:] ./ CaM_t)
###
###            push!(all_N0_t, tt)
###            push!(all_N1_t, tt)
###            push!(all_N2_t, tt)
###            push!(all_N3_t, tt)
###
###            push!(all_C0_t, tt)
###            push!(all_C1_t, tt)
###            push!(all_C2_t, tt)
###            push!(all_C3_t, tt)
###        end
###    end
###
###    return (n0 = all_N0_c, n1 = all_N1_c, n2 = all_N2_c, n3 = all_N3_c, 
###            c0 = all_C0_c, c1 = all_C1_c, c2 = all_C2_c, c3 = all_C3_c), 
###           (n0 = all_N0_t, n1 = all_N1_t, n2 = all_N2_t, n3 = all_N3_t, 
###            c0 = all_C0_t, c1 = all_C1_t, c2 = all_C2_t, c3 = all_C3_t), 
###            all_catot
###end
###
###
###function get_param_lists(fpms, symbols; 
###    sort=false, 
###    sort_C_names=nothing,
###    sort_N_names=nothing,
###    counternames=nothing)
###
###    params = []
###
###    for s in symbols
###        p_s = Float64[]
###
###        for f in fpms
###            pre = f.model.pre(coef(f), (η = [0.0; 0.0; 0.0],), f.data[1]).ret
###            push!(p_s, pre[s])
###        end
###
###        if length(unique(p_s)) == 1
###            p_s = p_s[1:1]
###        end
###
###        push!(params, p_s)
###    end
###
###    nt = NamedTuple{symbols}(params)
###
###    if sort
###        kd_c = nt[sort_C_names.off] ./ nt[sort_C_names.on]
###        kd_n = nt[sort_N_names.off] ./ nt[sort_N_names.on]
###        
###        repl_idx = kd_c .> kd_n
###
###        d = Dict()
###
###        for (n_i, n_j) in counternames
###            temp_n_i = nt[n_i]
###            temp_n_j = nt[n_j]
###
###            i_to_j = temp_n_i[repl_idx]
###            j_to_i = temp_n_j[repl_idx]
###
###            temp_n_i[repl_idx] .= j_to_i
###            temp_n_j[repl_idx] .= i_to_j
###
###            d[n_i] = temp_n_i
###            d[n_j] = temp_n_j
###        end
###
###        ret_nt = NamedTuple(d)
###    else
###        ret_nt = nt
###    end
###
###    return ret_nt
###end
###
###
##########################
###### PLOTING PROPPER ###
##########################
###
###
###function plot_by_condition(fpms; pops=nothing)
###
###    if pops === nothing
###        all_preds = vcat(predict.(fpms)...)
###    else
###        all_preds = vcat([predict(i,j) for (i,j) in zip(fpms, pops)]...)
###    end
###
###    uids = unique([split(i.id, "_")[1] * "_" * split(i.id, "_")[2] for i in fpms[1].data if occursin("#", i.id)])
###
###    fig = CairoMakie.Figure(size=(600, 800), fontsize=8)
###
###    for i in 1:length(uids)
###        idx = uids[i]
###
###        row_i = div(i,5) + 1
###        col_i = i - 4*(row_i - 1)
###
###        ax_i = CairoMakie.Axis(fig[row_i, col_i], title=idx)
###
###        preds_id = [j for j in all_preds if occursin(idx, j.subject.id)]
###
###        for pp in preds_id
###            CairoMakie.lines!(ax_i, pp.subject.time, pp.subject.observations.F_F0, color=:black)
###            CairoMakie.lines!(ax_i, pp.subject.time, pp.ipred.F_F0, color=(:red, 0.7), linewidth=3)
###        end
###    end
###
###    if pops === nothing
###        rmses = [mean(rmse(i)) for i in fpms]
###    else
###        rmses = [mean(rmse(i;pop=j)) for (i,j) in zip(fpms, pops)]
###    end
###
###    rmse_str = reduce(*, ["Run " * string(k) * " RMSE: " * string(rmses[k]) * "\n" for k in 1:length(rmses)])
###    ax_text = CairoMakie.Axis(fig[2,4], limits=(0, 1, 0, 1))
###    CairoMakie.hidedecorations!(ax_text)
###    CairoMakie.hidespines!(ax_text)
###    CairoMakie.text!(0.0, 1.0, text=rmse_str, align=(:left, :top), size=10)
###
###    return fig
###end
###
###
###
###
###

###
###
###
###function violin_plot(test_nt, gs; w=17, h=7.5, save=false)
###
###    size_cm     = (17, 7.5)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    
###    f = CairoMakie.Figure(resolution=size_units, fontsize=gs.fontsize)
###
###    ax11 = CairoMakie.Axis(f[1,1], ylabel="Test RMSE", 
###        xticks=(1.0:1.0:13.0, ["Scheme 1 + our rates";
###                      "Scheme 1 + Kim et al.";
###                      "Scheme 2 + our rates";
###                      "Scheme 2 + Bhalla and Iyengar";
###                      "Scheme 3 + our rates";
###                      "Scheme 3 + Shifman et al. + our rates";
###                      "Scheme 4 + our rates";
###                      "Scheme 4 + Pepke et al.";
###                      "Scheme 5 + our rates";
###                      "Scheme 5 + Faas et al.";
###                      "Scheme 5 + Pepke et al.";
###                      "Scheme 6 + our rates";
###                      "Scheme 6 + Byrne et al.";
###        ]),
###        xlabel="Scheme + parameters",
###        xticklabelrotation=pi/9,
###        xticklabelpad=5,
###        xlabelpadding=15,
###        xticklabelsize=8,
###        yticklabelsize=8,
###        xgridvisible=false,
###    )
###    h_111 = CairoMakie.violin!(ax11, 
###        [
###            1 * 1.0 *  ones(length(test_nt[:blackwell]));
###            2 * 1.0 *  ones(length(test_nt[:blackwell_fixed]));
###            3 * 1.0 *  ones(length(test_nt[:bhalla]));
###            4 * 1.0 *  ones(length(test_nt[:bhalla_fixed]));
###            5 * 1.0 *  ones(length(test_nt[:shifman]));
###            6 * 1.0 *  ones(length(test_nt[:shifman_fixed]));
###            7 * 1.0 *  ones(length(test_nt[:pepke_m2]));
###            8 * 1.0 *  ones(length(test_nt[:pepke_m2_fixed]));
###            9 * 1.0 *  ones(length(test_nt[:faas]));
###           10 * 1.0 *  ones(length(test_nt[:faas_fixed]));
###           11 * 1.0 *  ones(length(test_nt[:pepke_m1_fixed]));
###           12 * 1.0 *  ones(length(test_nt[:byrne]));
###           13 * 1.0 *  ones(length(test_nt[:byrne_fixed]))
###        ], 
###        [
###            test_nt[:blackwell];
###            test_nt[:blackwell_fixed];
###            test_nt[:bhalla];
###            test_nt[:bhalla_fixed];
###            test_nt[:shifman];
###            test_nt[:shifman_fixed];
###            test_nt[:pepke_m2];
###            test_nt[:pepke_m2_fixed];
###            test_nt[:faas];
###            test_nt[:faas_fixed];
###            test_nt[:pepke_m1_fixed];
###            test_nt[:byrne];
###            test_nt[:byrne_fixed]
###        ],
###        color = [
###            repeat([gs.our_col], length(test_nt[:blackwell]));
###            repeat([gs.kim_col], length(test_nt[:blackwell_fixed]));
###            repeat([gs.our_col], length(test_nt[:bhalla]));
###            repeat([gs.bhalla_col], length(test_nt[:bhalla_fixed]));
###            repeat([gs.our_col], length(test_nt[:shifman]));
###            repeat([gs.shifman_col], length(test_nt[:shifman_fixed]));
###            repeat([gs.our_col], length(test_nt[:pepke_m2]));
###            repeat([gs.pepke_col], length(test_nt[:pepke_m2_fixed]));
###            repeat([gs.our_col], length(test_nt[:faas]));
###            repeat([gs.faas_col], length(test_nt[:faas_fixed]));
###            repeat([gs.pepke_col], length(test_nt[:pepke_m1_fixed]));
###            repeat([gs.our_col], length(test_nt[:byrne]));
###            repeat([gs.byrne_col], length(test_nt[:byrne_fixed]));
###        ]
###    )
###    
###    if save
###        CairoMakie.save("plots/violin.png", f, pt_per_unit=1, 
###            px_per_unit=px_per_unit, size=size_units)
###    end
###
###    return f
###end
###
###
###function big_plot(fpms_nt, faas_data, val_idxs, gs; w=17, h=22.23, save=false)
###
###
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###
###    colors = [gs.our_col; gs.kim_col; gs.our_col; gs.bhalla_col; 
###              gs.our_col; gs.shifman_col; gs.our_col; gs.pepke_col; 
###              gs.our_col; gs.faas_col; gs.pepke_col; gs.our_col; gs.byrne_col];
###    rowslist = collect(1:13)
###    names = ["Our rates"; "Kim et al.";
###            "Our rates"; "Hayer and Bhalla";
###            "Our rates"; "Shifman et al.";
###            "Our rates"; "Peple et al.";
###            "Our rates"; "Faas et al."; "Pepke et al.";
###            "Our rates"; "Byrne et al."
###    ]
###
###    fpms = [fpms_nt[key][1] for key in keys(fpms_nt)]
###    pop = idx_to_groups(faas_data, val_idxs[1]);
###
###    title_dict = Dict("#1021_WT"=>"Group A", 
###                      "#0611_WTa"=>"Group B",
###                      "#0611_WTb"=>"Group C",
###                      "#1107_WTa"=>"Group D",
###                      "#1107_WTb"=>"Group E",
###                      "#1107_WTc"=>"Group F",
###                      "#1107_WTd"=>"Group G",
###    )
###
###    f = CairoMakie.Figure(resolution=size_units, fontsize=gs.fontsize)
###
###    for i in 1:length(fpms)
###
###        for j in 1:(length(pop))
###
###            pop_j = pop[j]
###
###            splits = split(pop_j[1].id, "_")
###            title_i = splits[1] * "_" * splits[2]
###            title_i = title_dict[title_i]
###
###            ax = CairoMakie.Axis(f[rowslist[i], j], 
###                limits=(nothing, 57, -2, nothing))
###            if i == 13
###                if j == 1
###                    hidedecorations!(ax, ticks=false, ticklabels=false)
###                    ax.xticksize=4
###                    ax.yticksize=4
###                    ax.xticks = [0, 35]
###                    ax.yticks = [0, 5]
###
###                else
###                    hidedecorations!(ax, ticks=false, ticklabels=false)
###                    ax.xticksize=4
###                    ax.yticksize=4
###                    ax.xticks = [0, 35]
###                    ax.yticks = [0, 5]
###                    ax.yticklabelsvisible = false
###                end
###                hidespines!(ax, :t, :r)
###            else
###                hidedecorations!(ax)
###                hidespines!(ax)
###            end
###
###            pre, sims = process_fpms(fpms[i:i], pop_j)
###
###            for sim_i in sims[1]
###                lines!(ax, sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", linewidth=1)
###            end
###
###            d = zip(pre[1], sims[1])
###            for (pre,sim) in d
###                pred = (sim.dynamics.OGB5 + 39.364 * sim.dynamics.CaOGB5) / (pre.OGB5₀ + 39.364 * pre.CaOGB5₀)
###                lines!(ax, sim.time, pred, color=colors[i], linewidth=1, label=names[i])
###            end
###        end
###    end
###
###    CairoMakie.Label(f[0,1], "Group A", fontsize=8, font=:bold)
###    CairoMakie.Label(f[0,2], "Group B", fontsize=8, font=:bold)
###    CairoMakie.Label(f[0,3], "Group C", fontsize=8, font=:bold)
###    CairoMakie.Label(f[0,4], "Group D", fontsize=8, font=:bold)
###    CairoMakie.Label(f[0,5], "Group E", fontsize=8, font=:bold)
###    CairoMakie.Label(f[0,6], "Group F", fontsize=8, font=:bold)
###    CairoMakie.Label(f[0,7], "Group G", fontsize=8, font=:bold)
###
###    CairoMakie.text!(f.content[1], 30.5, 5.5, text = "PCD (in μs)", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[1], 45, 4.0, text = "450", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[1], 45, 2.75, text = "410", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[1], 45, 1.5, text = "380", 
###        align = (:center, :center), fontsize=8)
###    
###    CairoMakie.text!(f.content[2], 47, 5.5, text = "460", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[2], 47, 8.5, text = "480", 
###        align = (:center, :center), fontsize=8)
###
###    CairoMakie.text!(f.content[3], 48, 13.5, text = "445", 
###        align = (:center, :center), fontsize=8)   
###
###    CairoMakie.text!(f.content[4], 48, 1.00, text = "410", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[4], 48, 4.00, text = "450", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[4], 48, 6.50, text = "470", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[4], 48, 9.00, text = "480", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[4], 48,11.10, text = "500", 
###        align = (:center, :center), fontsize=8)
###    
###    CairoMakie.text!(f.content[5], 48, 1.50, text = "400", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[5], 48, 7.00, text = "470", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[5], 48, 9.50, text = "480", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[5], 48,12.00, text = "500", 
###        align = (:center, :center), fontsize=8)
###    
###    CairoMakie.text!(f.content[6], 48, 0.00, text = "370", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[6], 48, 3.25, text = "380", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[6], 48, 6.00, text = "420", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[6], 48, 13.0, text = "480", 
###        align = (:center, :center), fontsize=8)
###    
###    CairoMakie.text!(f.content[7], 48, 1.00, text = "390", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[7], 48, 5.00, text = "420", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[7], 48,14.00, text = "460", 
###        align = (:center, :center), fontsize=8)
###    CairoMakie.text!(f.content[7], 48,18.00, text = "470", 
###        align = (:center, :center), fontsize=8)
###    
###    CairoMakie.Legend(f[1:2, 8], [f.content[1].scene.plots[6],  f.content[8].scene.plots[end]], 
###        ["Our rates", "Kim et al. \n rates"], "Scheme 1", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center)
###    CairoMakie.Legend(f[3:4, 8], [f.content[15].scene.plots[end], f.content[22].scene.plots[end]], 
###        ["Our rates", "Hayer and \n Bhalla rates"], "Scheme 2", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center)
###    CairoMakie.Legend(f[5:6, 8], [f.content[29].scene.plots[end], f.content[36].scene.plots[end]], 
###        ["Our rates", "Shifman et al. \n rates"], "Scheme 3", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center)
###    CairoMakie.Legend(f[7:8, 8], [f.content[43].scene.plots[end], f.content[50].scene.plots[end]], 
###        ["Our rates", "Pepke al. \n rates"], "Scheme 4", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center)
###    CairoMakie.Legend(f[9:11, 8], [f.content[57].scene.plots[end], f.content[64].scene.plots[end], f.content[71].scene.plots[end]], 
###        ["Our rates", "Faas et al. \n rates", "Pepke al. \n rates"], "Scheme 5", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center)
###    CairoMakie.Legend(f[12:13, 8], [f.content[78].scene.plots[end], f.content[85].scene.plots[end]], 
###        ["Our rates", "Byrne et al."], "Scheme 6", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center)
###    
###    f_ly = CairoMakie.Label(f[:,0], "ΔF/F₀", rotation=pi/2, fontsize=8)
###    f_lx = CairoMakie.Label(f[rowslist[end]+1,:], "Time (ms)", fontsize=8)
###    
###    
###    for i in 1:7
###        colsize!(f.layout, i, Relative(1/8))
###    end
###    colsize!(f.layout, 8, Relative(1/6))
###    
###    for i in rowslist
###        CairoMakie.rowsize!(f.layout, i, Relative(1/13))
###    end
###    
###    CairoMakie.colgap!(f.layout, 5)
###    CairoMakie.rowgap!(f.layout, 2)
###    
###    CairoMakie.rowgap!(f.layout, 14, 5)
###
###    if save
###        CairoMakie.save("plots/big_plot.png", f, pt_per_unit=1, 
###            px_per_unit=px_per_unit, size=size_units)
###    end
###
###    return f
###
###end
###
###
###
###function plot_CaM_eqs(ca_range, CaM_tup; i=1, j=1, 
###    xlabel = "Free Ca²⁺ (M)",
###    ylabel = L"Bound $\textrm{Ca}^{2+}$ per CaM",
###    title="Scheme 1", color=:red, line_alpha=0.6, ribbon_alpha=0.2, label="Our rates", 
###    f=nothing, new_axis=false, figsize=(800, 600), limits=(0, 1, 0, 1))
###
###    if f === nothing
###        f = CairoMakie.Figure(resolution=figsize, fontsize=8)
###        ax = CairoMakie.Axis(f[i, j], 
###            xticklabelsize=8, yticklabelsize=8,
###            )
###    elseif new_axis
###        ax = CairoMakie.Axis(f[i, j], 
###            xticklabelsize=8, yticklabelsize=8,
###            )
###    else
###        ax = f.content[end]
###    end
###    
###    if ylabel !== nothing
###        ax.ylabel = ylabel
###    end
###    if xlabel !== nothing
###        ax.xlabel = xlabel
###    end
###    
###    c_sum = CaM_tup[2] .+ CaM_tup[3]
###    c_mat = hcat(c_sum...)
###    
###    if allequal(eachcol(c_mat))
###        med_c = c_mat[:,1]
###        q_025 = c_mat[:,1]
###        q_975 = c_mat[:,1]
###    else
###        med_c = median.(eachrow(c_mat))
###        q_025 = quantile.(eachrow(c_mat), 0.025)
###        q_975 = quantile.(eachrow(c_mat), 0.975)
###    end
###    
###    CairoMakie.lines!(ax, ca_range, med_c, linewidth=2, color=(color, line_alpha), label=label)
###    CairoMakie.band!(ax, ca_range, q_025, q_975, linewidth=2, color=(color, ribbon_alpha))
###    
###    return f, ax 
###end
###
###
###function equilibrium_plot(fpms_nt, gs; w=17, h=22.23, save=false, shifman_dir="data/Shifman_2006.xlsx")
###
###
###    f_shifman = XLSX.readxlsx(shifman_dir);
###    X_shifman = f_shifman["Shifman_2006.csv"]["A"][2:end] .* 1e-6;
###    Y_shifman = [i for i in f_shifman["Shifman_2006.csv"]["B"][2:end]];
###    
###    ca_range    = 10 .^ LinRange(-9, log10(7e-5), 100);
###    pract_pop   = gen_dummy_data(ca_range);
###    
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    lims = (nothing, nothing, nothing, nothing)
###    
###    bwell_CaM_tup = get_multiple_eqs(fpms_nt[:blackwell], pract_pop, get_bwell_eqs);
###    bwell_fixed_CaM_tup = get_multiple_eqs(fpms_nt[:blackwell_fixed], pract_pop, get_bwell_eqs);
###   
###    f_eqs, f_eqs_ax11 = plot_CaM_eqs(ca_range, bwell_CaM_tup;       i=1, j=1, title="Scheme 1", 
###        color=gs.our_col,  f=nothing, xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates", figsize=size_units);
###    _, _      = plot_CaM_eqs(ca_range, bwell_fixed_CaM_tup; i=1, j=1, title="Scheme 1", 
###        color=gs.kim_col, f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Kim et al. rates");
###    CairoMakie.scatter!(f_eqs_ax11, X_shifman, Y_shifman, 
###        marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
###    CairoMakie.Legend(f_eqs[1, 2], f_eqs_ax11, "Scheme 1", unique=true, valign=:center, framevisible=false, labelsize=8)
###    CairoMakie.ylims!(f_eqs_ax11, -0.3, 4.3) 
###    
###    
###    bhalla_CaM_tup        = get_multiple_eqs(fpms_nt[:bhalla], pract_pop, get_bhalla_eqs);
###    bhalla_fixed_CaM_tup  = get_multiple_eqs(fpms_nt[:bhalla_fixed], pract_pop, get_bhalla_eqs);
###    
###    f_eqs, f_eqs_ax12 = plot_CaM_eqs(ca_range, bhalla_CaM_tup;       i=2, j=1, title="Scheme 2", 
###        color=gs.our_col,  f=f_eqs, new_axis=true,   xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
###    _, _      = plot_CaM_eqs(ca_range, bhalla_fixed_CaM_tup;         i=2, j=1, title="Scheme 2", 
###        color=gs.bhalla_col, f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Hayer and Bhalla rates");
###    CairoMakie.scatter!(f_eqs_ax12, X_shifman, Y_shifman, 
###        marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
###    CairoMakie.Legend(f_eqs[2, 2], f_eqs_ax12, "Scheme 2", unique=true, valign=:center, framevisible=false, labelsize=8)
###    CairoMakie.ylims!(f_eqs_ax12, -0.3, 4.3)
###    
###    
###    shifman_CaM_tup        = get_multiple_eqs(fpms_nt[:shifman], pract_pop, get_shifman_eqs);
###    shifman_fixed_CaM_tup  = get_multiple_eqs(fpms_nt[:shifman_fixed], pract_pop, get_shifman_eqs);
###    
###    f_eqs, f_eqs_ax13 = plot_CaM_eqs(ca_range, shifman_CaM_tup;       i=3, j=1, title="Scheme 3", 
###        color=gs.our_col,  f=f_eqs, new_axis=true,   xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
###    _, _      = plot_CaM_eqs(ca_range, shifman_fixed_CaM_tup;         i=3, j=1, title="Scheme 3", 
###        color=gs.shifman_col, f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Shifman et al. \n + our rates");
###    CairoMakie.scatter!(f_eqs_ax13, X_shifman, Y_shifman, 
###        marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
###    CairoMakie.Legend(f_eqs[3, 2], f_eqs_ax13, "Scheme 3", unique=true, valign=:center, framevisible=false, labelsize=8)
###    CairoMakie.ylims!(f_eqs_ax12, -0.3, 4.3)
###    
###    
###    pepke_m2_CaM_tup        = get_multiple_eqs(fpms_nt[:pepke_m2], pract_pop, get_pepke_m2_eqs);
###    pepke_m2_fixed_CaM_tup  = get_multiple_eqs(fpms_nt[:pepke_m2_fixed], pract_pop, get_pepke_m2_eqs);
###    
###    f_eqs, f_eqs_ax14 = plot_CaM_eqs(ca_range, pepke_m2_CaM_tup;       i=4, j=1, title="Scheme 4", 
###        color=gs.our_col,  f=f_eqs, new_axis=true,   xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
###    _, _      = plot_CaM_eqs(ca_range, pepke_m2_fixed_CaM_tup;         i=4, j=1, title="Scheme 4", 
###        color=gs.pepke_col, f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Pepke et al. rates");
###    CairoMakie.scatter!(f_eqs_ax14, X_shifman, Y_shifman, 
###        marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
###    CairoMakie.Legend(f_eqs[4, 2], f_eqs_ax14, "Scheme 4", unique=true, valign=:center, framevisible=false, labelsize=8)
###    CairoMakie.ylims!(f_eqs_ax14, -0.3, 4.3) 
###    
###    
###    faas_CaM_tup        = get_multiple_eqs(fpms_nt[:faas], pract_pop, get_faas_eqs);
###    faas_fixed_CaM_tup  = get_multiple_eqs(fpms_nt[:faas_fixed], pract_pop, get_faas_eqs);
###    pepke_m1_fixed_CaM_tup  = get_multiple_eqs(fpms_nt[:pepke_m1_fixed], pract_pop, get_faas_eqs);
###    
###    f_eqs, f_eqs_ax15 = plot_CaM_eqs(ca_range, faas_CaM_tup;       i=5, j=1, title="Scheme 5", 
###        color=gs.our_col,  f=f_eqs, new_axis=true, xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
###    _, _      = plot_CaM_eqs(ca_range, faas_fixed_CaM_tup;         i=5, j=1, title="Scheme 5", 
###        color=gs.faas_col, f=f_eqs,                xlabel=nothing, ylabel=nothing, limits=lims, label="Faas et al. rates");
###    _, _      = plot_CaM_eqs(ca_range, pepke_m1_fixed_CaM_tup;         i=2, j=1, title="Scheme 5", 
###        color=gs.pepke_col, f=f_eqs,               xlabel=nothing, ylabel=nothing, limits=lims, label="Pepke et al. rates");
###    CairoMakie.scatter!(f_eqs_ax15, X_shifman, Y_shifman, 
###        marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
###    CairoMakie.Legend(f_eqs[5, 2], f_eqs_ax15, "Scheme 5", unique=true, valign=:center, framevisible=false, labelsize=8)
###    CairoMakie.ylims!(f_eqs_ax15, -0.3, 4.3) 
###    
###    
###    byrne_CaM_tup        = get_multiple_eqs(fpms_nt[:byrne], pract_pop, get_byrne_eqs);
###    byrne_fixed_CaM_tup  = get_multiple_eqs(fpms_nt[:byrne_fixed], pract_pop, get_byrne_eqs);
###    
###    f_eqs, f_eqs_ax16 = plot_CaM_eqs(ca_range, byrne_CaM_tup;       i=6, j=1, title="Scheme 6", 
###        color=gs.our_col,  f=f_eqs, new_axis=true, ylabel=nothing, limits=lims, label="Our rates");
###    _, _      = plot_CaM_eqs(ca_range, byrne_fixed_CaM_tup;          i=6, j=1, title="Scheme 6", 
###        color=gs.byrne_col, f=f_eqs,                ylabel=nothing, limits=lims, label="Byrne et al. rates");
###    CairoMakie.scatter!(f_eqs_ax16, X_shifman, Y_shifman, 
###        marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
###    CairoMakie.Legend(f_eqs[6, 2], f_eqs_ax16, "Scheme 6", unique=true, valign=:center, framevisible=false, labelsize=8)
###    CairoMakie.ylims!(f_eqs_ax16, -0.3, 4.3) 
###    
###    
###    CairoMakie.Label(f_eqs[:,0], "# Bound Ca²⁺ per CaM", rotation=pi/2)
###    CairoMakie.colgap!(f_eqs.layout, 1, 10)
###    
###    if save
###        CairoMakie.save("plots/Shifman_joint.png", f_eqs, pt_per_unit=1, px_per_unit = px_per_unit, size=size_units)
###    end
###
###    return f_eqs
###end
###
###
###
###function aic_plot(fpms_nt, faas_data, idxs, gs; w=17, h=7.5, save=false)
###
###
###    pops = [faas_data[i] for i in idxs];
###    all_pops = repeat(pops, 13);
###
###
###    all_names = [repeat(["Scheme 1 + our rates"], length(fpms_nt[:blackwell])); 
###                 repeat(["Scheme 1 + Kim et al."], length(fpms_nt[:blackwell_fixed]));
###                 repeat(["Scheme 2 + our rates"], length(fpms_nt[:pepke_m2]));
###                 repeat(["Scheme 2 + Hayer and Bhalla"], length(fpms_nt[:pepke_m2_fixed]));
###                 repeat(["Scheme 3 + our rates"], length(fpms_nt[:pepke_m2]));
###                 repeat(["Scheme 3 + Shifman et al. + our rates"], length(fpms_nt[:pepke_m2_fixed]));
###                 repeat(["Scheme 4 + our rates"], length(fpms_nt[:pepke_m2]));
###                 repeat(["Scheme 4 + Pepke et al."], length(fpms_nt[:pepke_m2_fixed]));
###                 repeat(["Scheme 5 + our rates"], length(fpms_nt[:faas]));
###                 repeat(["Scheme 5 + Faas et al."], length(fpms_nt[:faas_fixed]));
###                 repeat(["Scheme 5 + Pepke et al."], length(fpms_nt[:pepke_m1_fixed]));
###                 repeat(["Scheme 6 + our rates"], length(fpms_nt[:byrne]));
###                 repeat(["Scheme 6 + Byrne et al."], length(fpms_nt[:byrne_fixed]))];
###    
###    aic_pairs = zip(vcat([fpms_nt[key] for key in keys(fpms_nt)]...), all_pops, all_names);
###
###    x, y = aic_boxplots(aic_pairs);
###
###    x_makie = [
###        repeat([1],  length(fpms_nt[:blackwell]));
###        repeat([2],  length(fpms_nt[:blackwell_fixed]));
###        repeat([3],  length(fpms_nt[:bhalla]));
###        repeat([4],  length(fpms_nt[:bhalla_fixed]));
###        repeat([5],  length(fpms_nt[:shifman]));
###        repeat([6],  length(fpms_nt[:shifman_fixed]));
###        repeat([7],  length(fpms_nt[:pepke_m2]));
###        repeat([8],  length(fpms_nt[:pepke_m2_fixed]));
###        repeat([9],  length(fpms_nt[:faas]));
###        repeat([10], length(fpms_nt[:faas_fixed]));
###        repeat([11], length(fpms_nt[:pepke_m1_fixed]));
###        repeat([12], length(fpms_nt[:byrne]));
###        repeat([13], length(fpms_nt[:byrne_fixed]))
###    ]
###
###    y_makie = [i for i in y];
###
###    boxplot_cols = [repeat([gs.our_col],       length(fpms_nt[:blackwell])); 
###                    repeat([gs.kim_col],       length(fpms_nt[:blackwell_fixed]));
###                    repeat([gs.our_col],       length(fpms_nt[:bhalla]));
###                    repeat([gs.bhalla_col],    length(fpms_nt[:bhalla_fixed]));
###                    repeat([gs.our_col],       length(fpms_nt[:shifman]));
###                    repeat([gs.shifman_col],   length(fpms_nt[:shifman_fixed]));
###                    repeat([gs.our_col],       length(fpms_nt[:pepke_m2]));
###                    repeat([gs.pepke_col],     length(fpms_nt[:pepke_m2_fixed]));
###                    repeat([gs.our_col],       length(fpms_nt[:faas]));
###                    repeat([gs.faas_col],      length(fpms_nt[:faas_fixed]));
###                    repeat([gs.pepke_col],     length(fpms_nt[:pepke_m1_fixed]));
###                    repeat([gs.our_col],       length(fpms_nt[:byrne]));
###                    repeat([gs.byrne_col],     length(fpms_nt[:byrne_fixed]))]
###
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    f = CairoMakie.Figure(resolution=size_units, fontsize=gs.fontsize)
###    ax_boxplot = CairoMakie.Axis(f[1,1], 
###        xticks=(unique(x_makie), unique(x)),
###        xlabel="Scheme + parameters",
###        ylabel="log10(AIC)",
###        xticklabelrotation=pi/9,
###        xticklabelpad=5,
###        xlabelpadding=15,
###        xticklabelsize=8,
###        yticklabelsize=8,
###        xgridvisible=false,
###        yscale=log10
###        )
###    CairoMakie.boxplot!(ax_boxplot, x_makie, y_makie, show_notch=false, color=boxplot_cols)
###
###    if save
###        CairoMakie.save("plots/AIC_boxplot.png", f, pt_per_unit=1, 
###            px_per_unit=px_per_unit, size=size_units)   
###    end
###
###    return f
###end
###
###
###function plot_int_lines(f::CairoMakie.Figure, freqs, lines, i, j, col, la, ra, lw, ticks, ticklabs, lims)
###
###    ax = CairoMakie.Axis(f[i, j],
###    xticklabelsize=8, yticklabelsize=8,
###    xticklabelsvisible=ticklabs[1], xticksvisible=ticks[1],
###    yticklabelsvisible=ticklabs[2], yticksvisible=ticks[2],
###    limits=lims[i], ##yscale=log10
###    )
###
###    mat = hcat([ff[i] for ff in lines]...)
###
###    CairoMakie.lines!(ax, freqs, median.(eachrow(mat)), color=(col, la), linewidth=lw)
###    CairoMakie.band!(ax, freqs, quantile.(eachrow(mat), 0.025), quantile.(eachrow(mat), 0.975), color=(col, ra))
###
###    return ax
###end
###
###
###function plot_int_lines(f::CairoMakie.Axis, freqs, lines, i, j, col, la, ra, lw)
###    mat = hcat([ff[i] for ff in lines]...)
###
###    CairoMakie.lines!(f, freqs, median.(eachrow(mat)), color=(col, la), linewidth=lw)
###    CairoMakie.band!(f, freqs, quantile.(eachrow(mat), 0.025), quantile.(eachrow(mat), 0.975), color=(col, ra))
###end
###
###
###function integration_plot(fpms_nt, gs; w=17, h=10, save=false)
###
###    f_range = 2:10:102
###    ca_amp = 0.7e-6
###
###    blackwell_lines       = calc_integration(fpms_nt[:blackwell], f_range, ca_amp, get_bwell_dyns);
###    blackwell_fixed_lines = calc_integration(fpms_nt[:blackwell_fixed], f_range, ca_amp, get_bwell_dyns);
###
###    bhalla_lines       = calc_integration(fpms_nt[:bhalla], f_range, ca_amp, get_bhalla_dyns);
###    bhalla_fixed_lines = calc_integration(fpms_nt[:bhalla_fixed], f_range, ca_amp, get_bhalla_dyns);
###
###    shifman_lines       = calc_integration(fpms_nt[:shifman], f_range, ca_amp, get_shifman_dyns);
###    shifman_fixed_lines = calc_integration(fpms_nt[:shifman_fixed], f_range, ca_amp, get_shifman_dyns);
###
###    pepke_m2_lines       = calc_integration(fpms_nt[:pepke_m2], f_range, ca_amp, get_bwell_CN_dyns);
###    pepke_m2_fixed_lines = calc_integration(fpms_nt[:pepke_m2_fixed], f_range, ca_amp, get_bwell_CN_dyns);
###
###    faas_lines       = calc_integration(fpms_nt[:faas], f_range, ca_amp, get_faas_dyns);
###    faas_fixed_lines = calc_integration(fpms_nt[:faas_fixed], f_range, ca_amp, get_faas_dyns);
###    pepke_m1_fixed_lines = calc_integration(fpms_nt[:pepke_m1_fixed], f_range, ca_amp, get_pepke_m1_dyns);
###
###    byrne_lines       = calc_integration(fpms_nt[:byrne], f_range, ca_amp, get_byrne_dyns);
###    byrne_fixed_lines = calc_integration(fpms_nt[:byrne_fixed], f_range, ca_amp, get_byrne_dyns);
###
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###
###    lw = 1
###    line_alpha = 0.6
###    ribbon_alpha = 0.2
###    int_limits = [(nothing, nothing, -1e-3, 0.016); (nothing, nothing, -1e-4, 0.0015)]
###    ###int_limits = [(nothing, nothing, nothing, nothing); (nothing, nothing, nothing, nothing)]
###
###    f_int = CairoMakie.Figure(resolution=size_units, fontsize=8)
###
###    int_ax_11 = plot_int_lines(f_int, collect(f_range), blackwell_lines, 1, 2, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (false, true), (false, true), int_limits);
###    plot_int_lines(int_ax_11, collect(f_range), blackwell_fixed_lines,   1, 2, gs.kim_col, line_alpha, ribbon_alpha, lw)
###    int_ax_21 = plot_int_lines(f_int, collect(f_range), blackwell_lines, 2, 2, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (true, true), (true, true), int_limits);
###    plot_int_lines(int_ax_21, collect(f_range), blackwell_fixed_lines,   2, 2, gs.kim_col, line_alpha, ribbon_alpha, lw);
###    elem_our = LineElement(color = gs.our_col, linestyle = nothing);
###    elem_kim = LineElement(color = gs.kim_col, linestyle = nothing);
###    CairoMakie.Legend(f_int[0, 2], [elem_our, elem_kim], 
###        ["Our \n rates", "Kim \n et al."], "Scheme 1", framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center, valign=:top)
###
###
###    int_ax_12 = plot_int_lines(f_int, collect(f_range), bhalla_lines, 1, 3, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
###    plot_int_lines(int_ax_12, collect(f_range), bhalla_fixed_lines,   1, 3, gs.bhalla_col, line_alpha, ribbon_alpha, lw);
###    int_ax_22 = plot_int_lines(f_int, collect(f_range), bhalla_lines, 2, 3, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
###    plot_int_lines(int_ax_22, collect(f_range), bhalla_fixed_lines,   2, 3, gs.bhalla_col, line_alpha, ribbon_alpha, lw);
###    elem_bhalla = LineElement(color = gs.bhalla_col, linestyle = nothing);
###    CairoMakie.Legend(f_int[0, 3], [elem_our, elem_bhalla], 
###        ["Our \n rates", "Hayer \n and Bhalla"], "Scheme 2", framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center, valign=:top)
###
###
###    int_ax_13 = plot_int_lines(f_int, collect(f_range), shifman_lines, 1, 4, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
###    plot_int_lines(int_ax_13, collect(f_range), shifman_fixed_lines,   1, 4, gs.shifman_col, line_alpha, ribbon_alpha, lw);
###    int_ax_23 = plot_int_lines(f_int, collect(f_range), shifman_lines, 2, 4, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
###    plot_int_lines(int_ax_23, collect(f_range), shifman_fixed_lines,   2, 4, gs.shifman_col, line_alpha, ribbon_alpha, lw);
###    elem_shifman = LineElement(color = gs.shifman_col, linestyle = nothing);
###    CairoMakie.Legend(f_int[0, 4], [elem_our, elem_shifman], 
###        ["Our \n rates", "Shifman \n et al."], "Scheme 3", framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center, valign=:top)
###    
###    
###    int_ax_14 = plot_int_lines(f_int, collect(f_range), pepke_m2_lines, 1, 5, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
###    plot_int_lines(int_ax_14, collect(f_range), pepke_m2_fixed_lines,   1, 5, gs.pepke_col, line_alpha, ribbon_alpha, lw);
###    int_ax_24 = plot_int_lines(f_int, collect(f_range), pepke_m2_lines, 2, 5, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
###    plot_int_lines(int_ax_24, collect(f_range), pepke_m2_fixed_lines,   2, 5, gs.pepke_col, line_alpha, ribbon_alpha, lw);
###    elem_pepke = LineElement(color = gs.pepke_col, linestyle = nothing);
###    CairoMakie.Legend(f_int[0, 5], [elem_our, elem_pepke], 
###        ["Our \n rates", "Pepke \n et al."], "Scheme 4", framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center, valign=:top)
###    
###    
###    int_ax_15 = plot_int_lines(f_int, collect(f_range), faas_lines,     1, 6, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
###    plot_int_lines(int_ax_15, collect(f_range), faas_fixed_lines,       1, 6, gs.faas_col, line_alpha, ribbon_alpha, lw);
###    plot_int_lines(int_ax_15, collect(f_range), pepke_m1_fixed_lines,   1, 6, gs.pepke_col, line_alpha, ribbon_alpha, lw);
###    int_ax_25 = plot_int_lines(f_int, collect(f_range), faas_lines,     2, 6, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
###    plot_int_lines(int_ax_25, collect(f_range), faas_fixed_lines,       2, 6, gs.faas_col, line_alpha, ribbon_alpha, lw);
###    plot_int_lines(int_ax_25, collect(f_range), pepke_m1_fixed_lines,   2, 6, gs.pepke_col, line_alpha, ribbon_alpha, lw);
###    elem_faas = LineElement(color = gs.faas_col, linestyle = nothing);
###    CairoMakie.Legend(f_int[0, 6], [elem_our, elem_faas, elem_pepke], 
###        ["Our rates", "Faas et al.", "Pepke et al."], "Scheme 5", framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center, valign=:top)
###    
###
###    int_ax_16 = plot_int_lines(f_int, collect(f_range), byrne_lines,     1, 7, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
###    plot_int_lines(int_ax_16, collect(f_range), byrne_fixed_lines,       1, 7, gs.byrne_col, line_alpha, ribbon_alpha, lw);
###    int_ax_26 = plot_int_lines(f_int, collect(f_range), byrne_lines,     2, 7, 
###        gs.our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
###    plot_int_lines(int_ax_26, collect(f_range), byrne_fixed_lines,       2, 7, gs.byrne_col, line_alpha, ribbon_alpha, lw);
###    elem_byrne = LineElement(color = gs.byrne_col, linestyle = nothing);
###    CairoMakie.Legend(f_int[0, 7], [elem_our, elem_byrne], 
###        ["Our \n rates", "Byrne \n et al."], "Scheme 6", framevisible=false, labelsize=8, nbanks=1, 
###        labeljustification=:center, valign=:top)
###    
###    CairoMakie.Label(f_int[1,1], "Partially bound CaM", rotation=pi/2)
###    CairoMakie.Label(f_int[2,1], "Fully bound CaM", rotation=pi/2)
###    CairoMakie.Label(f_int[1:2,0], "AUC", rotation=pi/2)
###    CairoMakie.Label(f_int[3,:], "Frequency (Hz)")
###    
###    for i in 1:2
###        rowsize!(f_int.layout, i, Relative(1/3))
###    end
###    rowsize!(f_int.layout, 0, Relative(1/2.5))
###    for i in 2:7
###        colsize!(f_int.layout, i, Relative(1/6))
###    end
###    
###    CairoMakie.colgap!(f_int.layout, 5)
###    CairoMakie.rowgap!(f_int.layout, 10)
###    
###    if save
###        CairoMakie.save("plots/integration.png", f_int, pt_per_unit=1, px_per_unit = px_per_unit, size=size_units)
###    end
###
###    return f_int
###end
###
###
###function add_heatmap(f, nt, i, j, dd, title, msize)
###
###    __, CC, PP = get_ppcor_matrix(nt)
###
###    PP_no_diag = PP .* .!I(first(size(PP)))
###    sig_mat  = PP_no_diag .< 0.05 .&& PP_no_diag .> 0
###    sig_dots = [(i[1], i[2]) for i in findall(sig_mat)]
###
###    ax = CairoMakie.Axis(f[i,j],
###        xticks=(collect(1:length(dd)), collect(map(k -> dd[k], keys(nt)))),
###        yticks=(collect(1:length(dd)), collect(map(k -> dd[k], keys(nt)))),
###        xticklabelrotation = pi/2,
###        title=title,
###        titlesize=12
###    )
###
###    hm = CairoMakie.heatmap!(ax, collect(1:length(dd)), collect(1:length(dd)), CC, colorrange = (-1, 1))
###    CairoMakie.scatter!(ax, sig_dots, color=:white, marker=:star8, markersize=msize)
###
###    return hm
###end
###
###
###function correlation_plot(fpms_nt, gs; w=17, h=22, save=false)
###
###    coef_bwell = get_coef_lists(fpms_nt[:blackwell], 
###        (:tv_kon_1, :tv_kon_2,
###         :tv_Kd_1,  :tv_Kd_2);
###    );
###    dd_bwell = Dict(:tv_kon_1 => "k₁",   :tv_kon_2 => "k₃",
###                    :tv_Kd_1  => "K_D₁", :tv_Kd_2  => "K_D₂")
###    
###    
###    coef_bhalla = get_coef_lists(fpms_nt[:bhalla], 
###        (:tv_kon_1, :tv_kon_2, :tv_kon_3,
###         :tv_Kd_1,  :tv_Kd_2,  :tv_Kd_3);
###    );
###    dd_bhalla = Dict(:tv_kon_1 => "k₁",   :tv_kon_2 => "k₃",   :tv_kon_3 => "k₅",
###                     :tv_Kd_1  => "K_D₁", :tv_Kd_2  => "K_D₂", :tv_Kd_3  => "K_D₃")
###    
###    
###    coef_shifman = get_coef_lists(fpms_nt[:shifman], 
###        (:tv_kon_1, :tv_kon_2, :tv_kon_3, :tv_kon_4,
###         :tv_Kd_1,  :tv_Kd_2,  :tv_Kd_3,  :tv_Kd_4);
###    );
###    dd_shifman = Dict(:tv_kon_1 => "k₁",   :tv_kon_2 => "k₃",   :tv_kon_3 => "k₅",   :tv_kon_4 => "k₇",
###                      :tv_Kd_1  => "K_D₁", :tv_Kd_2  => "K_D₂", :tv_Kd_3  => "K_D₃", :tv_Kd_4  => "K_D₄")
###    
###    
###    pepke_m2_coef_counternames = zip(
###        (:tv_kon_TC, :tv_kon_RC, :tv_Kd_TC, :tv_Kd_RC),
###        (:tv_kon_TN, :tv_kon_RN, :tv_Kd_TN, :tv_Kd_RN),
###    )
###    pepke_m2_keys_ordered = (:tv_kon_TC, :tv_kon_RC, :tv_kon_TN, :tv_kon_RN,
###                             :tv_Kd_TC,  :tv_Kd_RC,  :tv_Kd_TN,  :tv_Kd_RN)
###    pepke_m2_labs_ordered = ("k₁", "k₃", "k₅", "k₇",
###                             "K_D₁", "K_D₂", "K_D₃", "K_D₄")
###    coef_pepke_m2 = get_coef_lists(fpms_nt[:pepke_m2], 
###        (:tv_kon_RC, :tv_kon_RN, :tv_kon_TC, :tv_kon_TN,
###         :tv_Kd_RC,  :tv_Kd_RN,  :tv_Kd_TC,  :tv_Kd_TN);
###         sort = true,
###         kd_C_name=:tv_Kd_TC,
###         kd_N_name=:tv_Kd_TN,
###         counternames = pepke_m2_coef_counternames
###    );
###    coef_pepke_m2 = NamedTuple{pepke_m2_keys_ordered}(coef_pepke_m2)
###    dd_pepke_m2 = Dict(zip(pepke_m2_keys_ordered, pepke_m2_labs_ordered))
###    
###    
###    faas_coef_counternames = zip(
###        (:tv_kon_TC, :tv_kon_RC, :tv_Kd_TC, :tv_Kd_RC),
###        (:tv_kon_TN, :tv_kon_RN, :tv_Kd_TN, :tv_Kd_RN),
###    )
###    faas_keys_ordered = (:tv_kon_TC, :tv_kon_RC, :tv_kon_TN, :tv_kon_RN,
###                         :tv_Kd_TC,  :tv_Kd_RC,  :tv_Kd_TN,  :tv_Kd_RN)
###    faas_labs_ordered = ("k₁", "k₃", "k₅", "k₇",
###                         "K_D₁", "K_D₂", "K_D₃", "K_D₄")
###    coef_faas = get_coef_lists(fpms_nt[:faas], 
###        (:tv_kon_RC, :tv_kon_RN, :tv_kon_TC, :tv_kon_TN,
###         :tv_Kd_RC,  :tv_Kd_RN,  :tv_Kd_TC,  :tv_Kd_TN);
###         sort = true,
###         kd_C_name=:tv_Kd_TC,
###         kd_N_name=:tv_Kd_TN,
###         counternames = faas_coef_counternames
###    );
###    coef_faas = NamedTuple{faas_keys_ordered}(coef_faas)
###    dd_faas = Dict(zip(faas_keys_ordered, faas_labs_ordered))
###    
###    
###    byrne_coef_counternames = zip(
###        (:tv_k01_N, :tv_k02_N, :tv_k13_N, :tv_k23_N, :tv_K01d_N, :tv_K13d_N, :tv_K02d_N),
###        (:tv_k01_C, :tv_k02_C, :tv_k13_C, :tv_k23_C, :tv_K01d_C, :tv_K13d_C, :tv_K02d_C)
###    )
###    byrne_keys_ordered = (:tv_k01_N,  :tv_k02_N,  :tv_k13_N, :tv_k23_N,
###                          :tv_k01_C,  :tv_k02_C,  :tv_k13_C, :tv_k23_C,
###                          :tv_K01d_N, :tv_K13d_N, :tv_K02d_N, 
###                          :tv_K01d_C, :tv_K13d_C, :tv_K02d_C
###    ) 
###    byrne_labs_ordered = ("kⁿ₀₁", "kⁿ₀₂", "kⁿ₁₃", "kⁿ₂₃",
###                          "kᶜ₀₁", "kᶜ₀₂", "kᶜ₁₃", "kᶜ₂₃",
###                          "K_Dⁿ₀₁", "K_Dⁿ₁₃", "K_Dⁿ₀₂", 
###                          "K_Dᶜ₀₁", "K_Dᶜ₁₃", "K_Dᶜ₀₂"
###    )
###    coef_byrne = get_coef_lists(fpms_nt[:byrne], 
###        (:tv_k01_N,     :tv_k02_N,   :tv_k13_N, :tv_k23_N,
###         :tv_K01d_N,    :tv_K13d_N,  :tv_K02d_N, 
###         :tv_k01_C,     :tv_k02_C,   :tv_k13_C, :tv_k23_C,
###         :tv_K01d_C,    :tv_K13d_C,  :tv_K02d_C);
###         sort = true,
###         kd_C_name=:tv_K01d_C,
###         kd_N_name=:tv_K01d_N,
###         counternames = byrne_coef_counternames
###    );
###    coef_byrne = NamedTuple{byrne_keys_ordered}(coef_byrne)
###    dd_byrne = Dict(:tv_k01_N  => "kⁿ₀₁",   :tv_k02_N  => "kⁿ₀₂",   :tv_k13_N  => "kⁿ₁₃", :tv_k23_N => "kⁿ₂₃",
###                :tv_K01d_N => "K_Dⁿ₀₁", :tv_K13d_N => "K_Dⁿ₁₃", :tv_K02d_N => "K_Dⁿ₀₂", 
###                :tv_k01_C  => "kᶜ₀₁",   :tv_k02_C  => "kᶜ₀₂",   :tv_k13_C  => "kᶜ₁₃", :tv_k23_C => "kᶜ₂₃",
###                :tv_K01d_C => "K_Dᶜ₀₁", :tv_K13d_C => "K_Dᶜ₁₃", :tv_K02d_C => "K_Dᶜ₀₂")
###    
###    
###    
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    f = CairoMakie.Figure(resolution=size_units, fontsize=8)
###
###    add_heatmap(f, coef_bwell,   1, 1       , dd_bwell,    "Scheme 1", 05)
###    add_heatmap(f, coef_bhalla,  1, 2       , dd_bhalla,   "Scheme 2", 05)
###    add_heatmap(f, coef_shifman, 2, 1       , dd_shifman,  "Scheme 3", 05)
###    add_heatmap(f, coef_pepke_m2,2, 2       , dd_pepke_m2, "Scheme 4", 05)
###    add_heatmap(f, coef_faas,   3:4, 1:2    , dd_faas,     "Scheme 5", 05)
###    f_cbar = add_heatmap(f, coef_byrne,  5:6, 1:2    , dd_byrne,    "Scheme 6", 05)
###    CairoMakie.Colorbar(f[:, 3], f_cbar, label="Partial Correlation Coefficient", labelsize=10)
###    
###    CairoMakie.rowgap!(f.layout, 10)
###    CairoMakie.colgap!(f.layout, 10)
###    
###    if save
###        CairoMakie.save("plots/correlations.png", f, pt_per_unit=1, px_per_unit = px_per_unit, size=size_units)
###    end
###
###    return f
###end
###
###
###
###function byrne_c_n_plot(fpms, gs; w=17, h=17, save=false)
###
###
###    c_mat, t_mat, catot = byrne_N_C_lobe_binding_comp(fpms);
###
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###
###
###    linecol = (:black, 0.025)
###    f = CairoMakie.Figure(resolution=size_units, fontsize=8)
###
###    f_ax11 = CairoMakie.Axis(f[1,1], 
###        title = "k₀₁ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₁ᴰ = 33×10⁻⁶ M \n
###    k₀₂ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₂ᴰ = 229.88×10⁻⁶ M",
###        xscale=log10,
###        limits = (nothing, nothing, 0, 1.1),
###        xticks = [0.16, 1, 5, 35],
###        xticksvisible = false,
###        xticklabelsvisible = false
###        )
###    for (x, y) in zip(t_mat.n0, c_mat.n0)  
###        CairoMakie.lines!(f_ax11,
###        x, y, color=linecol )
###    end
###   
###    f_ax12 = CairoMakie.Axis(f[1,2], #ylabel="N1 / CaM_total", 
###        title = "k₀₁ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₁ᴰ = 33×10⁻⁶ M \n
###    k₁₃ = 5.072×10⁵ M⁻¹ms⁻¹ \n K₁₃ᴰ = 3.45×10⁻⁶ M",
###        xscale=log10,
###        limits = (nothing, nothing, 0, 1.1),
###        xticks = [0.16, 1, 5, 35],
###        xticksvisible = false,
###        xticklabelsvisible = false,
###        yticksvisible = false,
###        yticklabelsvisible = false
###        )
###    for (x, y) in zip(t_mat.n1, c_mat.n1)  
###        CairoMakie.lines!(f_ax12,
###        x, y, color=linecol )
###    end
###    
###    f_ax13 = CairoMakie.Axis(f[1,3], #ylabel="N2 / CaM_total", 
###        title = "k₀₂ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₂ᴰ = 229.88×10⁻⁶ M \n
###    k₂₃ = 5×10⁵ M⁻¹ms⁻¹ \n K₂₃ᴰ = 0.5×10⁻⁶ M",
###        xscale=log10,
###        limits = (nothing, nothing, 0, 1.1),
###        xticks = [0.16, 1, 5, 35],
###        xticksvisible = false,
###        xticklabelsvisible = false,
###        yticksvisible = false,
###        yticklabelsvisible = false
###        )
###    for (x, y) in zip(t_mat.n2, c_mat.n2)  
###        CairoMakie.lines!(f_ax13,
###        x, y, color=linecol)
###    end
###    
###    f_ax14 = CairoMakie.Axis(f[1,4], #ylabel="N3 / CaM_total", 
###        title = "k₁₃ = 5.072×10⁵ M⁻¹ms⁻¹ \n K₁₃ᴰ=3.45×10⁻⁶ M \n
###    k₂₃=5×10⁵ M⁻¹ms⁻¹ \n K₂₃ᴰ=0.5×10⁻⁶ M",
###        xscale=log10,
###        limits = (nothing, nothing, 0, 1.1),
###        xticks = [0.16, 1, 5, 35],
###        xticksvisible = false,
###        xticklabelsvisible = false,
###        yticksvisible = false,
###        yticklabelsvisible = false
###        )
###    for (x, y) in zip(t_mat.n3, c_mat.n3)  
###        CairoMakie.lines!(f_ax14,
###        x, y, color=linecol)
###    end
###    
###    
###    
###    f_ax21 = CairoMakie.Axis(f[2,1], #ylabel="C0 / CaM_total", 
###        title = "k₀₁ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₁ᴰ = 18.5×10⁻⁶ M \n
###    k₀₂ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₂ᴰ = 116×10⁻⁶ M",
###        xscale=log10,
###        limits = (nothing, nothing, 0, 1.1),
###        xlabel="log₁₀(t (ms))",
###        xticks = [0.16, 1, 5, 35]
###        )
###    for (x, y) in zip(t_mat.c0, c_mat.c0)  
###        CairoMakie.lines!(f_ax21,
###        x, y, color=linecol )
###    end
###    
###    f_ax22 = CairoMakie.Axis(f[2,2], #ylabel="C1 / CaM_total", 
###        title = "k₀₁ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₁ᴰ = 18.5×10⁻⁶ M \n
###    k₁₃ = 3.71×10³ M⁻¹ms⁻¹ \n K₁₃ᴰ = 0.38×10⁻⁶ M",
###        xscale=log10,
###        limits = (nothing, nothing, 0, 1.1),
###        xlabel="log₁₀(t (ms))",
###        xticks = [0.16, 1, 5, 35],
###        yticksvisible = false,
###        yticklabelsvisible = false
###    )
###    for (x, y) in zip(t_mat.c1, c_mat.c1)  
###        CairoMakie.lines!(f_ax22,
###        x, y, color=linecol)
###    end
###    
###    f_ax23 = CairoMakie.Axis(f[2,3], #ylabel="C2 / CaM_total", 
###        title = "k₀₂ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₂ᴰ = 116×10⁻⁶ M \n
###    k₂₃ = 1.18×10⁵ M⁻¹ms⁻¹ \n K₂₃ᴰ = 0.06×10⁻⁶ M",
###        xscale=log10,
###        limits = (nothing, nothing, 0, 1.1),
###        xticks = [0.16, 1, 5, 35],
###        yticksvisible = false,
###        yticklabelsvisible = false,
###        xlabel="log₁₀(t (ms))")
###    for (x, y) in zip(t_mat.c2, c_mat.c2)  
###        CairoMakie.lines!(f_ax23,
###        x, y, color=linecol)
###    end
###    
###    f_ax24 = CairoMakie.Axis(f[2,4], #ylabel="C3 / CaM_total", 
###        title = "k₁₃ = 3.71×10³ M⁻¹ms⁻¹ \n K₁₃ᴰ = 0.38×10⁻⁶ M \n
###    k₂₃ = 1.18×10⁵ M⁻¹ms⁻¹ \n K₂₃ᴰ = 0.06×10⁻⁶ M",
###        xscale=log10,
###        limits = (nothing, nothing, 0, 1.1),
###        xticks = [0.16, 1, 5, 35],
###        yticksvisible = false,
###        yticklabelsvisible = false,
###        xlabel="log₁₀(t (ms))")
###    for (x, y) in zip(t_mat.c3, c_mat.c3)  
###        CairoMakie.lines!(f_ax24,
###        x, y, color=linecol)
###    end
###    
###    CairoMakie.Label(f[1,0], "N lobe (fraction)", rotation=pi/2)
###    CairoMakie.Label(f[2,0], "C lobe (fraction)", rotation=pi/2)
###    
###    CairoMakie.Label(f[0,1], "CaM + 0Ca")
###    CairoMakie.Label(f[0,2], "CaM + 1Ca \n First site")
###    CairoMakie.Label(f[0,3], "CaM + 1Ca \n Second site")
###    CairoMakie.Label(f[0,4], "CaM + 2Ca")
###    
###    for i in 1:2
###        rowsize!(f.layout, i, Relative(1/2))
###    end
###    for i in 1:4
###        colsize!(f.layout, i, Relative(1/4))
###    end
###    colgap!(f.layout, 1, 5)
###    rowgap!(f.layout, 1, 5)
###
###    if save
###        CairoMakie.save("plots/byrne_lobe_analysis.png", f, pt_per_unit=1, px_per_unit=px_per_unit,
###            size=size_units)
###    end
###
###   return f
###end
###
###
###function pair_plots(all_params, par_names, cs, names, markers, s_widths;
###    lims = (0, 1, 0, 1),
###    ms=10, 
###    fig_w=17 * 28.3465,
###    fig_h=23 * 28.3465,
###    xvis = true,
###    yvis = false,
###    rot_xlab = 0
###    )
###
###    used_names = []    
###
###    f = CairoMakie.Figure(resolution=(fig_w, fig_h), fontsize=8)
###
###    for n_i in keys(all_params[1])
###        push!(used_names, n_i)
###
###        col_ctr=1
###
###        for n_j in setdiff(keys(all_params[1]), used_names)[end:-1:1]
###
###            ax = CairoMakie.Axis(f[length(used_names),col_ctr], 
###                xlabel=par_names[n_j], 
###                ylabel=par_names[n_i], 
###                xscale=log10, yscale=log10,
###                limits = lims
###            )
###
###
###            for i in 1:length(all_params)
###
###                xs = all_params[i][n_j]
###                ys = all_params[i][n_i]
###
###                CairoMakie.scatter!(ax, xs, ys, 
###                    markersize=ms, color=cs[i], label=names[i], 
###                    marker=markers[i], strokewidth=s_widths[i])
###            end
###
###            ax.xtickformat = values -> ["$(log10(value))" for value in values]
###            ax.ytickformat = values -> ["$(log10(value))" for value in values]
###            ax.xticklabelrotation = rot_xlab
###            ax.xlabelsize = 8
###            ax.ylabelsize = 8
###
###            if col_ctr == 1 && (length(used_names) == (length(par_names) - 1))
###                ax.xticklabelsvisible = true
###                ax.xticksvisible = true
###                ax.xlabelvisible = true
###                ax.yticklabelsvisible = true
###                ax.yticksvisible = true
###                ax.ylabelvisible = true
###            elseif (length(used_names) + col_ctr) == length(par_names)
###                ax.xticklabelsvisible = true
###                ax.xticksvisible = true 
###                ax.xlabelvisible = true 
###                ax.yticklabelsvisible = false
###                ax.yticksvisible = false 
###                ax.ylabelvisible = false
###            elseif col_ctr == 1 
###                ax.xticklabelsvisible = false
###                ax.xticksvisible = false
###                ax.xlabelvisible = false
###                ax.yticklabelsvisible = true
###                ax.yticksvisible = true
###                ax.ylabelvisible = true
###            else 
###                ax.xticklabelsvisible = false
###                ax.xticksvisible = false 
###                ax.xlabelvisible = false
###                ax.yticklabelsvisible = false
###                ax.yticksvisible = false
###                ax.ylabelvisible = false
###
###            end
###
###            col_ctr += 1
###
###        end
###    end
###    return f
###end
###
###
###function blackwell_parameter_pairplot(fpms_bwell, fpms_bwell_fixed, gs; w=17, h=17, save=false)
###
###    par_bwell         = get_param_lists(fpms_bwell, (:kon_1, :koff_1, :kon_2, :koff_2));
###    par_bwell_fixed   = get_param_lists(fpms_bwell_fixed, (:kon_1, :koff_1, :kon_2, :koff_2));    
###
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    bwell_lims = (10^-14, 10^14, 10^-14, 10^14)
###    f = pair_plots(
###        [par_bwell, par_bwell_fixed],
###        (; kon_1="log₁₀(k₁)", koff_1="log₁₀(k₂)", 
###           kon_2="log₁₀(k₃)", koff_2="log₁₀(k₄)"),
###        [(gs.our_col, 0.5), (gs.kim_col, 1.0)],
###        ["Our rates", "Kim et al."],
###        [:circle, :x],
###        [0, 0.5];
###        lims=bwell_lims,
###        fig_w = size_units[1],
###        fig_h = size_units[2],
###        yvis=true,
###        ms=6
###    )
###    CairoMakie.Legend(f[0, :], f.content[1], "Scheme 1",
###        unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
###        nbanks=3)
###    for i in 1:3
###        rowsize!(f.layout, i, Relative(1/3.3))
###        colsize!(f.layout, i, Relative(1/2.9))
###    end
###    colgap!(f.layout, 10)
###    rowgap!(f.layout, 1, 5)
###    rowgap!(f.layout, 2, -22.5)
###    rowgap!(f.layout, 3, -22.5)
###    
###
###    if save
###        CairoMakie.save("plots/blackwell_pairplot.png", f, pt_per_unit=1, px_per_unit=px_per_unit,
###            size=size_units)
###    end
###
###    return f
###end
###
###
###function bhalla_parameter_pairplot(fpms_bhalla, fpms_bhalla_fixed, gs; w=17, h=22.23, save=false)
###
###    par_bhalla        = get_param_lists(fpms_bhalla,       (:kon_1, :koff_1, :kon_2, :koff_2, :kon_3, :koff_3));
###    par_bhalla_fixed  = get_param_lists(fpms_bhalla_fixed, (:kon_1, :koff_1, :kon_2, :koff_2, :kon_3, :koff_3));
###
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    bhalla_lims = (10^-11, 10^11, 10^-11, 10^11)
###    
###    f = pair_plots(
###        [par_bhalla, par_bhalla_fixed],
###        (; kon_1="log₁₀(k₁)", koff_1="log₁₀(k₂)", 
###           kon_2="log₁₀(k₃)", koff_2="log₁₀(k₄)",
###           kon_3="log₁₀(k₅)", koff_3="log₁₀(k₆)"),
###        [(gs.our_col, 0.5), (gs.bhalla_col, 1.0)],
###        ["Our rates", "Bhalla and Iyengar."],
###        [:circle, :x],
###        [0, 0.5];
###        lims=bhalla_lims,
###        fig_w = size_units[1],
###        fig_h = size_units[2],
###        yvis=true,
###        ms=6,
###        rot_xlab = pi/2
###    )
###    CairoMakie.Legend(f[0, :], f.content[1], "Scheme 2",
###        unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
###        nbanks=2)
###    for i in 1:5
###        rowsize!(f.layout, i, Relative(1/5.3))
###        colsize!(f.layout, i, Relative(1/4.9))
###    end
###    colgap!(f.layout, 4)
###    rowgap!(f.layout, -32.5)
###    rowgap!(f.layout, 1, 5)
###    
###
###    if save
###        CairoMakie.save("plots/bhalla_pairplot.png", f, pt_per_unit=1, px_per_unit=px_per_unit,
###            size=size_units)
###    end
###
###    return f
###end
###
###
###function shifman_parameter_pairplot(fpms_shifman, fpms_shifman_fixed, gs; w=17, h=22.23, save=false)
###
###    par_shifman        = get_param_lists(fpms_shifman,       (:kon_1, :koff_1, :kon_2, :koff_2, :kon_3, :koff_3, :kon_4, :koff_4));
###    par_shifman_fixed  = get_param_lists(fpms_shifman_fixed, (:kon_1, :koff_1, :kon_2, :koff_2, :kon_3, :koff_3, :kon_4, :koff_4));
###
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    shifman_lims = (10^-11, 10^11, 10^-11, 10^11)
###    
###    f = pair_plots(
###        [par_shifman, par_shifman_fixed],
###        (; kon_1="log₁₀(k₁)", koff_1="log₁₀(k₂)", 
###           kon_2="log₁₀(k₃)", koff_2="log₁₀(k₄)",
###           kon_3="log₁₀(k₅)", koff_3="log₁₀(k₆)",
###           kon_4="log₁₀(k₇)", koff_4="log₁₀(k₈)"),
###        [(gs.our_col, 0.5), (gs.shifman_col, 1.0)],
###        ["Our rates", "Shifman et al."],
###        [:circle, :x],
###        [0, 0.5];
###        lims=shifman_lims,
###        fig_w = size_units[1],
###        fig_h = size_units[2],
###        yvis=true,
###        ms=6,
###        rot_xlab = pi/2
###    )
###    CairoMakie.Legend(f[0, :], f.content[1], "Scheme 3",
###        unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
###        nbanks=2)
###    for i in 1:7
###        rowsize!(f.layout, i, Relative(1/7.3))
###        colsize!(f.layout, i, Relative(1/6.9))
###    end
###    colgap!(f.layout, 4)
###    rowgap!(f.layout, -32.5)
###    rowgap!(f.layout, 1, 5)
###    
###
###    if save
###        CairoMakie.save("plots/shifman_pairplot.png", f, pt_per_unit=1, px_per_unit=px_per_unit,
###            size=size_units)
###    end
###
###    return f
###end
###
###
###function pepke_m2_parameter_pairplot(fpms_pepke_m2, fpms_pepke_m2_fixed, gs; w=17, h=22.23, save=false)
###
###    pepke_m2_c_names = (on = :kon_TC, off = :koff_TC)
###    pepke_m2_n_names = (on = :kon_TN, off = :koff_TN)
###    pepke_m2_counternames = zip(
###        (:kon_TN, :koff_TN, :kon_RN, :koff_RN),
###        (:kon_TC, :koff_TC, :kon_RC, :koff_RC))
###    par_pepke_m2      = get_param_lists(fpms_pepke_m2, (:kon_TN, :koff_TN, 
###                                                      :kon_RN, :koff_RN,
###                                                      :kon_TC, :koff_TC,
###                                                      :kon_RC, :koff_RC);
###                                                        sort = true,
###                                                        sort_C_names = pepke_m2_c_names,
###                                                        sort_N_names = pepke_m2_n_names,
###                                                        counternames = pepke_m2_counternames
###                                                        );
###    par_pepke_m2_fixed= get_param_lists(fpms_pepke_m2_fixed, (:kon_TN, :koff_TN, 
###                                                      :kon_RN, :koff_RN,
###                                                      :kon_TC, :koff_TC,
###                                                      :kon_RC, :koff_RC));
###
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    pepke_lims = (10^-5, 10^10, 10^-5, 10^10)
###    f = pair_plots(
###        [par_pepke_m2, par_pepke_m2_fixed],
###        (; kon_TC="log₁₀(k₊ᵗᶜ)", koff_TC="log₁₀(k₋ᵗᶜ)", 
###           kon_RC="log₁₀(k₊ʳᶜ)", koff_RC="log₁₀(k₋ʳᶜ)", 
###           kon_TN="log₁₀(k₊ᵗⁿ)", koff_TN="log₁₀(k₋ᵗⁿ)", 
###           kon_RN="log₁₀(k₊ʳⁿ)", koff_RN="log₁₀(k₋ʳⁿ)"),
###        [(gs.our_col, 0.5), (gs.pepke_col, 1.0)],
###        ["Our rates", "Pepke et al."],
###        [:circle, :x],
###        [0, 0.5];
###        lims = pepke_lims,
###        fig_w = size_units[1],
###        fig_h = size_units[2],
###        yvis=true,
###        ms = 6,
###        rot_xlab = pi/2
###    )
###    CairoMakie.Legend(f[0, :], f.content[1], "Scheme 4",
###        unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
###        nbanks=3)
###    for i in 1:7
###        rowsize!(f.layout, i, Relative(1/7.6))
###        colsize!(f.layout, i, Relative(1/6.5))
###    end
###    colgap!(f.layout, 10)
###    rowgap!(f.layout, 1, 5)
###    
###    for i in 2:7
###        rowgap!(f.layout, i, -20)
###    end
###    
###    if save
###        CairoMakie.save("plots/pepke_m2_pairplot.png", f, pt_per_unit=1, px_per_unit=px_per_unit,
###            size=size_units)
###    end
###
###    return f
###end
###
###
###function faas_parameter_pairplot(fpms_faas, fpms_faas_fixed, fpms_pepke_m1_fixed, gs; w=17, h=22.23, save=false)
###
###    faas_c_names = (on = :kon_TC, off = :koff_TC)
###    faas_n_names = (on = :kon_TN, off = :koff_TN)
###    faas_counternames = zip(
###        (:kon_TN, :koff_TN, :kon_RN, :koff_RN),
###        (:kon_TC, :koff_TC, :kon_RC, :koff_RC))
###    par_faas          = get_param_lists(fpms_faas, (:kon_TN, :koff_TN, 
###                                                      :kon_RN, :koff_RN,
###                                                      :kon_TC, :koff_TC,
###                                                      :kon_RC, :koff_RC);
###                                                      sort = true,
###                                                      sort_C_names = faas_c_names,
###                                                      sort_N_names = faas_n_names,
###                                                      counternames = faas_counternames
###                                                      );
###    par_faas_fixed    = get_param_lists(fpms_faas_fixed, (:kon_TN, :koff_TN, 
###                                                      :kon_RN, :koff_RN,
###                                                      :kon_TC, :koff_TC,
###                                                      :kon_RC, :koff_RC));
###    par_pepke_m1_fixed= get_param_lists(fpms_pepke_m1_fixed, (:kon_TN, :koff_TN, 
###                                                      :kon_RN, :koff_RN,
###                                                      :kon_TC, :koff_TC,
###                                                      :kon_RC, :koff_RC));
###    
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    faas_lims = (10^-5, 10^10, 10^-5, 10^10)
###    f = pair_plots(
###        [par_faas, par_faas_fixed, par_pepke_m1_fixed],
###        (; kon_TC="log₁₀(k₊ᵗᶜ)", koff_TC="log₁₀(k₋ᵗᶜ)", 
###           kon_RC="log₁₀(k₊ʳᶜ)", koff_RC="log₁₀(k₋ʳᶜ)", 
###           kon_TN="log₁₀(k₊ᵗⁿ)", koff_TN="log₁₀(k₋ᵗⁿ)", 
###           kon_RN="log₁₀(k₊ʳⁿ)", koff_RN="log₁₀(k₋ʳⁿ)"),
###        [(gs.our_col, 0.5), (gs.faas_col, 1.0), (gs.pepke_col, 1.0)],
###        ["Our rates", "Faas et al.", "Pepke et al."],
###        [:circle, :x, :utriangle],
###        [0, 0.5, 0.5];
###        lims = faas_lims,
###        fig_w = size_units[1],
###        fig_h = size_units[2],
###        ms=5,
###        rot_xlab = pi/2
###    )
###    CairoMakie.Legend(f[0, :], f.content[1], "Scheme 5",
###        unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
###        nbanks=3)
###    for i in 1:7
###        rowsize!(f.layout, i, Relative(1/7.6))
###        colsize!(f.layout, i, Relative(1/6.5))
###    end
###    colgap!(f.layout, 7)
###    rowgap!(f.layout, 1, 5)
###    
###    for i in 2:7
###        rowgap!(f.layout, i, -25)
###    end
###
###    if save
###        CairoMakie.save("plots/faas_pairplot.png", f, pt_per_unit=1, px_per_unit=px_per_unit,
###            size=size_units)
###    end
###
###    return f
###end
###
###
###function byrne_parameter_pairplot(fpms_byrne, fpms_byrne_fixed, gs; w=17, h=22.23, save=false)
###
###    byrne_c_names = (on = :k01_C, off = :k10_C)
###    byrne_n_names = (on = :k01_N, off = :k10_N)
###    byrne_counternames = zip(
###        (:k01_N, :k10_N, :k02_N, :k20_N, :k13_N, :k31_N, :k23_N, :k32_N),
###        (:k01_C, :k10_C, :k02_C, :k20_C, :k13_C, :k31_C, :k23_C, :k32_C))
###    par_byrne         = get_param_lists(fpms_byrne, (:k01_N, :k10_N, 
###                                                         :k02_N, :k20_N,
###                                                         :k13_N, :k31_N,
###                                                         :k23_N, :k32_N,
###                                                         :k01_C, :k10_C, 
###                                                         :k02_C, :k20_C,
###                                                         :k13_C, :k31_C,
###                                                         :k23_C, :k32_C);
###                                                         sort = true,
###                                                         sort_C_names = byrne_c_names,
###                                                         sort_N_names = byrne_n_names,
###                                                         counternames = byrne_counternames
###                                                         );                                                 
###    par_byrne_fixed   = get_param_lists(fpms_byrne_fixed, (:k01_N, :k10_N, 
###                                                         :k02_N, :k20_N,
###                                                         :k13_N, :k31_N,
###                                                         :k23_N, :k32_N,
###                                                         :k01_C, :k10_C, 
###                                                         :k02_C, :k20_C,
###                                                         :k13_C, :k31_C,
###                                                         :k23_C, :k32_C));
###    
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###    
###    byrne_lims = (10^-6, 10^11, 10^-6, 10^11)
###    f = pair_plots(
###        [par_byrne, par_byrne_fixed],
###        (; k01_N="k₀₁ⁿ", k10_N="k₁₀ⁿ", 
###           k02_N="k₀₂ⁿ", k20_N="k₂₀ⁿ",
###           k13_N="k₁₃ⁿ", k31_N="k₃₁ⁿ",
###           k23_N="k₂₃ⁿ", k32_N="k₃₂ⁿ",
###           k01_C="k₀₁ᶜ", k10_C="k₁₀ᶜ", 
###           k02_C="k₀₂ᶜ", k20_C="k₂₀ᶜ",
###           k13_C="k₁₃ᶜ", k31_C="k₃₁ᶜ",
###           k23_C="k₂₃ᶜ", k32_C="k₃₂ᶜ",
###        ),
###        [(gs.our_col, 0.5), (gs.byrne_col, 1.0)],
###        ["Our rates", "Byrne et al."],
###        [:circle, :x],
###        [0, 0.5];
###        lims = byrne_lims,
###        fig_w = size_units[1],
###        fig_h = size_units[2],
###        ms=3,
###        rot_xlab = pi/2
###    )
###    CairoMakie.Legend(f[0, :], f.content[1], "Scheme 6",
###        unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
###        nbanks=3)
###    
###    for i in 1:15
###        rowsize!(f.layout, i, Relative(1/15.5))
###        colsize!(f.layout, i, Relative(1/14))
###    end
###    
###    colgap!(f.layout, 5)
###    rowgap!(f.layout, 1, 7.5)
###    
###    for i in 2:15
###        rowgap!(f.layout, i, -32.5)
###    end
###    
###    if save
###        CairoMakie.save("plots/byrne_pairplot.png", f, pt_per_unit=1, px_per_unit=px_per_unit,
###            size=size_units)
###    end
###
###    return f
###end
###


function data_plots_fens(faas_data, shifman_dir, gs; w=17, h=7, save=false)

    size_cm     = (w, h)
    size_px     =  size_cm .* gs.dpc 
    size_units  =  size_cm .* gs.pt_per_cm
    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]

    lw  = 1
    ylim= [0, 22.5]

    groups = [faas_data[1:13], 
        faas_data[14:24], 
        faas_data[25:35], 
        faas_data[36:51], 
        faas_data[52:65], 
        faas_data[66:80], 
        faas_data[81:94]
    ]

    f = CairoMakie.Figure(resolution=size_units, fontsize=gs.fontsize)

    ax1 = CairoMakie.Axis(f[1,1], 
    )
    hidedecorations!(ax1)
    hidespines!(ax1)

    yticks = [1; 5; 10; 15; 20]

    ax2 = CairoMakie.Axis(f[1,2], 
        title="Group A \n DMn=5.56mM OGB-5N=50μM \n Ca²⁺=1.88μM CaM=123μM", 
        limits = (-3, 43, ylim[1], ylim[2]), 
        yticks = yticks,
        yticklabelsize=gs.fontsize, 
        xticklabelsize=gs.fontsize,
        xticksvisible=false,
        xticklabelsvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false,
    )
    for i in groups[1]
        CairoMakie.lines!(ax2, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    ax3 = CairoMakie.Axis(f[1,3], 
        title="Group B \n DMn=3.64mM OGB-5N=100μM \n Ca²⁺=1.10μM CaM=143μM", 
        limits = (-3, 43, ylim[1], ylim[2]), 
        yticks = yticks,
        yticklabelsize=gs.fontsize, 
        xticklabelsize=gs.fontsize,
        xticksvisible=false,
        xticklabelsvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false
    )
    for i in groups[2]
        CairoMakie.lines!(ax3, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    ax4 = CairoMakie.Axis(f[1,4], 
        title="Group C \n DMn=3.64mM OGB-5N=100μM \n Ca²⁺=0.73μM CaM=72μM", 
        limits = (-3, 43, ylim[1], ylim[2]), 
        yticks = yticks,
        yticklabelsize=gs.fontsize, 
        xticklabelsize=gs.fontsize,
        xticksvisible=false,
        xticklabelsvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false
    )
    for i in groups[3]
        CairoMakie.lines!(ax4, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    ax5 = CairoMakie.Axis(f[2,1], 
        title="Group D \n DMn=3.64mM OGB-5N=100μM \n Ca²⁺=0.255μM CaM=187μM", 
        ylabel="ΔF/F₀",
        limits = (-3, 43, ylim[1], ylim[2]), 
        yticks = yticks,
        yticklabelsize=gs.fontsize, 
        xticklabelsize=gs.fontsize
    )
    for i in groups[4]
        CairoMakie.lines!(ax5, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    ax6 = CairoMakie.Axis(f[2,2], 
        title="Group E \n DMn=3.64mM OGB-5N=100μM \n Ca²⁺=0.41μM CaM=140μM", 
        limits = (-3, 43, ylim[1], ylim[2]), 
        yticks = yticks,
        yticklabelsize=gs.fontsize, 
        xticklabelsize=gs.fontsize,
        yticksvisible=false,
        yticklabelsvisible=false
    )
    for i in groups[5]
        CairoMakie.lines!(ax6, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    ax7 = CairoMakie.Axis(f[2,3], 
        title="Group F \n DMn=3.64mM OGB-5N=100μM \n Ca²⁺=0.399μM CaM=94μM", 
        limits = (-3, 43, ylim[1], ylim[2]), 
        yticks = yticks,
        yticklabelsize=gs.fontsize, 
        xticklabelsize=gs.fontsize,
        yticksvisible=false,
        yticklabelsvisible=false
    )
    for i in groups[6]
        CairoMakie.lines!(ax7, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    ax8 = CairoMakie.Axis(f[2,4], 
        title="Group G \n DMn=3.64mM OGB-5N=100μM \n Ca²⁺=0.394μM CaM=47μM", 
        limits = (-3, 43, ylim[1], ylim[2]), 
        yticks = yticks,
        yticklabelsize=gs.fontsize, 
        xticklabelsize=gs.fontsize,
        yticksvisible=false,
        yticklabelsvisible=false
    )
    for i in groups[7]
        CairoMakie.lines!(ax8, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    f_shifman = XLSX.readxlsx(shifman_dir);
    X_shifman = f_shifman["Shifman_2006.csv"]["A"][2:end] .* 1e-6;
    Y_shifman = [i for i in f_shifman["Shifman_2006.csv"]["B"][2:end]];

    ax9 = CairoMakie.Axis(f[1:2,5:6], 
        ylabel="# Bound Ca²⁺ per CaM", 
        xlabel="Free Ca²⁺ (M)", 
        yticklabelsize=gs.fontsize, 
        xticklabelsize=gs.fontsize,
    )
    CairoMakie.scatter!(ax9, X_shifman, Y_shifman, color=(:black, 0.4), marker=:utriangle)

    Label(f[0,5:6], "Shifman et al. data", font=:bold, fontsize=gs.fontsize*2)
    Label(f[0,1:4], "Faas et al. data", font=:bold, fontsize=gs.fontsize*2)

    Label(f[3,1:4], "Time(ms)")

    rowsize!(f.layout, 1, Relative(1/2.0))
    rowsize!(f.layout, 2, Relative(1/2.0))

    for i in 1:6
        colsize!(f.layout, i, Relative(1/6))
    end

    rowgap!(f.layout, 10)
    colgap!(f.layout, 10)

    if save
        CairoMakie.save("plots_fens/data_plot.png", f, pt_per_unit=1, px_per_unit=px_per_unit,
        size=size_units)
    end

    return f
end


###function shifman_data_plot(shifman_dir, gs; w=17, h=7, save=false)
###
###    f_shifman = XLSX.readxlsx(shifman_dir);
###    X_shifman = f_shifman["Shifman_2006.csv"]["A"][2:end] .* 1e-6;
###    Y_shifman = [i for i in f_shifman["Shifman_2006.csv"]["B"][2:end]];
###
###    size_cm     = (w, h)
###    size_px     =  size_cm .* gs.dpc 
###    size_units  =  size_cm .* gs.pt_per_cm
###    px_per_unit = (size_cm .* gs.dpc ./ (size_cm .* gs.pt_per_cm))[1]
###
###    f = CairoMakie.Figure(resolution=size_units, fontsize=gs.fontsize)
###
###    ax = CairoMakie.Axis(f[1,1], 
###        ylabel="# Bound Ca²⁺ per CaM", 
###        xlabel="Free Ca²⁺ (M)", 
###        yticklabelsize=gs.fontsize, 
###        xticklabelsize=gs.fontsize,
###    )
###    CairoMakie.scatter!(ax, X_shifman, Y_shifman, color=(:black, 0.3), marker=:utriangle)
###
###    if save
###        CairoMakie.save("plots/data_shifman.png", f, pt_per_unit=1, px_per_unit=px_per_unit,
###            size=size_units)
###    end
###
###    return f
###end