function gen_dummy_data(Ca_range)

    pop = []

    for val in Ca_range
        sub = Subject(;
            id = string(val),
            covariates = (cov_f_frac    = 0, 
                            cov_τ_f     = 1,
                            cov_τ_s     = 1,
                            cov_kon_DMn = 0,
                            cov_koff_DMn= 0,
                            cov_Kd_DMn  = 0,
                            cov_koff_PP = 0,
                            cov_DMn_t   = 0,
                            cov_CaM_t   = 5e-6,
                            cov_OGB5_t  = 0,
                            cov_Fluo4FF_t = 5e-6,
                            cov_Ca_free = val,
                            PCD         = 0,
                            batch_id    = 1,
                            cov_τ_decay = 0,
                            cov_Ca_i    = 0,
                            CaM_equil   = true
                            
            ), 
            observations = (F_F0 = [1.0;],),
            time = [1.0;]
        )
        push!(pop, sub)
    end

    return [pop...]
end


function get_bwell_eqs(fpm, pop)
    no_bound    = []
    part_bound  = []
    full_bound  = []

    for p in pop
        pre = fpm.model.pre(coef(fpm), (η = [0.0; 0.0; 0.0],), p).ret

        push!(no_bound,     pre.CaM₀_all[1] / pre.CaM_t * 0)
        push!(part_bound,   pre.CaM₀_all[2] / pre.CaM_t * 2)
        push!(full_bound,   pre.CaM₀_all[3] / pre.CaM_t * 4)
    end

    return no_bound, part_bound, full_bound
end


function get_pepke_m2_eqs(fpm, pop)
    no_bound    = []
    part_bound  = []
    full_bound  = []

    for p in pop
        pre = fpm.model.pre(coef(fpm), (η = [0.0; 0.0; 0.0],), p).ret

        push!(no_bound,     pre.c0[1] / pre.CaM_t * 0)
        push!(part_bound,   (pre.c0[2] / pre.CaM_t + pre.c0[2] / pre.CaM_t))
        push!(full_bound,   pre.c0[4] / pre.CaM_t * 4)
    end

    return no_bound, part_bound, full_bound
end


function get_byrne_eqs(fpm, pop)
    no_bound    = []
    part_bound  = []
    full_bound  = []

    for p in pop
        pre = fpm.model.pre(coef(fpm), (η = [0.0; 0.0; 0.0],), p).ret

        push!(no_bound,     0)
        b1 = (pre.N₀[2] + pre.C₀[2]) / 2pre.CaM_t * 1
        b2 = (pre.N₀[3] + pre.C₀[3]) / 2pre.CaM_t * 2
        push!(part_bound,   b1+b2)
        push!(full_bound,   (pre.N₀[4] + pre.C₀[4]) / 2pre.CaM_t * 4)
    end

    return no_bound, part_bound, full_bound
end


function get_faas_eqs(fpm, pop)
    no_bound    = []
    part_bound  = []
    full_bound  = []

    Ca, CaM_T, k_off_C1, k_on_C1, k_off_N1, k_on_N1, k_off_C2, k_on_C2, k_off_N2, k_on_N2 =
        @variables Ca, CaM_T, k_off_C1, k_on_C1, k_off_N1, k_on_N1, k_off_C2, k_on_C2, k_off_N2, k_on_N2;

    A, _, b, _ = get_Pepke_eqs()

    for p in pop

        pre = fpm.model.pre(coef(fpm), (η = [0.0; 0.0; 0.0],), p).ret

        dd = Dict(Ca => pre.Ca_free, 
            CaM_T    => pre.CaM_t,
            k_off_C1 => pre.koff_TC,
            k_on_C1  => pre.kon_TC,
            k_off_C2 => pre.koff_RC,
            k_on_C2  => pre.kon_RC,
            k_off_N1 => pre.koff_TN,
            k_on_N1  => pre.kon_TN,
            k_off_N2 => pre.koff_RN,
            k_on_N2  => pre.kon_RN,
        )

        sub_A = Symbolics.value.(substitute.(A, (dd,)))
        sub_b = Symbolics.value.(substitute.(b, (dd,)))
        sol = sub_A\sub_b

        CaM₀  = sol[1] / pre.CaM_t * 0 
        CaMp₀ = sol[2] / pre.CaM_t * 1 + 
                sol[3] / pre.CaM_t * 1 +
                sol[4] / pre.CaM_t * 2 +
                sol[5] / pre.CaM_t * 2 +
                sol[6] / pre.CaM_t * 2 +
                sol[7] / pre.CaM_t * 3 +
                sol[8] / pre.CaM_t * 3 
        CaMf₀ = sol[9] / pre.CaM_t * 4 

        push!(no_bound,     CaM₀)
        push!(part_bound,   CaMp₀)
        push!(full_bound,   CaMf₀)
    end

    return no_bound, part_bound, full_bound
end


function get_multiple_eqs(fpms, pop, single_func)

    no_bound    = []
    part_bound  = []
    full_bound  = []

    for fpm in fpms
        n, p, f = single_func(fpm, pop)
        push!(no_bound, n)
        push!(part_bound, p)
        push!(full_bound, f)
    end

    return no_bound, part_bound, full_bound
end



function plot_CaM_eqs(ca_range, CaM_tup; i=1, j=1, 
    ###xlabel = L"$\textrm{Ca}^{2+}_{\textrm{Free}}$",
    xlabel = "Free Ca²⁺ (M)",
    ylabel = L"Bound $\textrm{Ca}^{2+}$ per CaM",
    title="Scheme 1", color=:red, label="Our rates", 
    f=nothing, new_axis=false, figsize=(800, 600), limits=(0, 1, 0, 1))

if f === nothing
    f = CairoMakie.Figure(size=figsize, fontsize=8)
    ax = CairoMakie.Axis(f[i, j], 
        xticklabelsize=8, yticklabelsize=8,
        xscale=log10
        )
elseif new_axis
    ax = CairoMakie.Axis(f[i, j], limits=limits, 
        xticklabelsize=8, yticklabelsize=8,
        xscale=log10
        )
else
    ax = f.content[end]
end

if ylabel !== nothing
    ax.ylabel = ylabel
###        ax.ylabelsize = 8
end
if xlabel !== nothing
    ax.xlabel = xlabel
###        ax.xlabelsize = 8
end

###    Plots.plot!(ca_range, CaM0_m, line=:dash,     lw=3, color=color, label=label*" CaM0" )
###    Plots.plot!(ca_range, CaMp_m, line=:dashdot,  lw=3, color=color, label=label*" CaM partial" )
###    Plots.plot!(ca_range, CaMf_m, line=:solid,    lw=3, color=color, label=label*" CaM full" )

for (c_n, c_4) in zip(CaM_tup[2], CaM_tup[3])
    CairoMakie.lines!(ax, ca_range, c_n + c_4, linewidth=2, color=color, label=label)
end

return f, ax 
end



function gen_dummy_data_total(Ca_range_total)

    pop = []

    for val in Ca_range_total
        sub = Subject(;
            id = string(val),
            covariates = (cov_f_frac    = 0, 
                            cov_τ_f     = 1,
                            cov_τ_s     = 1,
                            cov_kon_DMn = 0,
                            cov_koff_DMn= 0,
                            cov_Kd_DMn  = 0,
                            cov_koff_PP = 0,
                            cov_DMn_t   = 0,
                            cov_CaM_t   = 40e-6,
                            cov_OGB5_t  = 0,
                            cov_Ca_free = val,
                            PCD         = 0,
                            batch_id    = 1,
                            cov_τ_decay = 0,
                            cov_Ca_i    = 0,
                            CaM_equil   = false
            ), 
            observations = (F_F0 = repeat([missing;], 1000),),
            time = collect(LinRange(0, 1000, 1000))
        )
        push!(pop, sub)
    end

    return [pop...]
end



function get_bwell_dynamic_eqs(fpm, pop)

    CaM0 = []
    CaMp = []
    CaMf = []

    for p in pop

        obs = simobs(fpm.model, p, coef(fpm); diffeq_options=(;alg=RadauIIA5(), abstol=1e-14))

        push!(CaM0, abs.(obs.dynamics.CaM[end]))
        push!(CaMp, abs.(obs.dynamics.CaM2Ca[end]))
        push!(CaMf, abs.(obs.dynamics.CaM4Ca[end]))
    end

    return CaM0, CaMp, CaMf
end


function plot_CaM_dyn_eqs(ca_range, CaM_tup; 
        title="Scheme 1", color="red", label="Ours", 
        cam_conc=40e-6, new_plot=true, prev_plot=nothing,
        legloc=:best)

    if new_plot
        p = Plots.plot(title=title, xlabel="Log10(Free Calcium)", ylabel="Fraction at equilibrium", labelfontsize=12)
    else
        p = Plots.plot(prev_plot)
    end

    Plots.plot!(ca_range, CaM_tup[1][1], line=:dash,     lw=3, color=color, label=label*" CaM0" )
    for ln in CaM_tup[1][2:end]
        Plots.plot!(ca_range, ln, line=:dash,     lw=3, color=color, label=false)
    end

    Plots.plot!(ca_range, CaM_tup[2][1], line=:dashdot,     lw=3, color=color, label=label*" CaM partial" )
    for ln in CaM_tup[2][2:end]
        Plots.plot!(ca_range, ln, line=:dashdot,     lw=3, color=color, label=false)
    end

    Plots.plot!(ca_range, CaM_tup[3][1], line=:solid,     lw=3, color=color, label=label*" CaM full" )
    for ln in CaM_tup[3][2:end]
        Plots.plot!(ca_range, ln, line=:solid,     lw=3, color=color, label=false)
    end

    if new_plot 
        Plots.plot!(xscale=:log10, legend=legloc)
    end
    return p
end



function get_param_lists(fpms, symbols; 
    sort=false, 
    sort_C_names=nothing,
    sort_N_names=nothing,
    counternames=nothing)

    params = []

    for s in symbols
        p_s = Float64[]

        for f in fpms
            pre = f.model.pre(coef(f), (η = [0.0; 0.0; 0.0],), f.data[1]).ret
            push!(p_s, pre[s])
        end

        if length(unique(p_s)) == 1
            p_s = p_s[1:1]
        end

        push!(params, p_s)
    end

    nt = NamedTuple{symbols}(params)

    if sort
        kd_c = nt[sort_C_names.off] ./ nt[sort_C_names.on]
        kd_n = nt[sort_N_names.off] ./ nt[sort_N_names.on]
        
        repl_idx = kd_c .> kd_n

        d = Dict()

        for (n_i, n_j) in counternames
            temp_n_i = nt[n_i]
            temp_n_j = nt[n_j]

            i_to_j = temp_n_i[repl_idx]
            j_to_i = temp_n_j[repl_idx]

            temp_n_i[repl_idx] .= j_to_i
            temp_n_j[repl_idx] .= i_to_j

            d[n_i] = temp_n_i
            d[n_j] = temp_n_j
        end

        ret_nt = NamedTuple(d)
    else
        ret_nt = nt
    end

    return ret_nt
end

function pair_plots(all_params, par_names, cs, names;
    lims = (0, 1, 0, 1),
    ms=10, 
    fig_w=17 * 28.3465,
    fig_h=23 * 28.3465,
    xvis = true,
    yvis = false,
    rot_xlab = 0
    )

    used_names = []    

    f = CairoMakie.Figure(size=(fig_w, fig_h), fontsize=8)

    for n_i in keys(all_params[1])
        push!(used_names, n_i)

        col_ctr=1

        for n_j in setdiff(keys(all_params[1]), used_names)[end:-1:1]

            ax = CairoMakie.Axis(f[length(used_names),col_ctr], 
                xlabel=par_names[n_j], 
                ylabel=par_names[n_i], 
                xscale=log10, yscale=log10,
                limits = lims, #yticklabelsvisible=ylabs, yticksvisible=yticks,
                #xticks=([lims[1], lims[2]], [string(log10(lims[1])), string(log10(lims[2]))]), 
                #yticks=([lims[3], lims[4]], [string(log10(lims[3])), string(log10(lims[4]))])
            )


            for i in 1:length(all_params)

                xs = all_params[i][n_j]
                ys = all_params[i][n_i]

                CairoMakie.scatter!(ax, xs, ys, 
                    markersize=ms, color=cs[i], label=names[i])
            end

            ax.xtickformat = values -> ["$(log10(value))" for value in values]
            ax.ytickformat = values -> ["$(log10(value))" for value in values]
            ax.xticklabelrotation = rot_xlab
            ax.xlabelsize = 8
            ax.ylabelsize = 8

            if col_ctr == 1 && (length(used_names) == (length(par_names) - 1))
                ax.xticklabelsvisible = true
                ax.xticksvisible = true
                ax.xlabelvisible = true
                ax.yticklabelsvisible = true
                ax.yticksvisible = true
                ax.ylabelvisible = true
            elseif (length(used_names) + col_ctr) == length(par_names)
                ax.xticklabelsvisible = true
                ax.xticksvisible = true 
                ax.xlabelvisible = true 
                ax.yticklabelsvisible = false
                ax.yticksvisible = false 
                ax.ylabelvisible = false
            elseif col_ctr == 1 
                ax.xticklabelsvisible = false
                ax.xticksvisible = false
                ax.xlabelvisible = false
                ax.yticklabelsvisible = true
                ax.yticksvisible = true
                ax.ylabelvisible = true
            else 
                ax.xticklabelsvisible = false
                ax.xticksvisible = false 
                ax.xlabelvisible = false
                ax.yticklabelsvisible = false
                ax.yticksvisible = false
                ax.ylabelvisible = false

            end

            col_ctr += 1

        end
    end
    return f
end