function gen_dyn_data()

    sub = Subject(;
        id = "1",
        events = DosageRegimen([DosageRegimen(2e-6, time=i) for i in 20:20:1000]...),
        covariates = (cov_f_frac    = 0, 
                        cov_τ_f     = 1,
                        cov_τ_s     = 1,
                        cov_kon_DMn = 0,
                        cov_koff_DMn= 0,
                        cov_Kd_DMn  = 0,
                        cov_koff_PP = 0,
                        cov_DMn_t   = 0,
                        cov_CaM_t   = 20e-6,
                        cov_OGB5_t  = 0,
                        cov_Fluo4FF_t  = 0,
                        cov_Ca_free = 100e-9,
                        PCD         = 0,
                        batch_id    = 1,
                        cov_τ_decay = 0.0833,
                        cov_Ca_i    = 100e-9,
                        cov_kon_D   = 1.0,
                        cov_koff_D  = 1.0,
                        cov_Kd_Fluo4FF  = 23e-6,
                        cov_kon_Fluo4FF = 1.0,
                        CaM_equil   = true

        ), 
        observations = (F_F0 = repeat([missing], 10001),),
        time = collect(LinRange(0, 1500, 10001))
    )

    return [sub;]
end

function gen_f_data(f, amp)

    sub = Subject(;
        id = "1",
        ###events = DosageRegimen([DosageRegimen(35e-6/f, time=i) for i in LinRange(20, 1020, f)]...),
        events = DosageRegimen([DosageRegimen(amp, time=i) for i in LinRange(20, 1020, f)]...),
        covariates = (cov_f_frac    = 0, 
                        cov_τ_f     = 1,
                        cov_τ_s     = 1,
                        cov_kon_DMn = 0,
                        cov_koff_DMn= 0,
                        cov_Kd_DMn  = 0,
                        cov_koff_PP = 0,
                        cov_DMn_t   = 0,
                        cov_CaM_t   = 20e-6,
                        cov_OGB5_t  = 0,
                        cov_Fluo4FF_t  = 0,
                        cov_Ca_free = 100e-9,
                        PCD         = 0,
                        batch_id    = 1,
                        cov_τ_decay = 0.0833,
                        cov_Ca_i    = 100e-9,
                        cov_kon_D   = 1.0,
                        cov_koff_D  = 1.0,
                        cov_Kd_Fluo4FF  = 23e-6,
                        cov_kon_Fluo4FF = 1.0,
                        CaM_equil   = true

        ), 
        observations = (F_F0 = repeat([missing], 10001),),
        time = collect(LinRange(0, 1500, 10001))
    )

    return [sub;]
end

function gen_dyn_data_2Hz()

    sub = Subject(;
        id = "1",
        events = DosageRegimen([DosageRegimen(12e-6, time=i) for i in 100:500:100000]...),
        covariates = (cov_f_frac    = 0, 
                        cov_τ_f     = 1,
                        cov_τ_s     = 1,
                        cov_kon_DMn = 0,
                        cov_koff_DMn= 0,
                        cov_Kd_DMn  = 0,
                        cov_koff_PP = 0,
                        cov_DMn_t   = 0,
                        cov_CaM_t   = 20e-6,
                        cov_OGB5_t  = 0,
                        cov_Fluo4FF_t  = 0,
                        cov_Ca_free = 100e-9,
                        PCD         = 0,
                        batch_id    = 1,
                        cov_τ_decay = 0.0833,
                        cov_Ca_i    = 100e-9,
                        cov_kon_D   = 1.0,
                        cov_koff_D  = 1.0,
                        cov_Kd_Fluo4FF  = 23e-6,
                        cov_kon_Fluo4FF = 1.0,
                        CaM_equil   = true

        ), 
        observations = (F_F0 = repeat([missing], 10001),),
        time = collect(LinRange(0, 10000, 10001))
    )

    return [sub;]
end


function gen_dyn_eq_data(Ca_range)

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
                            cov_CaM_t   = 1.5e-5,
                            cov_OGB5_t  = 0,
                            cov_Fluo4FF_t = 0,
                            cov_Ca_free = val,
                            PCD         = 0,
                            batch_id    = 1,
                            cov_τ_decay = 0,
                            cov_Ca_i    = 0,
                            CaM_equil   = false 
                            
            ), 
            observations = (F_F0 = repeat([missing], 2),),
            time = collect(LinRange(0, 1e5, 2))
        )
        push!(pop, sub)
    end

    return [pop...]
end



function get_bwell_dyns(fpm, pop)

    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]

    return abs.(obs.dynamics.CaM), abs.(obs.dynamics.CaM2Ca), abs.(obs.dynamics.CaM4Ca)
end


function get_bhalla_dyns(fpm, pop)

    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]

    CaM0 = abs.(obs.dynamics.CaM)
    CaMp = abs.(obs.dynamics.CaM2Ca) + abs.(obs.dynamics.CaM3Ca)
    CaMf = abs.(obs.dynamics.CaM4Ca)

    return CaM0, CaMp, CaMf
end


function get_shifman_dyns(fpm, pop)

    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]

    CaM0 = abs.(obs.dynamics.CaM)
    CaMp = abs.(obs.dynamics.CaM1Ca) + abs.(obs.dynamics.CaM2Ca) + abs.(obs.dynamics.CaM3Ca)
    CaMf = abs.(obs.dynamics.CaM4Ca)

    return CaM0, CaMp, CaMf
end


function get_bwell_CN_dyns(fpm, pop)

    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]

    CaM0 = abs.(obs.dynamics.CaM0)
    CaMp = abs.(obs.dynamics.CaM2C) + abs.(obs.dynamics.CaM2N)
    CaMf = abs.(obs.dynamics.CaM4)

    return CaM0, CaMp, CaMf
end


function get_faas_dyns(fpm, pop)

    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]

    CaM0 = (abs.(obs.dynamics.NtNt) + abs.(obs.dynamics.CtCt)) / 2
    CaMp = (abs.(obs.dynamics.NtNr) + abs.(obs.dynamics.CtCr)) / 2
    CaMf = (abs.(obs.dynamics.NrNr) + abs.(obs.dynamics.CrCr)) / 2

    return CaM0, CaMp, CaMf
end


function get_pepke_m1_dyns(fpm, pop)

    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]

    CaM0 = (abs.(obs.dynamics.N0) + abs.(obs.dynamics.C0)) / 2
    CaMp = (abs.(obs.dynamics.N1) + abs.(obs.dynamics.C1)) / 2
    CaMf = (abs.(obs.dynamics.N2) + abs.(obs.dynamics.C2)) / 2

    return CaM0, CaMp, CaMf
end


function get_byrne_dyns(fpm, pop)

    obs = simobs(fpm.model, pop, coef(fpm); diffeq_options=(;alg=Rodas5P(), abstol=1e-16, reltol=1e-16))[1]
    d = obs.dynamics

    CaM0 = (abs.(d.N0) + abs.(d.C0)) / 2
    CaMp = (abs.(d.N1) + abs.(d.C1) + abs.(d.N2) + abs.(d.C2)) / 4
    CaMf = (abs.(d.N3) + abs.(d.C3)) / 2

    return CaM0, CaMp, CaMf
end


function get_multiple_dyns(fpms, pop, getter_func)

    no_cam = []
    p_cam = []
    f_cam = []

    for fpm in fpms
        n_i, p_i, f_i = getter_func(fpm, pop)

        push!(no_cam, n_i)
        push!(p_cam, p_i) 
        push!(f_cam, f_i)
    end

    return no_cam, p_cam, f_cam
end



function plot_CaM_dyns(t, CaM_tup; 
    j=1, title="Scheme 1", 
    color="red", f=nothing,
    figsize=(800, 600), lw=1, limits=[(0, 1, 0, 1); (0, 1, 0, 1); (0, 1, 0, 1)],
    ylabs=true, yticks=true,
    xlabs=true, xticks=true)

###    CaM0_m = mean(hcat(CaM_tup[1]...) ./ cam_conc, dims=2) 
###    CaMp_m = mean(hcat(CaM_tup[2]...) ./ cam_conc, dims=2)
###    CaMf_m = mean(hcat(CaM_tup[3]...) ./ cam_conc, dims=2)
    CaM0_m = hcat(CaM_tup[1]) * 1e6
    CaMp_m = hcat(CaM_tup[2]) * 1e6
    CaMf_m = hcat(CaM_tup[3]) * 1e6
 
    if f === nothing
        f = CairoMakie.Figure(size=figsize, fontsize=8)
    end


###    Plots.plot!(t, CaM0_m, alpha=0.20, lw=3, color=color, label=label*" CaM0" )
###    Plots.plot!(t, CaMp_m, alpha=0.60, lw=3, color=color, label=label*" CaM partial" )
###    Plots.plot!(t, CaMf_m, alpha=1.00, lw=3, color=color, label=label*" CaM full" )

    ax1 = CairoMakie.Axis(f[2, j], limits=limits[1], titlesize=8, #, titlealign=:left,
        ###title=title, xticks=[0, 1000],
        xticks=[0, 1000],
        xticklabelsize=8, yticklabelsize=8,
        yticklabelsvisible=ylabs, yticksvisible=yticks,
        xticklabelsvisible=false, xticksvisible=false,
        xgridvisible=false, xticksize=3)
##    hidedecorations!(ax1, ticklabels=false)
    hidespines!(ax1, :t, :r)
    for l in CaM0_m
        CairoMakie.lines!(ax1, t, l, linewidth=lw, color=color, label=title)
    end

    ax2 = CairoMakie.Axis(f[3, j], limits=limits[2],
        xticks=[0, 1000],
        xticklabelsize=8, yticklabelsize=8,
        yticklabelsvisible=ylabs, yticksvisible=yticks,
        xticklabelsvisible=false, xticksvisible=false,
        xgridvisible=false, xticksize=3)
##    hidedecorations!(ax2, ticks=false)
    hidespines!(ax2, :t, :r)
    for l in CaMp_m
        CairoMakie.lines!(ax2, t, l, linewidth=lw, color=color)
    end

    ax3 = CairoMakie.Axis(f[4, j], limits=limits[3],
        xticks=[0, 1000],
        xticklabelsize=8, yticklabelsize=8,
        yticklabelsvisible=ylabs, yticksvisible=yticks,
        xticklabelsvisible=xlabs, xticksvisible=xticks,
        xgridvisible=false, xticksize=3)
##    hidedecorations!(ax3, ticks=false)
    hidespines!(ax3, :t, :r)
    for l in CaMf_m
        CairoMakie.lines!(ax3, t, l, linewidth=lw, color=color)
    end

    return f
end



function calc_integration(fpms, f_range, amp, getter)

    lines = map(fpms) do j
       
        partial_j   = Vector{Float64}()
        full_j      = Vector{Float64}()

        for ff in f_range
            d_ff = gen_f_data(ff, amp)

            _, p_ff, f_ff = getter(j, d_ff)
            i_p = integrate(d_ff[1].time, p_ff)
            i_f = integrate(d_ff[1].time, f_ff)

            push!(partial_j, i_p)
            push!(full_j, i_f)
        end

        (partial_j, full_j)
    end

    return lines 
end



###function plot_int_lines(f::CairoMakie.Figure, freqs, lines, i, j, col, lw, ticks, ticklabs, lims)
###
###    ax = CairoMakie.Axis(f[i, j],
###    xticklabelsize=8, yticklabelsize=8,
###    xticklabelsvisible=ticklabs[1], xticksvisible=ticks[1],
###    yticklabelsvisible=ticklabs[2], yticksvisible=ticks[2],
###    limits=lims[i], ##yscale=log10
###    )
###
###    for l in lines
###        CairoMakie.lines!(ax, freqs, l[i], color=col, linewidth=lw)
###    end
###    return ax
###end
###
###function plot_int_lines(f::CairoMakie.Axis, freqs, lines, i, j, col, lw)
###    for l in lines
###        CairoMakie.lines!(f, freqs, l[i], color=col, linewidth=lw)
###    end
###end

function plot_int_lines(f::CairoMakie.Figure, freqs, lines, i, j, col, la, ra, lw, ticks, ticklabs, lims)

    ax = CairoMakie.Axis(f[i, j],
    xticklabelsize=8, yticklabelsize=8,
    xticklabelsvisible=ticklabs[1], xticksvisible=ticks[1],
    yticklabelsvisible=ticklabs[2], yticksvisible=ticks[2],
    limits=lims[i], ##yscale=log10
    )

    mat = hcat([ff[i] for ff in lines]...)

    CairoMakie.lines!(ax, freqs, median.(eachrow(mat)), color=(col, la), linewidth=lw)
    CairoMakie.band!(ax, freqs, quantile.(eachrow(mat), 0.025), quantile.(eachrow(mat), 0.975), color=(col, ra))

    return ax
end

function plot_int_lines(f::CairoMakie.Axis, freqs, lines, i, j, col, la, ra, lw)
    mat = hcat([ff[i] for ff in lines]...)

    CairoMakie.lines!(f, freqs, median.(eachrow(mat)), color=(col, la), linewidth=lw)
    CairoMakie.band!(f, freqs, quantile.(eachrow(mat), 0.025), quantile.(eachrow(mat), 0.975), color=(col, ra))
end