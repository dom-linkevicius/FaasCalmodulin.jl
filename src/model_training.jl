AVAILABLE_MODELS = (
    Blackwell = Blackwell_scheme,
    Bhalla = Bhalla_scheme,
    Shifman = Shifman_scheme,
    Pepke_m2 = Pepke_m2_scheme,
    Faas = Faas_scheme,
    Pepke_m1 = Pepke_m1_scheme,
    Byrne = Byrne_scheme,
)


FIXED_PARAMSETS = (
    Blackwell = (Blackwell_params, Blackwell_const_names),
    Bhalla = (Bhalla_params, Bhalla_const_names),
    Shifman = (Shifman_params, Shifman_const_names),
    Pepke_m2 = (Pepke_params, Faas_const_names),
    Faas = (Faas_params, Faas_const_names),
    Pepke_m1 = (Pepke_params, Faas_const_names),
    Byrne = (Byrne_params, Byrne_const_names),
)


function model_setup(scheme; fixed = false)

    if scheme in keys(AVAILABLE_MODELS)

        model = AVAILABLE_MODELS[scheme]

        if fixed
            init_p, constantcoefs = FIXED_PARAMSETS[scheme]
        else
            init_p = sample_params(model)
            constantcoefs = ()
        end
        return model, init_p, constantcoefs

    else
        ThrowError("Given model_type", scheme, " is not supported.")
    end
end


function faas_splitting()

    A_idx = shuffle(1:13)
    B_idx = shuffle(14:24)
    C_idx = shuffle([25:26; 29:35]) ### two excluded because they are numerically the same
    D_idx = shuffle(36:51)
    E_idx = shuffle(52:65)
    F_idx = shuffle(66:80)
    G_idx = shuffle(81:94)

    train_idx =
        [A_idx[1:7]; B_idx[1:7]; C_idx[1:7]; D_idx[1:7]; E_idx[1:7]; F_idx[1:7]; G_idx[1:7]]
    val_idx = [
        A_idx[8:10]
        B_idx[8:9]
        C_idx[8:8]
        D_idx[8:12]
        E_idx[8:11]
        F_idx[8:11]
        G_idx[8:11]
    ]
    test_idx = [
        A_idx[11:13]
        B_idx[10:11]
        C_idx[9:9]
        D_idx[13:16]
        E_idx[12:14]
        F_idx[12:15]
        G_idx[12:14]
    ]

    return train_idx, val_idx, test_idx
end


function train_n_models(
    scheme,
    faas_pop,
    shif_pop,
    n,
    seeds;
    fixed = false,
    plot = false,
    sub_range = nothing,
    alg = JointMAP(),
    kwargs...,
)

    n_fpms = []
    val_idxs = []
    test_idxs = []

    if sub_range !== nothing && faas_pop !== nothing
        sub_faas_pop = subsample_start(faas_pop, sub_range)
    else
        sub_faas_pop = faas_pop
    end

    for i = 1:n

        Random.seed!(seeds[i]) ### this is for sampling param reproducibility
        model, init_p, constantcoefs = model_setup(scheme; fixed = fixed)

        if faas_pop === nothing && shif_pop === nothing
            ThrowError("Both data sets are nothing")
        elseif faas_pop === nothing && shif_pop !== nothing
            train_pop = shif_pop
        end

        if faas_pop !== nothing

            Random.seed!(seeds[i]) ### this is for sampling train/val/test data set reproducibility
            train_idx, val_idx, test_idx = faas_splitting()

            push!(val_idxs, val_idx)
            push!(test_idxs, test_idx)

            if shif_pop === nothing
                train_pop = sub_faas_pop[train_idx]
            else
                train_pop = [sub_faas_pop[train_idx]; shif_pop]
            end
        end

        fpm_i = fit(model, train_pop, init_p, alg; constantcoef = constantcoefs, kwargs...)
        push!(n_fpms, fpm_i)
    end

    if plot && faas_pop !== nothing
        p_train = plot_by_condition(n_fpms)

        pops_val = [faas_pop[i] for i in val_idxs]
        p_val = plot_by_condition(n_fpms; pops = pops_val)

        return n_fpms, p_train, p_val, val_idxs, test_idxs
    end

    return n_fpms, nothing, nothing, val_idxs, test_idxs
end


function load_fpms(folder, fpm_names, suffix)

    fpms = map(string.(fpm_names)) do name
        fpms_i = Serialization.deserialize(folder * "fpms_" * name * suffix)
        if name == "blackwell"
            println("Removing a failed run from Blackwell")
            fpms_i = fpms_i[[1:12; 14:20]]
        end
        fpms_i
    end

    val_idxs = Serialization.deserialize(folder * "val_idxs_" * suffix)
    test_idxs = Serialization.deserialize(folder * "test_idxs_" * suffix)
    fpms_nt = NamedTuple{fpm_names}(fpms)

    return fpms_nt, val_idxs, test_idxs

end


function save_fpms(folder::String, fpm_nt::NamedTuple, val, test, suffix::String)

    map(keys(fpm_nt)) do key
        Serialization.serialize(folder * "fpms_" * string(key) * suffix, fpm_nt[key])
    end

    Serialization.serialize(folder * "val_idxs_" * suffix, val)
    Serialization.serialize(folder * "test_idxs_" * suffix, test)
end


function save_non_fpms_outputs(
    folder::String,
    fpm_nt::NamedTuple,
    test_loss_nt::NamedTuple,
    suffix::String,
)

    map(keys(fpm_nt)) do key
        params = [(params = i,) for i in coef.(fpm_nt[key])]
        losses = [(test_loss = i,) for i in test_loss_nt[key]]
        arr_nt = [merge(i, j) for (i, j) in zip(params, losses)]

        nt = NamedTuple{Tuple([Symbol("run" * string(i)) for i = 1:length(params)])}(
            Tuple(arr_nt),
        )
        Serialization.serialize(folder * "params_" * string(key) * suffix, nt)
    end
end

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
###

###
###
###function process_fpms(fpms, pop)
###
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
###function plot_preds_sims(preds, sims, label, color; 
###    legendfontsize=12,
###    labelfontsize=12,
###    titlefontsize=16,
###    tickfontsize=8,
###    rib=false,
###    alpha=1.0
###    )
###
###    for i in 1:length(preds[1])
###        d = zip([pre[i] for pre in preds] , [sim[i] for sim in sims])
###        arr = map(d) do (pre,sim)
###            pred = (sim.dynamics.OGB5 + 39.364 * sim.dynamics.CaOGB5) / (pre.OGB5₀ + 39.364 * pre.CaOGB5₀)
###        end
###
###        
######        if i == 1
######            Plots.plot!(sims[1][i].subject.time, arr[1], color=color, fillalpha=0.2, linewidth=2, alpha=alpha, label=label,
######                legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)    
######        else
######            Plots.plot!(sims[1][i].subject.time, arr[1], color=color, fillalpha=0.2, linewidth=2, alpha=alpha, label=false,
######                legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)    
######        end
######        for a in arr[2:end]
######            Plots.plot!(sims[1][i].subject.time, a, color=color, linewidth=2, alpha=alpha, label=false)       
######        end
###
###        arr = hcat(arr...)
###        med = median(arr, dims=2)
###        sd  = std(arr, dims=2)
###        if i == 1 && rib == true
###            Plots.plot!(sims[1][i].subject.time, med, ribbon=sd, color=color, fillalpha=0.2, linewidth=2, alpha=0.7, label=label,
###            legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)    
###        elseif i == 1 && rib == false
###            Plots.plot!(sims[1][i].subject.time, med, color=color, fillalpha=0.2, linewidth=2, alpha=0.7, label=label,
###            legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)           
###        elseif i !== 1 && rib == true
###            Plots.plot!(sims[1][i].subject.time, med, ribbon=sd, color=color, fillalpha=0.2, linewidth=2, alpha=0.7, label=false,
###            legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)
###        else
###            Plots.plot!(sims[1][i].subject.time, med, color=color, fillalpha=0.2, linewidth=2, alpha=0.7, label=false,
###            legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)    
###        end
###    end
###end
###
###
###
###function plot_mean_std(fpms, pops, colors; names=nothing, add_legend=false, std=false, im=nothing,
###    legendfontsize=12, titlefontsize=16, labelfontsize=12, tickfontsize=8, alpha_val=1.0, 
###    ylim=[0, 18])
###
###    if im !== nothing
###        im = Plots.plot(im, axis=false, grid=false, title="Scheme", titlefontsize=titlefontsize)
###    else
###        im = Plots.plot(axis=false, grid=false, title="Scheme", titlefontsize=titlefontsize)
###    end
###
###    plots = [im,]
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
###    for i in 1:length(pops)
###
###        pop_i = pops[i]
###
###        splits = split(pop_i[1].id, "_")
###        title_i = splits[1] * "_" * splits[2]
###        title_i = title_dict[title_i]
###
###        pre_o, sims_o = process_fpms(fpms[1], pop_i)
###        pre_p, sims_p = process_fpms(fpms[2], pop_i)
###
###        if i == 1 || i == 4
###            p_i = Plots.plot(ylabel="ΔF/F₀", title=title_i, titlefontsize=titlefontsize)
###            if i == 1
###                Plots.plot!(ylim=ylim)
###            end
###        elseif i in [4, 5, 6, 7]
###            p_i = Plots.plot(xlabel="Time (ms)", title=title_i, titlefontsize=titlefontsize)
###        else
###            p_i = Plots.plot(title=title_i, titlefontsize=titlefontsize)
###        end
###
###        if i == 1
###            Plots.plot!(sims_o[1][1].subject.time, sims_o[1][1].subject.observations.F_F0, color="black", label="Faas et al. (2011) data", linewidth=3)
###            for sim_i in sims_o[1][2:end]
###                Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false, linewidth=2)
###            end
###        else
###             for sim_i in sims_o[1]
###                Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false, linewidth=2)
###            end
###        end
###
###
###        if i == 1
###            plot_preds_sims(pre_o, sims_o, "Our parameters", colors[1];
###            legendfontsize=legendfontsize, labelfontsize=labelfontsize, tickfontsize=tickfontsize, titlefontsize=titlefontsize, alpha=alpha_val)
###            plot_preds_sims(pre_p, sims_p, "Published parameters", colors[2];
###            legendfontsize=legendfontsize, labelfontsize=labelfontsize, tickfontsize=tickfontsize, titlefontsize=titlefontsize, alpha=alpha_val)
###        else
###            plot_preds_sims(pre_o, sims_o, false, colors[1];
###            legendfontsize=legendfontsize, labelfontsize=labelfontsize, tickfontsize=tickfontsize, titlefontsize=titlefontsize, alpha=alpha_val)
###            plot_preds_sims(pre_p, sims_p, false, colors[2];
###            legendfontsize=legendfontsize, labelfontsize=labelfontsize, tickfontsize=tickfontsize, titlefontsize=titlefontsize, alpha=alpha_val)
###        end
###
###        push!(plots, p_i)
###    end
###
###    ret_plot = Plots.plot(plots..., 
###        size=(1000, 600), 
###        layout=(2,4), 
###        dpi=600, 
###        margin=-4mm, 
###        left_margin=8mm, 
###        right_margin=1mm,
###        bottom_margin=5mm,
###        top_margin=1mm)
###
###    return ret_plot
###end
###
###
###
###function plot_fens(fpms, pop, colors; names=nothing, add_legend=false)
###
###    im1 = Plots.plot(load("Simple_scheme_bigger.png"), axis=false, grid=false, title="Scheme 1")
###    im2 = Plots.plot(load("Faas_scheme.png"), axis=false, grid=false, title="Scheme 2")
###    im3 = Plots.plot(load("Byrne_scheme.png"), axis=false, grid=false, title="Scheme 3")
###
###    plots = []
###
###    pre1, sims1 = process_fpms(fpms[1], pop)
###    pre2, sims2 = process_fpms(fpms[2], pop)
###    pre3, sims3 = process_fpms(fpms[3], pop)
###    pre4, sims4 = process_fpms(fpms[4], pop)
###    pre5, sims5 = process_fpms(fpms[5], pop)
###    pre6, sims6 = process_fpms(fpms[6], pop)
###
###
###    p1 = Plots.plot(xlabel="Time (ms)", ylabel="ΔF/F0", ylim=[1, 14])
###    for sim_i in sims1[1]
###        Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false)
###    end
###
###    plot_preds_sims(pre1, sims1, false, colors[1])
###    plot_preds_sims(pre2, sims2, false, colors[2])
###
###
###    p2 = Plots.plot(xlabel="Time (ms)", ylim=[1, 14])
###
###    for i in 1:length(sims1[1])
###        sim_i = sims1[1][i]
###        if i == 1
###            Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label="Data")
###        else
###            Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false)
###        end
###    end
###    plot_preds_sims(pre3, sims3, "Ours", colors[1])
###    plot_preds_sims(pre4, sims4, "Published", colors[2])
###
###
###    p3 = Plots.plot(xlabel="Time (ms)", ylim=[1, 14])
###    for sim_i in sims1[1]
###        Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false)
###    end
###    plot_preds_sims(pre5, sims5, false, colors[1])
###    plot_preds_sims(pre6, sims6, false, colors[2])
###
###    ret_plot = Plots.plot(im1, im2, im3, p1, p2, p3, 
###        size=(600, 800), 
###        layout=(2,3), 
###        dpi=600, 
###        margin=-4mm, 
###        left_margin=2mm, 
###        right_margin=1mm,
###        bottom_margin=-3mm)
###
###    return ret_plot
###end
###
###
###function idx_to_val_vec(full_pop, idxs)
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
###
###
###function plot_exp_data(groups; figsize=(800, 600), ylim=[0, 15], lw=2)
###
###
###
###    f = CairoMakie.Figure(size=figsize, fontsize=8)
###
###    ax1 = CairoMakie.Axis(f[1,1], 
###        ylabel="ΔF/F₀", 
###        title="Group A", 
###        limits = (-3, 43, ylim[1], ylim[2]), 
###        yticklabelsize=8, 
###        xticklabelsize=8,
###        )
###    for i in groups[1]
###        CairoMakie.lines!(ax1, i.time, i.observations.F_F0, color="black", linewidth=lw)
###    end
###
###    ax2 = CairoMakie.Axis(f[1,2], title="Group B", xlabel="Time (ms)", limits = (-3, 43, ylim[1], ylim[2]), yticklabelsize=8, xticklabelsize=8)
###    for i in groups[2]
###        CairoMakie.lines!(ax2, i.time, i.observations.F_F0, color="black", linewidth=lw)
###    end
###
###    ax3 = CairoMakie.Axis(f[1,3], title="Group G", limits = (-3, 43, ylim[1], ylim[2]), yticklabelsize=8, xticklabelsize=8)
###    for i in groups[3]
###        CairoMakie.lines!(ax3, i.time, i.observations.F_F0, color="black", linewidth=lw)
###    end
###
###    rowsize!(f.layout, 1, Relative(1))
###
###    return f
###end
###
###
###function plot_subsampled_data(groups_full, groups_subsampled; figsize=(800, 600), ylim=[0, 15], lw=2, la=0.5)
###
###    f = CairoMakie.Figure(size=figsize, fontsize=8)
###
###    ax1 = CairoMakie.Axis(f[1,1], 
###        title="Group A", 
###        yticklabelsize=8, 
###        xticklabelsize=8,
###        )
###    for (i,j) in zip(groups_full[1], groups_subsampled[1])
###        CairoMakie.lines!(ax1, i.time, i.observations.F_F0, color=(:black, la), linewidth=lw, label="Full data")
###        CairoMakie.lines!(ax1, j.time, j.observations.F_F0, color=(:red, la), linewidth=lw, label="Subsampled data")
###    end
###    ax2 = CairoMakie.Axis(f[1,2], 
###        title="Group B", 
###        yticklabelsize=8, 
###        xticklabelsize=8,
###        )
###    for (i,j) in zip(groups_full[2], groups_subsampled[2])
###        CairoMakie.lines!(ax2, i.time, i.observations.F_F0, color=(:black, la), linewidth=lw)
###        CairoMakie.lines!(ax2, j.time, j.observations.F_F0, color=(:red, la), linewidth=lw)
###    end
###    ax3 = CairoMakie.Axis(f[2,1], 
###        title="Group C", 
###        yticklabelsize=8, 
###        xticklabelsize=8,
###        )
###    for (i,j) in zip(groups_full[3], groups_subsampled[3])
###        CairoMakie.lines!(ax3, i.time, i.observations.F_F0, color=(:black, la), linewidth=lw)
###        CairoMakie.lines!(ax3, j.time, j.observations.F_F0, color=(:red, la), linewidth=lw)
###    end
###    ax4 = CairoMakie.Axis(f[2,2], 
###        title="Group D", 
###        yticklabelsize=8, 
###        xticklabelsize=8,
###        )
###    for (i,j) in zip(groups_full[4], groups_subsampled[4])
###        CairoMakie.lines!(ax4, i.time, i.observations.F_F0, color=(:black, la), linewidth=lw)
###        CairoMakie.lines!(ax4, j.time, j.observations.F_F0, color=(:red, la), linewidth=lw)
###    end
###    ax5 = CairoMakie.Axis(f[3,1], 
###        title="Group E", 
###        yticklabelsize=8, 
###        xticklabelsize=8,
###        )
###    for (i,j) in zip(groups_full[5], groups_subsampled[5])
###        CairoMakie.lines!(ax5, i.time, i.observations.F_F0, color=(:black, la), linewidth=lw)
###        CairoMakie.lines!(ax5, j.time, j.observations.F_F0, color=(:red, la), linewidth=lw)
###    end
###    ax6 = CairoMakie.Axis(f[3,2], 
###        title="Group F", 
###        yticklabelsize=8, 
###        xticklabelsize=8,
###        )
###    for (i,j) in zip(groups_full[6], groups_subsampled[6])
###        CairoMakie.lines!(ax6, i.time, i.observations.F_F0, color=(:black, la), linewidth=lw)
###        CairoMakie.lines!(ax6, j.time, j.observations.F_F0, color=(:red, la), linewidth=lw)
###    end   
###    ax7 = CairoMakie.Axis(f[4,1], 
###        title="Group G", 
###        yticklabelsize=8, 
###        xticklabelsize=8,
###        )
###    for (i,j) in zip(groups_full[7], groups_subsampled[7])
###        CairoMakie.lines!(ax7, i.time, i.observations.F_F0, color=(:black, la), linewidth=lw)
###        CairoMakie.lines!(ax7, j.time, j.observations.F_F0, color=(:red, la), linewidth=lw)
###    end
###
###    CairoMakie.Legend(f[4, 2], ax1, unique=true)
###
###    CairoMakie.Label(f[:,0], "ΔF/F₀", rotation=pi/2)
###    CairoMakie.Label(f[5,:], "Time (ms)")
###
###    colsize!(f.layout, 1, Relative(1/2))
###    colsize!(f.layout, 2, Relative(1/2))
###
###    return f
###end
###
###
###
###function bigplot(fpms, pop, colors, rowslist; 
###    names=nothing, figsize=(800, 600))
###
###    plots = []
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
###    f = CairoMakie.Figure(size=figsize, fontsize=8)
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
###    return f
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
###
######            push!(all_N1_c, maximum(ss[end-5,:])/CaM_t)
######            push!(all_N2_c, maximum(ss[end-4,:])/CaM_t)
######            push!(all_N3_c, maximum(ss[end-3,:])/CaM_t)
######
######            push!(all_C1_c, maximum(ss[end-2,:])/CaM_t)
######            push!(all_C2_c, maximum(ss[end-1,:])/CaM_t)
######            push!(all_C3_c, maximum(ss[end-0,:])/CaM_t)
######
######            push!(all_N1_t, tt[argmax(ss[end-5,:])])
######            push!(all_N2_t, tt[argmax(ss[end-4,:])])
######            push!(all_N3_t, tt[argmax(ss[end-3,:])])
######
######            push!(all_C1_t, tt[argmax(ss[end-2,:])])
######            push!(all_C2_t, tt[argmax(ss[end-1,:])])
######            push!(all_C3_t, tt[argmax(ss[end-0,:])])
###        end
###    end
###
###    ###return hcat(all_N1_c, all_N2_c, all_N3_c, all_C1_c, all_C2_c, all_C3_c), hcat(all_N1_t, all_N2_t, all_N3_t, all_C1_t, all_C2_t, all_C3_t), all_catot
###    return (n0 = all_N0_c, n1 = all_N1_c, n2 = all_N2_c, n3 = all_N3_c, 
###            c0 = all_C0_c, c1 = all_C1_c, c2 = all_C2_c, c3 = all_C3_c), 
###           (n0 = all_N0_t, n1 = all_N1_t, n2 = all_N2_t, n3 = all_N3_t, 
###            c0 = all_C0_t, c1 = all_C1_t, c2 = all_C2_t, c3 = all_C3_t), 
###            all_catot
###end
###
###


###function fpms_to_table_str(fpms)
###
###    pnames = setdiff(keys(coef(fpms[1])), (:μ, :Ω, :σ_faas, :σ_shifman))
###    str_names = string.(pnames)
###
###    str_table = map(1:length(fpms)) do i
###        f = fpms[i]
###        str_j = "\n \\hline \n Seed " * string(i) * " & "
###
###        map(pnames[1:end-1]) do p
###            str_j *= string(round(coef(f)[p]; digits=2)) * " & "
###        end
###        str_j *= string(round(coef(f)[pnames[end]]; digits=2)) * " "
###        str_j *= "\\\\"
###    end
###    return pnames, reduce(*, str_table)
###end
