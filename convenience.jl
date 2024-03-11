function model_setup(model_type)    

    if model_type == "Blackwell"
        model = Blackwell_scheme
        init_p = sample_params(Blackwell_scheme)
        ###constantcoefs = (:tv_mα, :tv_α₀)
        constantcoefs = ()

    elseif model_type == "Blackwell_fixed"
        model = Blackwell_scheme
        init_p = Blackwell_params
        constantcoefs = Blackwell_const_names

    elseif model_type == "Pepke_m2"
        model = Pepke_m2_scheme
        init_p = sample_params(model)
        ###constantcoefs = (:tv_mα, :tv_α₀)
        constantcoefs = ()

    elseif model_type == "Pepke_m2_fixed"
        model = Pepke_m2_scheme 
        init_p = Pepke_M2_params
        constantcoefs = Pepke_M2_const_names

    elseif model_type == "Blackwell_TR"
        model = Blackwell_scheme_TR
        init_p = sample_params(model)
        ###constantcoefs = (:tv_mα, :tv_α₀)
        constantcoefs = ()
    
    elseif  model_type == "Faas"
        model = Faas_scheme
        init_p = sample_params(model)
        ###constantcoefs = (:tv_mα, :tv_α₀)
        constantcoefs = ()

    elseif  model_type == "Faas_fixed"
        model = Faas_scheme
        init_p = Faas_params
        constantcoefs = Faas_const_names

    elseif  model_type == "Pepke_m1_fixed"
        model = Faas_scheme
        init_p = Pepke_M1_params
        constantcoefs = Faas_const_names

    elseif  model_type == "Byrne"
        model = Byrne_scheme
        init_p = sample_params(model)
        ###constantcoefs = (:tv_mα, :tv_α₀)
        constantcoefs = ()
    
    elseif  model_type == "Byrne_fixed"
        model = Byrne_scheme
        init_p = Byrne_params
        constantcoefs = Byrne_const_names   
    end

    return model, init_p, constantcoefs

end


function custom_derived(fpm, idx, color)

    all_pre  = [fpm.model.pre(coef(fpm), empirical_bayes(fpm)[i], fpm.data[i]).ret for i in 1:length(fpm.data)]
    all_sims = simobs(fpm);

    Plots.plot()
    for (pr, sm) in zip(all_pre[idx], all_sims[idx])
        Plots.plot!(sm.subject.time, sm.subject.observations.F_F0, color="black", label=false)

        pred = (sm.dynamics.OGB5 + 39.364 * sm.dynamics.CaOGB5) / (pr.OGB5₀ + 39.364 * pr.CaOGB5₀)
        Plots.plot!(sm.subject.time, pred, ls=:dash, color=color, label=false)
    end
    Plots.plot!()
end


function plot_training_by_condition(fpms, colors; names=nothing, add_legend=false)

    all_pre = []
    all_sims = []

    for fpm in fpms
        pre_i = [fpm.model.pre(coef(fpm), empirical_bayes(fpm)[i], fpm.data[i]).ret for i in 1:length(fpm.data)]
        sim_i = simobs(fpm)

        push!(all_pre, pre_i)
        push!(all_sims, sim_i)
    end

    uids = unique([split(i.id, "_")[1] * "_" * split(i.id, "_")[2] for i in fpms[1].data])
    plots = []

    for idx in uids
        p = Plots.plot(title=idx, xlabel="Time (ms)", ylabel="F/F0")

        for sims in all_sims

            true_data = [i for i in sims if occursin(idx, i.subject.id)]
            for sm in true_data
                Plots.plot!(sm.subject.time, sm.subject.observations.F_F0, color="black", label=false)
            end
    
            for i in 1:length(fpms)
                sub_ids = findall(occursin.(idx, [i.subject.id for i in all_sims[i]]))
                
                _pr = all_pre[i][sub_ids]
                _sm = all_sims[i][sub_ids]
                for (pr, sm) in zip(_pr, _sm)
                    pred = (sm.dynamics.OGB5 + 39.364 * sm.dynamics.CaOGB5) / (pr.OGB5₀ + 39.364 * pr.CaOGB5₀)
                    Plots.plot!(sm.subject.time, pred, color=colors[i], label=false, linewidth=3)
                end
            end
        end
        push!(plots, p)
    end

    if add_legend
        rmses = [mean(training_rmse(i)) for i in fpms]

        p0 = Plots.plot([], [], label="True data", color="black", legend=:inside, legendcolumns=1, legendfontsize=15, showaxis=false, grid=false)
        for (n, c, m) in zip(names, colors, rmses) 
            Plots.plot!([], [], label=n * " RMSE="*string(round(m; sigdigits=5)), color=c) #, ls=:dash)
        end
        push!(plots, p0)
    end

    ret_plot = Plots.plot(plots..., size=(1600, 1200), layout=(2,4))
    return ret_plot
end


function plot_unseen_by_condition(fpms, pops, colors; names=nothing, add_legend=false)

    all_pre = []
    all_sims = []

    for (fpm,pop) in zip(fpms, pops)
        ebes = empirical_bayes(fpm.model, pop, coef(fpm), LaplaceI(), diffeq_options=(;alg=Rodas5P(), abstol=1e-10))
        pre_i = [fpm.model.pre(coef(fpm), ebes[i], pop[i]).ret for i in 1:length(pop)]
        sim_i = simobs(fpm.model, 
            pop,
            coef(fpm),
            ebes,
            diffeq_options=(;alg=Rodas5P(), abstol=1e-10),
            )

        push!(all_pre, pre_i)
        push!(all_sims, sim_i)
    end

    uids = unique([split(i.id, "_")[1] * "_" * split(i.id, "_")[2] for i in pops[1]])
    plots = []

    for idx in uids
        p = Plots.plot(title=idx, xlabel="Time (ms)", ylabel="F/F0")
        for sims in all_sims
    
            true_data = [i for i in sims if occursin(idx, i.subject.id)]
            for sm in true_data
                Plots.plot!(sm.subject.time, sm.subject.observations.F_F0, color="black", label=false)
            end
    
            for i in 1:length(fpms)
                sub_ids = findall(occursin.(idx, [i.subject.id for i in all_sims[i]]))
                
                _pr = all_pre[i][sub_ids]
                _sm = all_sims[i][sub_ids]
                for (pr, sm) in zip(_pr, _sm)
                    pred = (sm.dynamics.OGB5 + 39.364 * sm.dynamics.CaOGB5) / (pr.OGB5₀ + 39.364 * pr.CaOGB5₀)
                    Plots.plot!(sm.subject.time, pred, color=colors[i], label=false, linewidth=3)
                end
            end
        end
        push!(plots, p)
    end

    if add_legend

        rmses = [mean(unseen_rmse(i, j)) for (i,j) in zip(fpms, pops)]

        p0 = Plots.plot([], [], label="True data", color="black", legend=:inside, legendcolumns=1, legendfontsize=15, showaxis=false, grid=false)
        for (n, c, m) in zip(names, colors, rmses) 
            Plots.plot!([], [], label=n * " RMSE="*string(round(m; sigdigits=5)), color=c) ###, ls=:dash)
            ###Plots.plot!([], [], label=n, color=c) ###, ls=:dash)
        end
        push!(plots, p0)
    end

    ret_plot = Plots.plot(plots..., size=(1600, 1200), layout=(2,4))
    return ret_plot
end


function training_rmse(fpm)

   all_pre  = [fpm.model.pre(coef(fpm), empirical_bayes(fpm)[i], fpm.data[i]).ret for i in 1:length(fpm.data)];
   all_sims = simobs(fpm);

   preds = map(zip(all_pre, all_sims)) do j
        F_F0 = (j[2].dynamics.OGB5 + 39.364 * j[2].dynamics.CaOGB5) / (j[1].OGB5₀ + 39.364 * j[1].CaOGB5₀)
   end

   RMSEs = map(zip(preds, fpm.data)) do d
        rmse_i = sqrt.(mean((d[1] .- d[2].observations.F_F0).^2))
   end

   return RMSEs
end

function unseen_rmse(fpm, pop)

    ebes = empirical_bayes(fpm.model, pop, coef(fpm), LaplaceI(), diffeq_options=(;alg=Rodas5P(), abstol=1e-10))
    all_pre  = [fpm.model.pre(coef(fpm), ebes[i], pop[i]).ret for i in 1:length(pop)]
    all_sims = simobs(fpm.model, 
            pop,
            coef(fpm),
            ebes,
            diffeq_options=(;alg=Rodas5P(), abstol=1e-10)
    )
   preds = map(zip(all_pre, all_sims)) do j
        F_F0 = (j[2].dynamics.OGB5 + 39.364 * j[2].dynamics.CaOGB5) / (j[1].OGB5₀ + 39.364 * j[1].CaOGB5₀)
   end

   RMSEs = map(zip(preds, pop)) do d
        rmse_i = sqrt.(mean((d[1] .- d[2].observations.F_F0).^2))
   end

   return RMSEs
end


function my_aic(fpm::Pumas.AbstractFittedPumasModel; pop=nothing)
    if pop == nothing
        return 2*length(coef(fpm)) - 2*loglikelihood(fpm)
    else
        ll = loglikelihood(fpm.model, pop, coef(fpm), LaplaceI(), diffeq_options=(;alg=Rodas5P(), abstol=1e-16))
        return 2*length(coef(fpm)) - 2*ll
    end
end


###function aic_boxplots(fpms, names; pop=nothing)
function aic_boxplots(triples)

    all_x = []
    all_y = []

    for trip in triples
        fpm = trip[1]
        pop = vcat(trip[2]...)
        name = trip[3]
        
        y_i = my_aic(fpm; pop=pop)

        push!(all_x, name)
        push!(all_y, y_i)
    end

    return all_x, all_y

end



function train_plot(model_type, 
    full_pop, 
    n,
    seeds;
    sub_range=nothing,
    alg=JointMAP, 
    optim_options=(; iterations=2000, store_trace=true),
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
###    diffeq_options=(;alg=RadauIIA5(), abstol=1e-14),
)

    n_fpms = []

    if sub_range != nothing
        sub_pop = subsample_start(full_pop, sub_range)
    else
        sub_pop = full_pop
    end

    val_idxs  = []
    test_idxs = []

    for i in 1:n

        Random.seed!(seeds[i]) ### this is for sampling param reproducibility

        model, init_p, constantcoefs = model_setup(model_type)

        Random.seed!(seeds[i]) ### this is for sampling train/val/test data set reproducibility

        A_idx = shuffle(1:13)
        B_idx = shuffle(14:24)
        C_idx = shuffle(25:35)
        D_idx = shuffle(36:51)
        E_idx = shuffle(52:65)
        F_idx = shuffle(66:80)
        G_idx = shuffle(81:94)

        train_idx = [A_idx[1:7]; B_idx[1:7]; D_idx[1:7]; E_idx[1:7]; F_idx[1:7]; G_idx[1:7]]
        val_idx   = [A_idx[8:10]; B_idx[8:9]; D_idx[8:12]; E_idx[8:11]; F_idx[8:11]; G_idx[8:11]]
        test_idx  = [A_idx[11:13]; B_idx[10:11]; D_idx[13:16]; E_idx[12:14]; F_idx[12:15]; G_idx[12:14]]

        push!(val_idxs, val_idx)
        push!(test_idxs, test_idx)

        fpm_i = fit(model, 
        sub_pop[train_idx], 
        init_p, 
        ###alg(), 
        MAP(FOCE()), 
        optim_options=optim_options,
        diffeq_options=diffeq_options,
        constantcoef=constantcoefs
    );
        push!(n_fpms, fpm_i)
    end

    p_train = plot_training_by_condition(
        n_fpms, 
        palette([:purple, :green], n+1);
        names=["Run " * string(i) for i in 1:n],
        add_legend=true
    )

    pops_val = [full_pop[i] for i in val_idxs]

    p_val = plot_unseen_by_condition(
        n_fpms,    
        pops_val,
        palette([:purple, :green], n+1);
        names=["Run " * string(i) for i in 1:n],
        add_legend=true
    )

    return n_fpms, p_train, p_val, val_idxs, test_idxs
end


function process_fpms(fpms, pop)

    all_pre = []
    all_sims = []

    for j in 1:length(fpms)

        fpm = fpms[j]

        ebes = empirical_bayes(fpm.model, pop, coef(fpm), LaplaceI(), diffeq_options=(;alg=Rodas5P(), abstol=1e-16))
        pre_i = [fpm.model.pre(coef(fpm), ebes[i], pop[i]).ret for i in 1:length(pop)]
        sim_i = simobs(fpm.model, 
            pop,
            coef(fpm),
            ebes,
            diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
        )

        push!(all_pre, pre_i)
        push!(all_sims, sim_i)
    end
    return all_pre, all_sims
end

function plot_preds_sims(preds, sims, label, color; 
    legendfontsize=12,
    labelfontsize=12,
    titlefontsize=16,
    tickfontsize=8,
    rib=false,
    alpha=1.0
    )

    for i in 1:length(preds[1])
        d = zip([pre[i] for pre in preds] , [sim[i] for sim in sims])
        arr = map(d) do (pre,sim)
            pred = (sim.dynamics.OGB5 + 39.364 * sim.dynamics.CaOGB5) / (pre.OGB5₀ + 39.364 * pre.CaOGB5₀)
        end

        
###        if i == 1
###            Plots.plot!(sims[1][i].subject.time, arr[1], color=color, fillalpha=0.2, linewidth=2, alpha=alpha, label=label,
###                legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)    
###        else
###            Plots.plot!(sims[1][i].subject.time, arr[1], color=color, fillalpha=0.2, linewidth=2, alpha=alpha, label=false,
###                legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)    
###        end
###        for a in arr[2:end]
###            Plots.plot!(sims[1][i].subject.time, a, color=color, linewidth=2, alpha=alpha, label=false)       
###        end

        arr = hcat(arr...)
        med = median(arr, dims=2)
        sd  = std(arr, dims=2)
        if i == 1 && rib == true
            Plots.plot!(sims[1][i].subject.time, med, ribbon=sd, color=color, fillalpha=0.2, linewidth=2, alpha=0.7, label=label,
            legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)    
        elseif i == 1 && rib == false
            Plots.plot!(sims[1][i].subject.time, med, color=color, fillalpha=0.2, linewidth=2, alpha=0.7, label=label,
            legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)           
        elseif i !== 1 && rib == true
            Plots.plot!(sims[1][i].subject.time, med, ribbon=sd, color=color, fillalpha=0.2, linewidth=2, alpha=0.7, label=false,
            legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)
        else
            Plots.plot!(sims[1][i].subject.time, med, color=color, fillalpha=0.2, linewidth=2, alpha=0.7, label=false,
            legendfontsize=legendfontsize, labelfontsize=labelfontsize, titlefontsize=titlefontsize, tickfontsize=tickfontsize)    
        end
    end
end



function plot_mean_std(fpms, pops, colors; names=nothing, add_legend=false, std=false, im=nothing,
    legendfontsize=12, titlefontsize=16, labelfontsize=12, tickfontsize=8, alpha_val=1.0, 
    ylim=[0, 18])

    if im !== nothing
        im = Plots.plot(im, axis=false, grid=false, title="Scheme", titlefontsize=titlefontsize)
    else
        im = Plots.plot(axis=false, grid=false, title="Scheme", titlefontsize=titlefontsize)
    end

    plots = [im,]

    title_dict = Dict("#1021_WT"=>"Group A", 
                      "#0611_WTa"=>"Group B",
                      "#0611_WTb"=>"Group C",
                      "#1107_WTa"=>"Group D",
                      "#1107_WTb"=>"Group E",
                      "#1107_WTc"=>"Group F",
                      "#1107_WTd"=>"Group G",
    )

    for i in 1:length(pops)

        pop_i = pops[i]

        splits = split(pop_i[1].id, "_")
        title_i = splits[1] * "_" * splits[2]
        title_i = title_dict[title_i]

        pre_o, sims_o = process_fpms(fpms[1], pop_i)
        pre_p, sims_p = process_fpms(fpms[2], pop_i)

        if i == 1 || i == 4
            p_i = Plots.plot(ylabel="ΔF/F₀", title=title_i, titlefontsize=titlefontsize)
            if i == 1
                Plots.plot!(ylim=ylim)
            end
        elseif i in [4, 5, 6, 7]
            p_i = Plots.plot(xlabel="Time (ms)", title=title_i, titlefontsize=titlefontsize)
        else
            p_i = Plots.plot(title=title_i, titlefontsize=titlefontsize)
        end

        if i == 1
            Plots.plot!(sims_o[1][1].subject.time, sims_o[1][1].subject.observations.F_F0, color="black", label="Faas et al. (2011) data", linewidth=3)
            for sim_i in sims_o[1][2:end]
                Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false, linewidth=2)
            end
        else
             for sim_i in sims_o[1]
                Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false, linewidth=2)
            end
        end


        if i == 1
            plot_preds_sims(pre_o, sims_o, "Our parameters", colors[1];
            legendfontsize=legendfontsize, labelfontsize=labelfontsize, tickfontsize=tickfontsize, titlefontsize=titlefontsize, alpha=alpha_val)
            plot_preds_sims(pre_p, sims_p, "Published parameters", colors[2];
            legendfontsize=legendfontsize, labelfontsize=labelfontsize, tickfontsize=tickfontsize, titlefontsize=titlefontsize, alpha=alpha_val)
        else
            plot_preds_sims(pre_o, sims_o, false, colors[1];
            legendfontsize=legendfontsize, labelfontsize=labelfontsize, tickfontsize=tickfontsize, titlefontsize=titlefontsize, alpha=alpha_val)
            plot_preds_sims(pre_p, sims_p, false, colors[2];
            legendfontsize=legendfontsize, labelfontsize=labelfontsize, tickfontsize=tickfontsize, titlefontsize=titlefontsize, alpha=alpha_val)
        end

        push!(plots, p_i)
    end

    ret_plot = Plots.plot(plots..., 
        size=(1000, 600), 
        layout=(2,4), 
        dpi=600, 
        margin=-4mm, 
        left_margin=8mm, 
        right_margin=1mm,
        bottom_margin=5mm,
        top_margin=1mm)

    return ret_plot
end



function plot_fens(fpms, pop, colors; names=nothing, add_legend=false)

    im1 = Plots.plot(load("Simple_scheme_bigger.png"), axis=false, grid=false, title="Scheme 1")
    im2 = Plots.plot(load("Faas_scheme.png"), axis=false, grid=false, title="Scheme 2")
    im3 = Plots.plot(load("Byrne_scheme.png"), axis=false, grid=false, title="Scheme 3")

    plots = []

    pre1, sims1 = process_fpms(fpms[1], pop)
    pre2, sims2 = process_fpms(fpms[2], pop)
    pre3, sims3 = process_fpms(fpms[3], pop)
    pre4, sims4 = process_fpms(fpms[4], pop)
    pre5, sims5 = process_fpms(fpms[5], pop)
    pre6, sims6 = process_fpms(fpms[6], pop)


    p1 = Plots.plot(xlabel="Time (ms)", ylabel="ΔF/F0", ylim=[1, 14])
    for sim_i in sims1[1]
        Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false)
    end

    plot_preds_sims(pre1, sims1, false, colors[1])
    plot_preds_sims(pre2, sims2, false, colors[2])


    p2 = Plots.plot(xlabel="Time (ms)", ylim=[1, 14])

    for i in 1:length(sims1[1])
        sim_i = sims1[1][i]
        if i == 1
            Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label="Data")
        else
            Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false)
        end
    end
    plot_preds_sims(pre3, sims3, "Ours", colors[1])
    plot_preds_sims(pre4, sims4, "Published", colors[2])


    p3 = Plots.plot(xlabel="Time (ms)", ylim=[1, 14])
    for sim_i in sims1[1]
        Plots.plot!(sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", label=false)
    end
    plot_preds_sims(pre5, sims5, false, colors[1])
    plot_preds_sims(pre6, sims6, false, colors[2])

    ret_plot = Plots.plot(im1, im2, im3, p1, p2, p3, 
        size=(600, 800), 
        layout=(2,3), 
        dpi=600, 
        margin=-4mm, 
        left_margin=2mm, 
        right_margin=1mm,
        bottom_margin=-3mm)

    return ret_plot
end


function idx_to_val_vec(full_pop, idxs)

    A = full_pop[idxs[idxs .> 0  .&& idxs .<= 13]]
    B = full_pop[idxs[idxs .> 13 .&& idxs .<= 24]]
    C = full_pop[idxs[idxs .> 24 .&& idxs .<= 35]]
    D = full_pop[idxs[idxs .> 35 .&& idxs .<= 51]]
    E = full_pop[idxs[idxs .> 51 .&& idxs .<= 65]]
    F = full_pop[idxs[idxs .> 65 .&& idxs .<= 80]]
    G = full_pop[idxs[idxs .> 81 .&& idxs .<= 94]]

    return [A, B, D, E, F, G]
end



function plot_exp_data(groups; figsize=(800, 600), ylim=[0, 15], lw=2)



    f = CairoMakie.Figure(size=figsize, fontsize=12)

    ax1 = CairoMakie.Axis(f[1,1], 
        ylabel=L"\frac{\Delta F}{F_0}", 
        title="Group A", 
        limits = (-3, 43, ylim[1], ylim[2]), 
        yticklabelsize=8, 
        xticklabelsize=8,
        protrusions=0)
    for i in groups[1]
        CairoMakie.lines!(ax1, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    ax2 = CairoMakie.Axis(f[1,2], title="Group B", xlabel="Time (ms)", limits = (-3, 43, ylim[1], ylim[2]), yticklabelsize=8, xticklabelsize=8)
    for i in groups[2]
        CairoMakie.lines!(ax2, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    ax3 = CairoMakie.Axis(f[1,3], title="Group G", limits = (-3, 43, ylim[1], ylim[2]), yticklabelsize=8, xticklabelsize=8)
    for i in groups[3]
        CairoMakie.lines!(ax3, i.time, i.observations.F_F0, color="black", linewidth=lw)
    end

    rowsize!(f.layout, 1, Relative(1))

    return f
end



function bigplot(fpms, pop; 
    names=nothing, figsize=(800, 600))

    plots = []

    title_dict = Dict("#1021_WT"=>"Group A", 
                      "#0611_WTa"=>"Group B",
                      "#0611_WTb"=>"Group C",
                      "#1107_WTa"=>"Group D",
                      "#1107_WTb"=>"Group E",
                      "#1107_WTc"=>"Group F",
                      "#1107_WTd"=>"Group G",
    )

    f = CairoMakie.Figure(size=figsize, fontsize=12)

    for i in 1:length(fpms)

        for j in 1:(length(pop)+2)

            if j in 1:length(pop)

                pop_j = pop[j]

                splits = split(pop_j[1].id, "_")
                title_i = splits[1] * "_" * splits[2]
                title_i = title_dict[title_i]

                ax = CairoMakie.Axis(f[i, j])
                hidedecorations!(ax)
                hidespines!(ax)

                if length(f.content) in 1:6
                    ax.title=title_i
                    ax.titlesize=8
                end


                pre, sims = process_fpms(fpms[i:i], pop_j)

                for sim_i in sims[1]
                    lines!(ax, sim_i.subject.time, sim_i.subject.observations.F_F0, color="black", linewidth=1)
                end

                d = zip(pre[1], sims[1])
                for (pre,sim) in d
                    pred = (sim.dynamics.OGB5 + 39.364 * sim.dynamics.CaOGB5) / (pre.OGB5₀ + 39.364 * pre.CaOGB5₀)
                    lines!(ax, sim.time, pred, color="red", linewidth=1)
                end


            elseif j == length(pop)+1
##                ax = CairoMakie.Axis(f[i, j])
##                hidedecorations!(ax)
##                hidespines!(ax)
            end
        end
    end

    for i in 1:length(fpms) 
        CairoMakie.Label(f[i,7], names[i], fontsize=8)
    end

    CairoMakie.Label(f[:,0], L"\frac{\Delta F}{F_0}", rotation=pi/2)
    CairoMakie.Label(f[length(fpms)+1,:], "Time (ms)")

    for i in 1:length(fpms)
        rowsize!(f.layout, i, Relative(1/length(fpms)))
    end

    return f
end