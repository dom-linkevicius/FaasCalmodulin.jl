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