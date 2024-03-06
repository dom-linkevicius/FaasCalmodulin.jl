function gen_Stefan_TR_data(n)

    kf_AT_s = rand(Distributions.Uniform(1e2, 1e8), n)
    kb_AT_s = rand(Distributions.Uniform(1e-3, 1e5),n)

    kf_BT_s = rand(Distributions.Uniform(1e2, 1e8), n)
    kb_BT_s = rand(Distributions.Uniform(1e-3, 1e5), n)

    kf_CT_s = rand(Distributions.Uniform(1e2, 1e8), n)
    kb_CT_s = rand(Distributions.Uniform(1e-3, 1e5), n)

    kf_DT_s = rand(Distributions.Uniform(1e2, 1e8), n)
    kb_DT_s = rand(Distributions.Uniform(1e-3, 1e5), n)

    Ca_s    = rand([1.88e-6; 1.1e-6; 0.730e-6; 255e-6;
                    410e-6; 399e-6; 394e-6], n)
    CaM_T_s = rand([123e-6; 143e-6; 72e-6; 187e-6;
                    140e-6; 94e-6 ; 47e-6], n)
    
    all_inps = cat(kf_AT_s, kb_AT_s, kf_BT_s, kb_BT_s, kf_CT_s, kb_CT_s, kf_DT_s, kb_DT_s, Ca_s, CaM_T_s, dims=2)'
    all_outs = Stefan_TR_model_equilibriums(all_inps)

    return Flux.normalise(all_inps), all_outs
end


function train_model(in, out, batch_size, n_its, LR)

    model = Flux.f64(Flux.Chain(
        Flux.Dense(10 => 100, relu),
        Flux.Dense(100 => 100, relu),
        Flux.Dense(100 => 16, sigmoid)
    )) |> gpu

    loader = Flux.DataLoader((in, out) |>gpu, batchsize=batch_size, shuffle=true);
    
    optim = Flux.setup(Flux.Adam(LR), model) # will store optimiser momentum, etc.
    
    # Training loop, using the whole data set 1000 times:
    losses = []
    @showprogress for epoch in 1:n_its
        for (x, y) in loader
            loss, grads = Flux.withgradient(model) do m
                # Evaluate model and loss inside gradient context:
                y_hat = m(x)
                Flux.mae(y_hat, y)
            end
            Flux.update!(optim, model, grads[1])
            push!(losses, loss)  # logging, outside gradient context
        end
    end
    return losses, model
end