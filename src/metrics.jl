function rmse(fpm; pop=nothing)

    if pop === nothing
        preds = predict(fpm)
    else
        preds = predict(fpm, pop)
    end

    RMSEs = map(preds) do d
        if ismissing(d.subject.observations.F_F0[1])
            missing
        else
            rmse_i = sqrt.(mean((d.ipred.F_F0 .- d.subject.observations.F_F0).^2))
        end
   end

   return skipmissing(RMSEs)
end


function shifman_errors(fpms)

    all_errors = Vector{Float64}()

    idxs = .!occursin.("#", [i.id for i in fpms[1].data])

    for f in fpms
        p = predict(f)[idxs]

        errs_i = map(p) do j
            (j.subject.observations.Ca_CaM[1] - j.ipred.Ca_CaM[1])^2
        end
        push!(all_errors, sqrt(mean(errs_i)))
    end

    return all_errors
end


function my_aic(fpm::Pumas.AbstractFittedPumasModel; pop=nothing)
    if pop == nothing
        return 2*length(reduce(vcat, coef(fpm))) - 2*loglikelihood(fpm)
    else
        ll = loglikelihood(fpm.model, pop, coef(fpm), LaplaceI(), diffeq_options=(;alg=Rodas5P(), abstol=1e-16))
        return 2*length(reduce(vcat, coef(fpm))) - 2*ll
    end
end


function variance(x)
    x̄ = mean(x)
    n = length(x)

    s² = 1 / (n - 1) * sum((x .- x̄).^2)
end


function cohen_d(x1, x2)

    x̄₁  = mean(x1)
    x̄₂  = mean(x2)
    n₁  = length(x1)
    n₂  = length(x2)
    s²₁ = variance(x1)
    s²₂ = variance(x1)

    s = sqrt(((n₁ - 1)s²₁ + (n₂ - 1)s²₂)/(n₁ + n₂ - 2))

    d = abs((x̄₁ - x̄₂)/s)
end


function print_mses_training(nt)

    blackwell_mses_tr       = [mean(rmse(i)) for i in nt[:blackwell]];
    blackwell_fixed_mses_tr = [mean(rmse(i)) for i in nt[:blackwell_fixed]];
    bhalla_mses_tr          = [mean(rmse(i)) for i in nt[:bhalla]];
    bhalla_fixed_mses_tr    = [mean(rmse(i)) for i in nt[:bhalla_fixed]];
    shifman_mses_tr         = [mean(rmse(i)) for i in nt[:shifman]];
    shifman_fixed_mses_tr   = [mean(rmse(i)) for i in nt[:shifman_fixed]];
    pepke_m2_mses_tr        = [mean(rmse(i)) for i in nt[:pepke_m2]];
    pepke_m2_fixed_mses_tr  = [mean(rmse(i)) for i in nt[:pepke_m2_fixed]];
    faas_mses_tr            = [mean(rmse(i)) for i in nt[:faas]];
    faas_fixed_mses_tr      = [mean(rmse(i)) for i in nt[:faas_fixed]];
    pepke_m1_fixed_mses_tr  = [mean(rmse(i)) for i in nt[:pepke_m1_fixed]];
    byrne_mses_tr           = [mean(rmse(i)) for i in nt[:byrne]];
    byrne_fixed_mses_tr     = [mean(rmse(i)) for i in nt[:byrne_fixed]];

    println("Blackwell training RMSE: ", mean(blackwell_mses_tr), " STD: ", std(blackwell_mses_tr), "\n",
        "Blackwell Fixed training RMSE: ", mean(blackwell_fixed_mses_tr), " STD: ", std(blackwell_fixed_mses_tr), "\n",
         "Bhalla training RMSE: ", mean(bhalla_mses_tr), " STD: ", std(bhalla_mses_tr), "\n",
         "Bhalla Fixed training RMSE: ", mean(bhalla_fixed_mses_tr), " STD: ", std(bhalla_fixed_mses_tr), "\n",
         "Shifman training RMSE: ", mean(shifman_mses_tr), " STD: ", std(shifman_mses_tr), "\n",
         "Shifman Fixed training RMSE: ", mean(shifman_fixed_mses_tr), " STD: ", std(shifman_fixed_mses_tr), "\n",
         "Pepke M2 training RMSE: ", mean(pepke_m2_mses_tr), " STD: ", std(pepke_m2_mses_tr), "\n",
         "Pepke M2 Fixed training RMSE: ", mean(pepke_m2_fixed_mses_tr), " STD: ", std(pepke_m2_fixed_mses_tr), "\n",
         "Faas training RMSE: ", mean(faas_mses_tr), " STD: ", std(faas_mses_tr), "\n",
         "Faas Fixed training RMSE: ", mean(faas_fixed_mses_tr), " STD: ", std(faas_fixed_mses_tr), "\n",
         "Pepke M1 Fixed training RMSE: ", mean(pepke_m1_fixed_mses_tr), " STD: ", std(pepke_m1_fixed_mses_tr), "\n",
         "Byrne training RMSE: ", mean(byrne_mses_tr), " STD: ", std(byrne_mses_tr), "\n",
         "Byrne Fixed training RMSE: ", mean(byrne_fixed_mses_tr), " STD: ", std(byrne_fixed_mses_tr)
    )

    return (
        blackwell        = blackwell_mses_tr,
        blackwell_fixed  = blackwell_fixed_mses_tr,
        bhalla           = bhalla_mses_tr,
        bhalla_fixed     = bhalla_fixed_mses_tr,
        shifman          = shifman_mses_tr,
        shifman_fixed    = shifman_fixed_mses_tr,
        pepke_m2         = pepke_m2_mses_tr,
        pepke_m2_fixed   = pepke_m2_fixed_mses_tr,
        faas             = faas_mses_tr,
        faas_fixed       = faas_fixed_mses_tr,
        pepke_m1_fixed   = pepke_m1_fixed_mses_tr,
        byrne            = byrne_mses_tr,
        byrne_fixed      = byrne_fixed_mses_tr,
    )
end


function print_mses_validation(nt, faas_pop, val_idxs)

    blackwell_mses_val       = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:blackwell], val_idxs)];
    blackwell_fixed_mses_val = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:blackwell_fixed], val_idxs)];
    bhalla_mses_val          = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:bhalla], val_idxs)];
    bhalla_fixed_mses_val    = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:bhalla_fixed], val_idxs)];
    shifman_mses_val         = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:shifman], val_idxs)];
    shifman_fixed_mses_val   = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:shifman_fixed], val_idxs)];
    pepke_m2_mses_val        = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:pepke_m2], val_idxs)];
    pepke_m2_fixed_mses_val  = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:pepke_m2_fixed], val_idxs)];
    faas_mses_val            = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:faas], val_idxs)];
    faas_fixed_mses_val      = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:faas_fixed], val_idxs)];
    pepke_m1_fixed_mses_val  = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:pepke_m1_fixed], val_idxs)];
    byrne_mses_val           = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:byrne], val_idxs)];
    byrne_fixed_mses_val     = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:byrne_fixed], val_idxs)];

    println("Blackwell validation RMSE: ", mean(blackwell_mses_val), " STD: ", std(blackwell_mses_val), "\n",
        "Blackwell Fixed validation RMSE: ", mean(blackwell_fixed_mses_val), " STD: ", std(blackwell_fixed_mses_val), "\n",
         "Bhalla validation RMSE: ", mean(bhalla_mses_val), " STD: ", std(bhalla_mses_val), "\n",
         "Bhalla Fixed validation RMSE: ", mean(bhalla_fixed_mses_val), " STD: ", std(bhalla_fixed_mses_val), "\n",
         "Shifman validation RMSE: ", mean(shifman_mses_val), " STD: ", std(shifman_mses_val), "\n",
         "Shifman Fixed validation RMSE: ", mean(shifman_fixed_mses_val), " STD: ", std(shifman_fixed_mses_val), "\n",
         "Pepke M2 validation RMSE: ", mean(pepke_m2_mses_val), " STD: ", std(pepke_m2_mses_val), "\n",
         "Pepke M2 Fixed validation RMSE: ", mean(pepke_m2_fixed_mses_val), " STD: ", std(pepke_m2_fixed_mses_val), "\n",
         "Faas validation RMSE: ", mean(faas_mses_val), " STD: ", std(faas_mses_val), "\n",
         "Faas Fixed validation RMSE: ", mean(faas_fixed_mses_val), " STD: ", std(faas_fixed_mses_val), "\n",
         "Pepke M1 Fixed validation RMSE: ", mean(pepke_m1_fixed_mses_val), " STD: ", std(pepke_m1_fixed_mses_val), "\n",
         "Byrne validation RMSE: ", mean(byrne_mses_val), " STD: ", std(byrne_mses_val), "\n",
         "Byrne Fixed validation RMSE: ", mean(byrne_fixed_mses_val), " STD: ", std(byrne_fixed_mses_val)
    )

    return (
        blackwell        = blackwell_mses_val,
        blackwell_fixed  = blackwell_fixed_mses_val,
        bhalla           = bhalla_mses_val,
        bhalla_fixed     = bhalla_fixed_mses_val,
        shifman          = shifman_mses_val,
        shifman_fixed    = shifman_fixed_mses_val,
        pepke_m2         = pepke_m2_mses_val,
        pepke_m2_fixed   = pepke_m2_fixed_mses_val,
        faas             = faas_mses_val,
        faas_fixed       = faas_fixed_mses_val,
        pepke_m1_fixed   = pepke_m1_fixed_mses_val,
        byrne            = byrne_mses_val,
        byrne_fixed      = byrne_fixed_mses_val,
    )
end


function print_mses_testing(nt, faas_pop, test_idxs)

    blackwell_mses_test       = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:blackwell], test_idxs)];
    blackwell_fixed_mses_test = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:blackwell_fixed], test_idxs)];
    bhalla_mses_test          = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:bhalla], test_idxs)];
    bhalla_fixed_mses_test    = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:bhalla_fixed], test_idxs)];
    shifman_mses_test         = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:shifman], test_idxs)];
    shifman_fixed_mses_test   = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:shifman_fixed], test_idxs)];
    pepke_m2_mses_test        = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:pepke_m2], test_idxs)];
    pepke_m2_fixed_mses_test  = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:pepke_m2_fixed], test_idxs)];
    faas_mses_test            = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:faas], test_idxs)];
    faas_fixed_mses_test      = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:faas_fixed], test_idxs)];
    pepke_m1_fixed_mses_test  = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:pepke_m1_fixed], test_idxs)];
    byrne_mses_test           = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:byrne], test_idxs)];
    byrne_fixed_mses_test     = [mean(rmse(i; pop = faas_pop[j])) for (i,j) in zip(nt[:byrne_fixed], test_idxs)];

    println("Blackwell testing RMSE: ", mean(blackwell_mses_test), " STD: ", std(blackwell_mses_test), "\n",
        "Blackwell Fixed testing RMSE: ", mean(blackwell_fixed_mses_test), " STD: ", std(blackwell_fixed_mses_test), "\n",
         "Bhalla testing RMSE: ", mean(bhalla_mses_test), " STD: ", std(bhalla_mses_test), "\n",
         "Bhalla Fixed testing RMSE: ", mean(bhalla_fixed_mses_test), " STD: ", std(bhalla_fixed_mses_test), "\n",
         "Shifman testing RMSE: ", mean(shifman_mses_test), " STD: ", std(shifman_mses_test), "\n",
         "Shifman Fixed testing RMSE: ", mean(shifman_fixed_mses_test), " STD: ", std(shifman_fixed_mses_test), "\n",
         "Pepke M2 testing RMSE: ", mean(pepke_m2_mses_test), " STD: ", std(pepke_m2_mses_test), "\n",
         "Pepke M2 Fixed testing RMSE: ", mean(pepke_m2_fixed_mses_test), " STD: ", std(pepke_m2_fixed_mses_test), "\n",
         "Faas testing RMSE: ", mean(faas_mses_test), " STD: ", std(faas_mses_test), "\n",
         "Faas Fixed testing RMSE: ", mean(faas_fixed_mses_test), " STD: ", std(faas_fixed_mses_test), "\n",
         "Pepke M1 Fixed testing RMSE: ", mean(pepke_m1_fixed_mses_test), " STD: ", std(pepke_m1_fixed_mses_test), "\n",
         "Byrne testing RMSE: ", mean(byrne_mses_test), " STD: ", std(byrne_mses_test), "\n",
         "Byrne Fixed testing RMSE: ", mean(byrne_fixed_mses_test), " STD: ", std(byrne_fixed_mses_test)
    )

    return (
        blackwell        = blackwell_mses_test,
        blackwell_fixed  = blackwell_fixed_mses_test,
        bhalla           = bhalla_mses_test,
        bhalla_fixed     = bhalla_fixed_mses_test,
        shifman          = shifman_mses_test,
        shifman_fixed    = shifman_fixed_mses_test,
        pepke_m2         = pepke_m2_mses_test,
        pepke_m2_fixed   = pepke_m2_fixed_mses_test,
        faas             = faas_mses_test,
        faas_fixed       = faas_fixed_mses_test,
        pepke_m1_fixed   = pepke_m1_fixed_mses_test,
        byrne            = byrne_mses_test,
        byrne_fixed      = byrne_fixed_mses_test,
    )
end