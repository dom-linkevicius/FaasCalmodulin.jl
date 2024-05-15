using Pumas
using DeepPumas
using CairoMakie
using Plots

using XLSX
using DifferentialEquations
using Random
using JLD2
using Symbolics
using LaTeXStrings
using Serialization
using HypothesisTests
using NumericalIntegration


include("models.jl");
include("SymbolicsForwardDiffExt.jl")
include("process_data.jl");
include("convenience.jl");
include("symbolic_manipulations.jl");
include("eq_comparisons.jl");
include("dyn_comparisons.jl");

f_xlsx = readxlsx("data/dataset.xlsx");
full_pop     = fetch_faas(f_xlsx);

f_shifman = XLSX.readxlsx("data/Shifman_2006.xlsx");
shifman_pop  = fetch_shifman(f_shifman);

##f_peersen   = XLSX.readxlsx("Peersen_1997.xlsx")
###peersen_pop = fetch_peersen(f_peersen);

A = full_pop[1:13];
B = full_pop[14:24];
C = full_pop[25:35];
D = full_pop[36:51];
E = full_pop[52:65];
F = full_pop[66:80];
G = full_pop[81:94];


n_models = 1#20
seeds = 1:n_models
sub_range = 1:20:201


fpms_bwell, p_bwell_train, p_bwell_val, val_idxs_bwell, test_idxs_bwell = 
train_plot("Blackwell", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds;
    sub_range = sub_range,
);
###serialize("Data/fpms_bwell_shifman_all_data_subsample.jls",fpms_bwell)
###serialize("Data/val_idxs_bwell_shifman_all_data_subsample.jls",val_idxs_bwell)
###serialize("Data/test_idxs_bwell_shifman_all_data_subsample.jls",test_idxs_bwell)

###fpms_bwell = deserialize("Data/fpms_bwell_shifman_all_data_subsample.jls")
###val_idxs_bwell = deserialize("Data/val_idxs_bwell_shifman_all_data_subsample.jls")
###test_idxs_bwell = deserialize("Data/test_idxs_bwell_shifman_all_data_subsample.jls")
###fpms_bwell = fpms_bwell[[1:11;13:20]]

fpms_bwell_fixed, p_bwell_train_fixed, p_bwell_val_fixed, val_idxs_bwell_fixed, test_idxs_bwell_fixed = 
train_plot("Blackwell_fixed", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_bwell_fixed_shifman_all_data_subsample.jls",fpms_bwell_fixed)
###fpms_bwell_fixed = deserialize("Data/fpms_bwell_fixed_shifman_all_data_subsample.jls")



fpms_bhalla, p_bhalla_train, p_bhalla_val, val_idxs_bhalla, test_idxs_bhalla = 
train_plot("Bhalla", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds;
    sub_range = sub_range,
);
###serialize("Data/fpms_bhalla_shifman_all_data_subsample.jls",fpms_bhalla)
###fpms_bhalla = deserialize("Data/fpms_bhalla_shifman_all_data_subsample.jls")

fpms_bhalla_fixed, p_bhalla_train_fixed, p_bhalla_val_fixed, val_idxs_bhalla_fixed, test_idxs_bhalla_fixed = 
train_plot("Bhalla_fixed", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds;
    sub_range = sub_range,
);
###serialize("Data/fpms_bhalla_fixed_shifman_all_data_subsample.jls",fpms_bhalla_fixed)
###fpms_bhalla_fixed = deserialize("Data/fpms_bhalla_fixed_shifman_all_data_subsample.jls")



fpms_shifman, p_shifman_train, p_shifman_val, val_idxs_shifman, test_idxs_shifman = 
train_plot("Shifman", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds;
    sub_range = sub_range,
);
###serialize("Data/fpms_shifman_shifman_all_data_subsample.jls",fpms_shifman)
###fpms_shifman = deserialize("Data/fpms_shifman_shifman_all_data_subsample.jls")

fpms_shifman_fixed, p_shifman_fixed_train, p_shifman_fixed_val, val_idxs_shifman_fixed, test_idxs_shifman_fixed = 
train_plot("Shifman_fixed", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds;
    sub_range = sub_range,
);
###serialize("Data/fpms_shifman_fixed_shifman_all_data_subsample.jls",fpms_shifman_fixed)
###fpms_shifman_fixed = deserialize("Data/fpms_shifman_fixed_shifman_all_data_subsample.jls")



fpms_pepke_m2, p_pepke_m2_train, p_pepke_m2_val, val_idxs_pepke_m2, test_idxs_pepke_m2 = 
train_plot("Pepke_m2", 
    full_pop, 
    shifman_pop,
    n_models, 
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_pepke_m2_shifman_all_data_subsample.jls",fpms_pepke_m2)
###fpms_pepke_m2 = deserialize("Data/fpms_pepke_m2_shifman_all_data_subsample.jls")

fpms_pepke_m2_fixed, p_pepke_m2_train_fixed, p_pepke_m2_val_fixed, val_idxs_pepke_m2_fixed, test_idxs_pepke_m2_fixed = 
train_plot("Pepke_m2_fixed", 
    full_pop, 
    shifman_pop,
    n_models, 
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_pepke_m2_fixed_shifman_all_data_subsample.jls",fpms_pepke_m2_fixed)
###fpms_pepke_m2_fixed = deserialize("Data/fpms_pepke_m2_fixed_shifman_all_data_subsample.jls")



fpms_faas, p_faas_train, p_faas_val, val_idxs_faas, test_idxs_faas = 
train_plot("Faas", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds;
    sub_range = sub_range
);
###serialize("Data/fpms_faas_shifman_all_data_subsample.jls",fpms_faas)
###fpms_faas = deserialize("Data/fpms_faas_shifman_all_data_subsample.jls")

fpms_faas_fixed, p_faas_train_fixed, p_faas_val_fixed, val_idxs_faas_fixed, test_idxs_faas_fixed = 
train_plot("Faas_fixed", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_faas_fixed_shifman_all_data_subsample.jls",fpms_faas_fixed)
###fpms_faas_fixed = deserialize("Data/fpms_faas_fixed_shifman_all_data_subsample.jls")

fpms_pepke_m1_fixed, p_pepke_m1_train_fixed, p_pepke_m1_val_fixed, val_idxs_pepke_m1_fixed, test_idxs_pepke_m1_fixed = 
train_plot("Pepke_m1_fixed", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_pepke_m1_fixed_shifman_all_data_subsample.jls",fpms_pepke_m1_fixed)
###fpms_pepke_m1_fixed = deserialize("Data/fpms_pepke_m1_fixed_shifman_all_data_subsample.jls")



fpms_byrne, p_byrne_train, p_byrne_val, val_idxs_byrne, test_idxs_byrne = 
train_plot("Byrne", 
    full_pop,
    shifman_pop,
    n_models,
    seeds; 
    sub_range = sub_range,
);
###serialize("Data/fpms_byrne_shifman_all_data_subsample.jls",fpms_byrne)
###fpms_byrne = deserialize("Data/fpms_byrne_shifman_all_data_subsample.jls")

fpms_byrne_fixed, p_byrne_train_fixed, p_byrne_val_fixed, val_idxs_byrne_fixed, test_idxs_byrne_fixed = 
train_plot("Byrne_fixed", 
    full_pop, 
    shifman_pop,
    n_models,
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_byrne_fixed_shifman_all_data_subsample.jls",fpms_byrne_fixed)
###fpms_byrne_fixed = deserialize("Data/fpms_byrne_fixed_shifman_all_data_subsample.jls")


fpms_bwell          = deserialize("Data/fpms_bwell_shifman_all_data_subsample.jls");
val_idxs_bwell      = deserialize("Data/val_idxs_bwell_shifman_all_data_subsample.jls");
test_idxs_bwell     = deserialize("Data/test_idxs_bwell_shifman_all_data_subsample.jls");
fpms_bwell          = fpms_bwell[[1:11;13:20]];
fpms_bwell_fixed    = deserialize("Data/fpms_bwell_fixed_shifman_all_data_subsample.jls");
fpms_bhalla         = deserialize("Data/fpms_bhalla_shifman_all_data_subsample.jls");
fpms_bhalla_fixed   = deserialize("Data/fpms_bhalla_fixed_shifman_all_data_subsample.jls");
fpms_shifman        = deserialize("Data/fpms_shifman_shifman_all_data_subsample.jls");
fpms_shifman_fixed  = deserialize("Data/fpms_shifman_fixed_shifman_all_data_subsample.jls");
fpms_pepke_m2       = deserialize("Data/fpms_pepke_m2_shifman_all_data_subsample.jls");
fpms_pepke_m2_fixed = deserialize("Data/fpms_pepke_m2_fixed_shifman_all_data_subsample.jls");
fpms_faas           = deserialize("Data/fpms_faas_shifman_all_data_subsample.jls");
fpms_faas_fixed     = deserialize("Data/fpms_faas_fixed_shifman_all_data_subsample.jls");
fpms_pepke_m1_fixed = deserialize("Data/fpms_pepke_m1_fixed_shifman_all_data_subsample.jls");
fpms_byrne          = deserialize("Data/fpms_byrne_shifman_all_data_subsample.jls");
fpms_byrne_fixed    = deserialize("Data/fpms_byrne_fixed_shifman_all_data_subsample.jls");



###nn = 1
###sims_faas     = simobs(fpms_faas_fixed[1])
###sims_pepke_m1 = simobs(fpms_pepke_m1_fixed[1])
###sims_pepke_m2 = simobs(fpms_pepke_m2_fixed[1])
###
###preds_faas     = predict(fpms_faas_fixed[1])
###preds_pepke_m1 = predict(fpms_pepke_m1_fixed[1])
###preds_pepke_m2 = predict(fpms_pepke_m2_fixed[1])
###
###fg = CairoMakie.Figure()
###
###ax1 = CairoMakie.Axis(fg[2, 1], title="CaM C0")
###l1 = CairoMakie.lines!(ax1, sims_faas[nn].time, sims_faas[nn].sol[end-2,:], label="Faas", color=our_col)
###l2 = CairoMakie.lines!(ax1, sims_pepke_m1[nn].time, sims_pepke_m1[nn].sol[end-2,:], label="Pepke M1", color=kim_col)
###l3 = CairoMakie.lines!(ax1, sims_pepke_m2[nn].time, sims_pepke_m2[nn].sol[end-3,:], label="Pepke M2", color=faas_col)
###hidedecorations!(ax1)
###
###ax2 = CairoMakie.Axis(fg[3, 1], title="CaM C1")
###CairoMakie.lines!(ax2, sims_faas[nn].time, sims_faas[nn].sol[end-1,:], color=our_col)
###CairoMakie.lines!(ax2, sims_pepke_m1[nn].time, sims_pepke_m1[nn].sol[end-1,:], color=kim_col)
###hidedecorations!(ax2)
###
###ax3 = CairoMakie.Axis(fg[4, 1], title="CaM C2")
###CairoMakie.lines!(ax3, sims_faas[nn].time, sims_faas[nn].sol[end-0,:], color=our_col)
###CairoMakie.lines!(ax3, sims_pepke_m1[nn].time, sims_pepke_m1[nn].sol[end-0,:], color=kim_col)
###CairoMakie.lines!(ax3, sims_pepke_m2[nn].time, sims_pepke_m2[nn].sol[end-2,:], color=faas_col)
###hidedecorations!(ax3)
###
###ax4 = CairoMakie.Axis(fg[2, 2], title="CaM N0")
###CairoMakie.lines!(ax4, sims_faas[nn].time, sims_faas[nn].sol[end-5,:], color=our_col)
###CairoMakie.lines!(ax4, sims_pepke_m1[nn].time, sims_pepke_m1[nn].sol[end-5,:], color=kim_col)
###CairoMakie.lines!(ax4, sims_pepke_m2[nn].time, sims_pepke_m2[nn].sol[end-3,:], color=faas_col)
###hidedecorations!(ax4)
###
###ax5 = CairoMakie.Axis(fg[3, 2], title="CaM N1")
###CairoMakie.lines!(ax5, sims_faas[nn].time, sims_faas[nn].sol[end-4,:], color=our_col)
###CairoMakie.lines!(ax5, sims_pepke_m1[nn].time, sims_pepke_m1[nn].sol[end-4,:], color=kim_col)
###hidedecorations!(ax5)
###
###ax6 = CairoMakie.Axis(fg[4, 2], title="CaM N2")
###CairoMakie.lines!(ax6, sims_faas[nn].time, sims_faas[nn].sol[end-3,:], color=our_col)
###CairoMakie.lines!(ax6, sims_pepke_m1[nn].time, sims_pepke_m1[nn].sol[end-3,:], color=kim_col)
###CairoMakie.lines!(ax6, sims_pepke_m2[nn].time, sims_pepke_m2[nn].sol[end-1,:], color=faas_col)
###hidedecorations!(ax6)
###
###ax7 = CairoMakie.Axis(fg[5, 1:2], title="CaM 4")
###CairoMakie.lines!(ax7, sims_pepke_m2[nn].time, sims_pepke_m2[nn].sol[end,:], color=faas_col)
###hidedecorations!(ax7)
###
###CairoMakie.Legend(fg[1, :], [l1, l2, l3], ["Faas", "Pepke M1", "Pepke M2"], nbanks=3) 
###
###ax8 = CairoMakie.Axis(fg[:, 3:4], title="ΔF/F₀")
###
###faas_pred = (sims_faas[nn].dynamics.OGB5 + 39.364 * sims_faas[nn].dynamics.CaOGB5) / (sims_faas[nn].icoefs.OGB5₀[1] + 39.364 * sims_faas[nn].icoefs.CaOGB5₀[1])
###pepke_m1_pred = (sims_pepke_m1[nn].dynamics.OGB5 + 39.364 * sims_pepke_m1[nn].dynamics.CaOGB5) / (sims_pepke_m1[nn].icoefs.OGB5₀[1] + 39.364 * sims_pepke_m1[nn].icoefs.CaOGB5₀[1])
###pepke_m2_pred = (sims_pepke_m2[nn].dynamics.OGB5 + 39.364 * sims_pepke_m2[nn].dynamics.CaOGB5) / (sims_pepke_m2[nn].icoefs.OGB5₀[1] + 39.364 * sims_pepke_m2[nn].icoefs.CaOGB5₀[1])
###
###CairoMakie.lines!(ax8, sims_faas[nn].time, sims_faas[nn].subject.observations.F_F0, color=:black)
###CairoMakie.lines!(ax8, sims_faas[nn].time,      faas_pred, color=our_col)
###CairoMakie.lines!(ax8, sims_pepke_m1[nn].time,  pepke_m1_pred, color=kim_col)
###CairoMakie.lines!(ax8, sims_pepke_m2[nn].time,  pepke_m2_pred, color=faas_col)


######full_pop[1].observations.F_F0
######
###fpm_ms = msfit(
###    Faas_scheme,
###    fpms_bwell[1].data,
###    sample_params(Faas_scheme),
###    ###MAP(FOCE()),
###    JointMAP(),
###    2,
###    5e-16, #or standard deviation of shooting penalty hyper prior if below argument is not nothing
###    :L1; #:L1 or :L2, :L1 should be highly favoured because it performed much better in experiments
###    hyper_prior_type = nothing, #nothing/:L1/:L2,
###    return_with_shooting = true, #true/false, if false, does a 0 iteration fit and returns a FittePumasModel with the original model
###    #regular Pumas.fit kwargs, e.g. optim_options, diffeq_options, constantcoef, etc.
###    optim_options=(; iterations=10, store_trace=true),
###    diffeq_options=(;alg=Rodas5P(), abstol=1e-16)
###)

############################
### GLOBAL PLOT SETTINGS ###
############################

regular_font_string = "Fonts/Arial Regular.ttf"
bold_font_string = "Fonts/Arial Bold.ttf"
set_theme!(fonts=(;regular = regular_font_string, bold= bold_font_string))
dpc = 236.22 ### dots per cm, this is equal to 600dpi
pt_per_cm = 28.3465

our_col     = RGB( 34/255,113/255,178/255);
bhalla_col  = RGB(061/255,183/255,233/255);
kim_col     = RGB(247/255,072/255,165/255);
pepke_col   = RGB(053/255,155/255,115/255);
faas_col    = RGB(213/255,094/255,000/255);
shifman_col = RGB(230/255,159/255,000/255);
byrne_col   = RGB(240/255,228/255,066/255);

###########################
### AIC BOX PLOTS START ###
###########################

all_fpms = [fpms_bwell;    fpms_bwell_fixed;
            fpms_bhalla;    fpms_bhalla_fixed;
            fpms_shifman;    fpms_shifman_fixed;
            fpms_pepke_m2; fpms_pepke_m2_fixed;
            fpms_faas;     fpms_faas_fixed; fpms_pepke_m1_fixed;
            fpms_byrne;    fpms_byrne_fixed];


val_pops = [idx_to_val_vec(full_pop, idxs) for idxs in val_idxs_bwell];
test_pops = [idx_to_val_vec(full_pop, idxs) for idxs in test_idxs_bwell];
all_val_pops = repeat(val_pops, 13);
all_test_pops = repeat(test_pops, 13);


all_names = [repeat(["Scheme 1 + our rates"], length(fpms_bwell)); 
             repeat(["Scheme 1 + Kim et al."], length(fpms_bwell_fixed));
             repeat(["Scheme 2 + our rates"], length(fpms_pepke_m2));
             repeat(["Scheme 2 + Hayer and Bhalla"], length(fpms_pepke_m2_fixed));
             repeat(["Scheme 3 + our rates"], length(fpms_pepke_m2));
             repeat(["Scheme 3 + Shifman et al. + our rates"], length(fpms_pepke_m2_fixed));
             repeat(["Scheme 4 + our rates"], length(fpms_pepke_m2));
             repeat(["Scheme 4 + Pepke et al."], length(fpms_pepke_m2_fixed));
             repeat(["Scheme 5 + our rates"], length(fpms_faas));
             repeat(["Scheme 5 + Faas et al."], length(fpms_faas_fixed));
             repeat(["Scheme 5 + Pepke et al."], length(fpms_pepke_m1_fixed));
             repeat(["Scheme 6 + our rates"], length(fpms_byrne));
             repeat(["Scheme 6 + Byrne et al."], length(fpms_byrne_fixed))];

aic_pairs = zip(all_fpms, all_test_pops, all_names);

x, y = aic_boxplots(aic_pairs);
###x_makie = sort(repeat(1:length(unique(x)), Int(length(x)/length(unique(x)))));
x_makie = [
    repeat([1],  length(fpms_bwell));
    repeat([2],  length(fpms_bwell_fixed));
    repeat([3],  length(fpms_bhalla));
    repeat([4],  length(fpms_bhalla_fixed));
    repeat([5],  length(fpms_shifman));
    repeat([6],  length(fpms_shifman_fixed));
    repeat([7],  length(fpms_pepke_m2));
    repeat([8],  length(fpms_pepke_m2_fixed));
    repeat([9],  length(fpms_faas))
    repeat([10], length(fpms_faas_fixed))
    repeat([11], length(fpms_pepke_m1_fixed))
    repeat([12], length(fpms_byrne))
    repeat([13], length(fpms_byrne_fixed))
]
y_makie = [i for i in y];
boxplot_cols = [repeat([our_col],       length(fpms_bwell)); 
                repeat([kim_col],       length(fpms_bwell_fixed));
                repeat([our_col],       length(fpms_bhalla));
                repeat([bhalla_col],    length(fpms_bhalla_fixed));
                repeat([our_col],       length(fpms_shifman));
                repeat([shifman_col],   length(fpms_shifman_fixed));
                repeat([our_col],       length(fpms_pepke_m2));
                repeat([pepke_col],     length(fpms_pepke_m2_fixed));
                repeat([our_col],       length(fpms_faas));
                repeat([faas_col],      length(fpms_faas_fixed));
                repeat([pepke_col],     length(fpms_pepke_m1_fixed));
                repeat([our_col],       length(fpms_byrne));
                repeat([byrne_col],     length(fpms_byrne_fixed))];

boxplot_size_cm     = (17, 7.5)
boxplot_size_px     =  boxplot_size_cm .* dpc 
boxplot_size_units  =  boxplot_size_cm .* pt_per_cm
boxplot_px_per_unit = (boxplot_size_cm .* dpc ./ (boxplot_size_cm .* pt_per_cm))[1]

boxplot = CairoMakie.Figure(size=boxplot_size_units, fontsize=8)
ax_boxplot = CairoMakie.Axis(boxplot[1,1], 
    xticks=(unique(x_makie), unique(x)),
    xlabel="Scheme + parameters",
    ylabel="log10(AIC)",
    xticklabelrotation=pi/9,
    xticklabelpad=5,
    xlabelpadding=15,
    xticklabelsize=8,
    yticklabelsize=8,
    xgridvisible=false,
    yscale=log10
    )
CairoMakie.boxplot!(ax_boxplot, x_makie, y_makie, show_notch=false, color=boxplot_cols)
###CairoMakie.save("Plots/AIC_boxplot.svg", boxplot, pt_per_unit=1, 
###    px_per_unit = boxplot_px_per_unit, size=boxplot_size_units)
CairoMakie.save("Plots/AIC_boxplot.png", boxplot, pt_per_unit=1, 
    px_per_unit = boxplot_px_per_unit, size=boxplot_size_units)   

#########################
### AIC BOX PLOTS END ###
#########################
#############################
### MSE HISTS PLOTS START ###
#############################

bwell_mses_tr           = [mean(rmse(i)) for i in fpms_bwell];
bwell_fixed_mses_tr     = [mean(rmse(i)) for i in fpms_bwell_fixed];

bhalla_mses_tr          = [mean(rmse(i)) for i in fpms_bhalla];
bhalla_fixed_mses_tr    = [mean(rmse(i)) for i in fpms_bhalla_fixed];

shifman_mses_tr         = [mean(rmse(i)) for i in fpms_shifman];
shifman_fixed_mses_tr   = [mean(rmse(i)) for i in fpms_shifman_fixed];

pepke_m2_mses_tr        = [mean(rmse(i)) for i in fpms_pepke_m2];
pepke_m2_fixed_mses_tr  = [mean(rmse(i)) for i in fpms_pepke_m2_fixed];

faas_mses_tr            = [mean(rmse(i)) for i in fpms_faas];
faas_fixed_mses_tr      = [mean(rmse(i)) for i in fpms_faas_fixed];
pepke_m1_fixed_mses_tr  = [mean(rmse(i)) for i in fpms_pepke_m1_fixed];

byrne_mses_tr           = [mean(rmse(i)) for i in fpms_byrne];
byrne_fixed_mses_tr     = [mean(rmse(i)) for i in fpms_byrne_fixed];

println("Blackwell training RMSE: ", mean(bwell_mses_tr), " STD: ", std(bwell_mses_tr), "\n",
"Blackwell Fixed training RMSE: ", mean(bwell_fixed_mses_tr), " STD: ", std(bwell_fixed_mses_tr), "\n",
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

bwell_mses           = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_bwell, val_idxs_bwell)];
bwell_fixed_mses     = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_bwell_fixed, val_idxs_bwell)];

bhalla_mses          = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_bhalla, val_idxs_bwell)];
bhalla_fixed_mses    = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_bhalla_fixed, val_idxs_bwell)];

shifman_mses         = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_shifman, val_idxs_bwell)];
shifman_fixed_mses   = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_shifman_fixed, val_idxs_bwell)];

pepke_m2_mses        = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_pepke_m2, val_idxs_bwell)];
pepke_m2_fixed_mses  = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_pepke_m2_fixed, val_idxs_bwell)];

faas_mses            = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_faas, val_idxs_bwell)];
faas_fixed_mses      = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_faas_fixed, val_idxs_bwell)];
pepke_m1_fixed_mses  = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_pepke_m1_fixed, val_idxs_bwell)];

byrne_mses           = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_byrne, val_idxs_bwell)];
byrne_fixed_mses     = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_byrne_fixed, val_idxs_bwell)];

println("Blackwell validation RMSE: ", mean(bwell_mses), " STD: ", std(bwell_mses), "\n",
"Blackwell Fixed validation RMSE: ", mean(bwell_fixed_mses), " STD: ", std(bwell_fixed_mses), "\n",
"Bhalla validation RMSE: ", mean(bhalla_mses), " STD: ", std(bhalla_mses), "\n",
"Bhalla Fixed validation RMSE: ", mean(bhalla_fixed_mses), " STD: ", std(bhalla_fixed_mses), "\n",
"Shifman validation RMSE: ", mean(shifman_mses), " STD: ", std(shifman_mses), "\n",
"Shifman Fixed validation RMSE: ", mean(shifman_fixed_mses), " STD: ", std(shifman_fixed_mses), "\n",
"Pepke M2 validation RMSE: ", mean(pepke_m2_mses), " STD: ", std(pepke_m2_mses), "\n",
"Pepke M2 Fixed validation RMSE: ", mean(pepke_m2_fixed_mses), " STD: ", std(pepke_m2_fixed_mses), "\n",
"Faas validation RMSE: ", mean(faas_mses), " STD: ", std(faas_mses), "\n",
"Faas Fixed validation RMSE: ", mean(faas_fixed_mses), " STD: ", std(faas_fixed_mses), "\n",
"Pepke M1 Fixed validation RMSE: ", mean(pepke_m1_fixed_mses), " STD: ", std(pepke_m1_fixed_mses), "\n",
"Byrne validation RMSE: ", mean(byrne_mses), " STD: ", std(byrne_mses), "\n",
"Byrne Fixed validation RMSE: ", mean(byrne_fixed_mses), " STD: ", std(byrne_fixed_mses))

bwell_mses_test           = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_bwell, test_idxs_bwell)];
bwell_fixed_mses_test     = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_bwell_fixed, test_idxs_bwell)];

bhalla_mses_test          = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_bhalla, test_idxs_bwell)];
bhalla_fixed_mses_test    = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_bhalla_fixed, test_idxs_bwell)];

shifman_mses_test         = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_shifman, test_idxs_bwell)];
shifman_fixed_mses_test   = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_shifman_fixed, test_idxs_bwell)];

pepke_m2_mses_test        = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_pepke_m2, test_idxs_bwell)];
pepke_m2_fixed_mses_test  = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_pepke_m2_fixed, test_idxs_bwell)];

faas_mses_test            = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_faas, test_idxs_bwell)];
faas_fixed_mses_test      = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_faas_fixed, test_idxs_bwell)];
pepke_m1_fixed_mses_test  = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_pepke_m1_fixed, test_idxs_bwell)];

byrne_mses_test           = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_byrne, test_idxs_bwell)];
byrne_fixed_mses_test     = [mean(rmse(i; pop = full_pop[j])) for (i,j) in zip(fpms_byrne_fixed, test_idxs_bwell)];

println("Blackwell test RMSE: ", mean(bwell_mses_test), " STD: ", std(bwell_mses_test), "\n",
"Blackwell Fixed test RMSE: ", mean(bwell_fixed_mses_test), " STD: ", std(bwell_fixed_mses_test), "\n",
"Bhalla test RMSE: ", mean(bhalla_mses_test), " STD: ", std(bhalla_mses_test), "\n",
"Bhalla Fixed test RMSE: ", mean(bhalla_fixed_mses_test), " STD: ", std(bhalla_fixed_mses_test), "\n",
"Shifman test RMSE: ", mean(shifman_mses_test), " STD: ", std(shifman_mses_test), "\n",
"Shifman Fixed test RMSE: ", mean(shifman_fixed_mses_test), " STD: ", std(shifman_fixed_mses_test), "\n",
"Pepke M2 test RMSE: ", mean(pepke_m2_mses_test), " STD: ", std(pepke_m2_mses_test), "\n",
"Pepke M2 Fixed test RMSE: ", mean(pepke_m2_fixed_mses_test), " STD: ", std(pepke_m2_fixed_mses_test), "\n",
"Faas test RMSE: ", mean(faas_mses_test), " STD: ", std(faas_mses_test), "\n",
"Faas Fixed test RMSE: ", mean(faas_fixed_mses_test), " STD: ", std(faas_fixed_mses_test), "\n",
"Pepke M1 Fixed test RMSE: ", mean(pepke_m1_fixed_mses_test), " STD: ", std(pepke_m1_fixed_mses_test), "\n",
"Byrne test RMSE: ", mean(byrne_mses_test), " STD: ", std(byrne_mses_test), "\n",
"Byrne Fixed test RMSE: ", mean(byrne_fixed_mses_test), " STD: ", std(byrne_fixed_mses_test))

cohen_d_bwell    = cohen_d(bwell_mses_test, bwell_fixed_mses_test);
cohen_d_bhalla   = cohen_d(bhalla_mses_test, bhalla_fixed_mses_test);
cohen_d_shifman  = cohen_d(shifman_mses_test, shifman_fixed_mses_test);
cohen_d_pepke_m2 = cohen_d(pepke_m2_mses_test, pepke_m2_fixed_mses_test);
cohen_d_faas     = cohen_d(faas_mses_test, faas_fixed_mses_test);
cohen_d_pepke_m1 = cohen_d(faas_mses_test, pepke_m1_fixed_mses_test);
cohen_d_byrne    = cohen_d(byrne_mses_test, byrne_fixed_mses_test);

println("Blackwell Cohen's d: ", cohen_d_bwell, "\n",
"Bhalla Cohen's d: ", cohen_d_bhalla, "\n",
"Shifman Cohen's d: ", cohen_d_shifman, "\n",
"Pepke M2 Cohen's d: ", cohen_d_pepke_m2, "\n",
"Faas Fixed Cohen's d: ", cohen_d_faas, "\n",
"Pepke M1 Fixed Cohen's d: ", cohen_d_pepke_m1, "\n",
"Byrne Fixed Cohen's d: ", cohen_d_byrne)


bwell_mses_eq           = shifman_errors(fpms_bwell);
bwell_fixed_mses_eq     = shifman_errors(fpms_bwell_fixed);

bhalla_mses_eq          = shifman_errors(fpms_bhalla);
bhalla_fixed_mses_eq    = shifman_errors(fpms_bhalla_fixed);

shifman_mses_eq         = shifman_errors(fpms_shifman);
shifman_fixed_mses_eq   = shifman_errors(fpms_shifman_fixed);

pepke_m2_mses_eq        = shifman_errors(fpms_pepke_m2);
pepke_m2_fixed_mses_eq  = shifman_errors(fpms_pepke_m2_fixed);

faas_mses_eq            = shifman_errors(fpms_faas);
faas_fixed_mses_eq      = shifman_errors(fpms_faas_fixed);
pepke_m1_fixed_mses_eq  = shifman_errors(fpms_pepke_m1_fixed);

byrne_mses_eq           = shifman_errors(fpms_byrne);
byrne_fixed_mses_eq     = shifman_errors(fpms_byrne_fixed);

println("Blackwell Shifman RMSE: ", mean(bwell_mses_eq),            " STD: ", std(bwell_mses_eq), "\n",
"Blackwell Fixed Shifman RMSE: ",   mean(bwell_fixed_mses_eq),      " STD: ", std(bwell_fixed_mses_eq), "\n",
"Bhalla Shifman RMSE: ",            mean(bhalla_mses_eq),           " STD: ", std(bhalla_mses_eq), "\n",
"Bhalla Fixed Shifman RMSE: ",      mean(bhalla_fixed_mses_eq),     " STD: ", std(bhalla_fixed_mses_eq), "\n",
"Shifman Shifman RMSE: ",           mean(shifman_mses_eq),          " STD: ", std(shifman_mses_eq), "\n",
"Shifman Fixed Shifman RMSE: ",     mean(shifman_fixed_mses_eq),    " STD: ", std(shifman_fixed_mses_eq), "\n",
"Pepke M2 Shifman RMSE: ",          mean(pepke_m2_mses_eq),         " STD: ", std(pepke_m2_mses_eq), "\n",
"Pepke M2 Fixed Shifman RMSE: ",    mean(pepke_m2_fixed_mses_eq),   " STD: ", std(pepke_m2_fixed_mses_eq), "\n",
"Faas Shifman RMSE: ",              mean(faas_mses_eq),             " STD: ", std(faas_mses_eq), "\n",
"Faas Fixed Shifman RMSE: ",        mean(faas_fixed_mses_eq),       " STD: ", std(faas_fixed_mses_eq), "\n",
"Pepke M1 Fixed Shifman RMSE: ",    mean(pepke_m1_fixed_mses_eq),   " STD: ", std(pepke_m1_fixed_mses_eq), "\n",
"Byrne Shifman RMSE: ",             mean(byrne_mses_eq),            " STD: ", std(byrne_mses_eq), "\n",
"Byrne Fixed Shifman RMSE: ",       mean(byrne_fixed_mses_eq),      " STD: ", std(byrne_fixed_mses_eq)
)



hist_size_cm     = (17, 7.5)
hist_size_px     =  hist_size_cm .* dpc 
hist_size_units  =  hist_size_cm .* pt_per_cm
hist_px_per_unit = (hist_size_cm .* dpc ./ (hist_size_cm .* pt_per_cm))[1]


hist = CairoMakie.Figure(size=hist_size_units, fontsize=8)
ax11 = CairoMakie.Axis(hist[1,1], ylabel="Test RMSE", 
    xticks=(1.0:1.0:13.0, ["Scheme 1 + our rates";
                  "Scheme 1 + Kim et al.";
                  "Scheme 2 + our rates";
                  "Scheme 2 + Bhalla and Iyengar";
                  "Scheme 3 + our rates";
                  "Scheme 3 + Shifman et al. + our rates";
                  "Scheme 4 + our rates";
                  "Scheme 4 + Pepke et al.";
                  "Scheme 5 + our rates";
                  "Scheme 5 + Faas et al.";
                  "Scheme 5 + Pepke et al.";
                  "Scheme 6 + our rates";
                  "Scheme 6 + Byrne et al.";
    ]),
    xlabel="Scheme + parameters",
    xticklabelrotation=pi/9,
    xticklabelpad=5,
    xlabelpadding=15,
    xticklabelsize=8,
    yticklabelsize=8,
    xgridvisible=false,
)
h_111 = CairoMakie.violin!(ax11, 
    [
        1 * 1.0 *  ones(length(bwell_mses_test));
        2 * 1.0 *  ones(length(bwell_fixed_mses_test));
        3 * 1.0 *  ones(length(bhalla_mses_test));
        4 * 1.0 *  ones(length(bhalla_fixed_mses_test));
        5 * 1.0 *  ones(length(shifman_mses_test));
        6 * 1.0 *  ones(length(shifman_fixed_mses_test));
        7 * 1.0 *  ones(length(pepke_m2_mses_test));
        8 * 1.0 *  ones(length(pepke_m2_fixed_mses_test));
        9 * 1.0 *  ones(length(faas_mses_test));
       10 * 1.0 *  ones(length(faas_fixed_mses_test));
       11 * 1.0 *  ones(length(pepke_m1_fixed_mses_test));
       12 * 1.0 *  ones(length(byrne_mses_test));
       13 * 1.0 *  ones(length(byrne_fixed_mses_test))
    ], 
    [
        bwell_mses_test;
        bwell_fixed_mses_test;
        bhalla_mses_test;
        bhalla_fixed_mses_test;
        shifman_mses_test;
        shifman_fixed_mses_test;
        pepke_m2_mses_test;
        pepke_m2_fixed_mses_test;
        faas_mses_test;
        faas_fixed_mses_test;
        pepke_m1_fixed_mses_test;
        byrne_mses_test;
        byrne_fixed_mses_test
    ],
    color = [
        repeat([our_col], length(fpms_bwell));
        repeat([kim_col], length(fpms_bwell_fixed));
        repeat([our_col], length(fpms_bhalla));
        repeat([bhalla_col], length(fpms_bhalla_fixed));
        repeat([our_col], length(fpms_shifman));
        repeat([shifman_col], length(fpms_shifman_fixed));
        repeat([our_col], length(fpms_pepke_m2));
        repeat([pepke_col], length(fpms_pepke_m2_fixed));
        repeat([our_col], length(fpms_faas));
        repeat([faas_col], length(fpms_faas_fixed));
        repeat([pepke_col], length(fpms_pepke_m1_fixed));
        repeat([our_col], length(fpms_byrne));
        repeat([byrne_col], length(fpms_byrne_fixed));
    ]
)

###CairoMakie.save("Plots/violin.svg", hist, pt_per_unit=1, 
###    px_per_unit=hist_px_per_unit, size=hist_size_units)
CairoMakie.save("Plots/violin.png", hist, pt_per_unit=1, 
    px_per_unit=hist_px_per_unit, size=hist_size_units)



###h_111 = CairoMakie.hist!(ax11, bwell_mses,       bins=15, color=(our_col, 0.5), strokewidth=1, strokecolor=:black, label="Our rates")
###h_112 = CairoMakie.hist!(ax11, bwell_fixed_mses, bins=5, color=(kim_col, 0.5), strokewidth=1, strokecolor=:black, label="Kim et al. rates")
###CairoMakie.Legend(hist[1, 1], [h_111, h_112], ["Our rates", "Kim et al. rates"], "Scheme 1", 
###    titleposition=:left, framevisible=false, labelsize=8, nbanks=1)
###
###
###ax12 = CairoMakie.Axis(hist[2,2])
###h_121 = CairoMakie.hist!(ax12, pepke_m2_mses,       bins=15, color=(our_col, 0.5), strokewidth=1, strokecolor=:black, label="Our rates")
###h_122 = CairoMakie.hist!(ax12, pepke_m2_fixed_mses, bins=6, color=(pepke_col, 0.5), strokewidth=1, strokecolor=:black, label="Pepke et al. rates")
###CairoMakie.Legend(hist[1, 2], [h_121, h_122], ["Our rates", "Pepke et al. rates"], "Scheme 2", 
###    titleposition=:left, framevisible=false, labelsize=8, nbanks=1)
###
###
###ax21 = CairoMakie.Axis(hist[4,1])
###h_211 = CairoMakie.hist!(ax21, faas_mses,           bins=7, color=(our_col, 0.5),   strokewidth=1, strokecolor=:black, label="Our rates")
###h_212 = CairoMakie.hist!(ax21, faas_fixed_mses,     bins=3, color=(faas_col, 0.5),  strokewidth=1, strokecolor=:black, label="Faas et al. rates")
###h_213 = CairoMakie.hist!(ax21, pepke_m1_fixed_mses, bins=4, color=(pepke_col, 0.5), strokewidth=1, strokecolor=:black, label="Pepke et al. rates")
###CairoMakie.Legend(hist[3, 1], [h_211, h_212, h_213], ["Our rates", "Faas et al. rates", "Pepke et al. rates"], "Scheme 3", 
###    titleposition=:left, framevisible=false, labelsize=8, nbanks=2)
###
###
###ax22 = CairoMakie.Axis(hist[4,2])
###h_221 = CairoMakie.hist!(ax22, byrne_mses,       bins=15, color=(our_col, 0.5),    strokewidth=1, strokecolor=:black, label="Our rates")
###h_222 = CairoMakie.hist!(ax22, byrne_fixed_mses, bins=5, color=(byrne_col, 0.5),  strokewidth=1, strokecolor=:black, label="Byrne et al. rates")
###CairoMakie.Legend(hist[3, 2], [h_221, h_222], ["Our rates", "Byrne et al. rates"], "Scheme 4", 
###    titleposition=:left, framevisible=false, labelsize=8, nbanks=1)
###
###CairoMakie.Label(hist[:,0], "Count", rotation=pi/2)
###CairoMakie.Label(hist[5,:], "RMSE")
###
###CairoMakie.rowsize!(hist.layout, 1, Relative(0.125))
###CairoMakie.rowsize!(hist.layout, 3, Relative(0.125))
###CairoMakie.rowsize!(hist.layout, 5, Relative(0.00))
###
###CairoMakie.colsize!(hist.layout, 1, Relative(0.5))
###CairoMakie.colsize!(hist.layout, 2, Relative(0.5))
###
###CairoMakie.rowgap!(hist.layout, 4, 5)
###CairoMakie.colgap!(hist.layout, 1, 10)
###
###
###CairoMakie.save("Plots/hist.svg", hist, pt_per_unit=1, 
###    px_per_unit=hist_px_per_unit, size=hist_size_units)
###########################
### MSE HISTS PLOTS END ###
###########################
#######################
### VAL PLOTS START ###
#######################

val_pops_bigplot = idx_to_val_vec(full_pop, val_idxs_bwell[1]);

bp_size_cm     = (17, 22.23)
bp_size_px     =  bp_size_cm .* dpc 
bp_size_units  =  bp_size_cm .* pt_per_cm
bp_px_per_unit = (bp_size_cm .* dpc ./ (bp_size_cm .* pt_per_cm))[1]
bp_colors = [our_col; kim_col; our_col; bhalla_col; our_col; shifman_col; our_col; pepke_col; our_col; faas_col; pepke_col; our_col; byrne_col];
###bp_rowlist = [2;3;5;6;8;9;11;12;14;15;16;18;19]
bp_rowlist = collect(1:13)
bp = bigplot([fpms_bwell[1];    fpms_bwell_fixed[1];
         fpms_bhalla[1];        fpms_bhalla_fixed[1];
         fpms_shifman[1];       fpms_shifman_fixed[1];
         fpms_pepke_m2[1];      fpms_pepke_m2_fixed[1];
         fpms_faas[1];          fpms_faas_fixed[1]; fpms_pepke_m1_fixed[1];
         fpms_byrne[1];         fpms_byrne_fixed[1]],
         val_pops_bigplot,
         bp_colors,
         bp_rowlist,
         names = ["Our rates"; "Kim et al.";
                  "Our rates"; "Hayer and Bhalla";
                  "Our rates"; "Shifman et al.";
                  "Our rates"; "Peple et al.";
                  "Our rates"; "Faas et al."; "Pepke et al.";
                  "Our rates"; "Byrne et al."],
        figsize=bp_size_units
)

CairoMakie.Legend(bp[1:2, 8], [bp.content[1].scene.plots[6],  bp.content[8].scene.plots[end]], 
    ["Our rates", "Kim et al. \n rates"], "Scheme 1", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
    colors=[our_col, kim_col], labeljustification=:center)
CairoMakie.Legend(bp[3:4, 8], [bp.content[15].scene.plots[end], bp.content[22].scene.plots[end]], 
    ["Our rates", "Hayer and \n Bhalla rates"], "Scheme 2", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
    colors=[our_col, bhalla_col], labeljustification=:center)
CairoMakie.Legend(bp[5:6, 8], [bp.content[29].scene.plots[end], bp.content[36].scene.plots[end]], 
    ["Our rates", "Shifman et al. \n rates"], "Scheme 3", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
    colors=[our_col, shifman_col], labeljustification=:center)
CairoMakie.Legend(bp[7:8, 8], [bp.content[43].scene.plots[end], bp.content[50].scene.plots[end]], 
    ["Our rates", "Pepke al. \n rates"], "Scheme 4", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
    colors=[our_col, pepke_col], labeljustification=:center)
CairoMakie.Legend(bp[9:11, 8], [bp.content[57].scene.plots[end], bp.content[64].scene.plots[end], bp.content[71].scene.plots[end]], 
    ["Our rates", "Faas et al. \n rates", "Pepke al. \n rates"], "Scheme 5", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
    colors=[our_col, faas_col, pepke_col], labeljustification=:center)
CairoMakie.Legend(bp[12:13, 8], [bp.content[78].scene.plots[end], bp.content[85].scene.plots[end]], 
    ["Our rates", "Byrne et al."], "Scheme 6", titleposition=:top, framevisible=false, labelsize=8, nbanks=1, 
    colors=[our_col, byrne_col], labeljustification=:center)

bp_ly = CairoMakie.Label(bp[:,0], "ΔF/F₀", rotation=pi/2, fontsize=8)
bp_lx = CairoMakie.Label(bp[bp_rowlist[end]+1,:], "Time (ms)", fontsize=8)


for i in 1:7
    colsize!(bp.layout, i, Relative(1/8))
end
colsize!(bp.layout, 8, Relative(1/6))

for i in bp_rowlist
    CairoMakie.rowsize!(bp.layout, i, Relative(1/13))
end

CairoMakie.colgap!(bp.layout, 5)
CairoMakie.rowgap!(bp.layout, 2)

CairoMakie.rowgap!(bp.layout, 14, 5)

###CairoMakie.save("Plots/bigplot.svg", bp, pt_per_unit=1, 
###    px_per_unit=bp_px_per_unit, size=bp_size_units)
CairoMakie.save("Plots/bigplot.png", bp, pt_per_unit=1, 
    px_per_unit=bp_px_per_unit, size=bp_size_units)

#####################
### VAL PLOTS END ###
#####################
###############################
### PRACTICAL EQ PLOT START ###
###############################

### loading in Shifman et al. 2006 data
f_shifman = XLSX.readxlsx("Shifman_2006.xlsx");
X_shifman = f_shifman["Shifman_2006.csv"]["A"][2:end] .* 1e-6;
Y_shifman = [i for i in f_shifman["Shifman_2006.csv"]["B"][2:end]];

ca_range    = 10 .^ LinRange(-9, log10(7e-5), 100);
pract_pop   = gen_dummy_data(ca_range);

eq_size_cm     = (17, 22.23)
eq_size_px     =  eq_size_cm .* dpc 
eq_size_units  =  eq_size_cm .* pt_per_cm
eq_px_per_unit = (eq_size_cm .* dpc ./ (eq_size_cm .* pt_per_cm))[1]

lims = (nothing, nothing, nothing, nothing)

bwell_CaM_tup = get_multiple_eqs(fpms_bwell, pract_pop, get_bwell_eqs);
bwell_fixed_CaM_tup = get_multiple_eqs(fpms_bwell_fixed, pract_pop, get_bwell_eqs);

f_eqs, f_eqs_ax11 = plot_CaM_eqs(ca_range, bwell_CaM_tup;       i=1, j=1, title="Scheme 1", 
    color=our_col,  f=nothing, xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates", figsize=eq_size_units);
_, _      = plot_CaM_eqs(ca_range, bwell_fixed_CaM_tup; i=1, j=1, title="Scheme 1", 
    color=kim_col, f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Kim et al. rates");
CairoMakie.scatter!(f_eqs_ax11, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
CairoMakie.Legend(f_eqs[1, 2], f_eqs_ax11, "Scheme 1", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax11, -0.3, 4.3) 



bhalla_CaM_tup        = get_multiple_eqs(fpms_bhalla, pract_pop, get_bhalla_eqs);
bhalla_fixed_CaM_tup  = get_multiple_eqs(fpms_bhalla_fixed, pract_pop, get_bhalla_eqs);

f_eqs, f_eqs_ax12 = plot_CaM_eqs(ca_range, bhalla_CaM_tup;       i=2, j=1, title="Scheme 2", 
    color=our_col,  f=f_eqs, new_axis=true,   xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, bhalla_fixed_CaM_tup;         i=2, j=1, title="Scheme 2", 
    color=bhalla_col, f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Hayer and Bhalla rates");
CairoMakie.scatter!(f_eqs_ax12, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
CairoMakie.Legend(f_eqs[2, 2], f_eqs_ax12, "Scheme 2", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax12, -0.3, 4.3)



shifman_CaM_tup        = get_multiple_eqs(fpms_shifman, pract_pop, get_shifman_eqs);
shifman_fixed_CaM_tup  = get_multiple_eqs(fpms_shifman_fixed, pract_pop, get_shifman_eqs);

f_eqs, f_eqs_ax13 = plot_CaM_eqs(ca_range, shifman_CaM_tup;       i=3, j=1, title="Scheme 3", 
    color=our_col,  f=f_eqs, new_axis=true,   xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, shifman_fixed_CaM_tup;         i=3, j=1, title="Scheme 3", 
    color=shifman_col, f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Shifman et al. \n + our rates");
CairoMakie.scatter!(f_eqs_ax13, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
CairoMakie.Legend(f_eqs[3, 2], f_eqs_ax13, "Scheme 3", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax12, -0.3, 4.3)



pepke_m2_CaM_tup        = get_multiple_eqs(fpms_pepke_m2, pract_pop, get_pepke_m2_eqs);
pepke_m2_fixed_CaM_tup  = get_multiple_eqs(fpms_pepke_m2_fixed, pract_pop, get_pepke_m2_eqs);

f_eqs, f_eqs_ax14 = plot_CaM_eqs(ca_range, pepke_m2_CaM_tup;       i=4, j=1, title="Scheme 4", 
    color=our_col,  f=f_eqs, new_axis=true,   xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, pepke_m2_fixed_CaM_tup;         i=4, j=1, title="Scheme 4", 
    color=pepke_col, f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Pepke et al. rates");
CairoMakie.scatter!(f_eqs_ax14, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
CairoMakie.Legend(f_eqs[4, 2], f_eqs_ax14, "Scheme 4", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax14, -0.3, 4.3) 




faas_CaM_tup        = get_multiple_eqs(fpms_faas, pract_pop, get_faas_eqs);
faas_fixed_CaM_tup  = get_multiple_eqs(fpms_faas_fixed, pract_pop, get_faas_eqs);
pepke_m1_fixed_CaM_tup  = get_multiple_eqs(fpms_pepke_m1_fixed, pract_pop, get_faas_eqs);

f_eqs, f_eqs_ax15 = plot_CaM_eqs(ca_range, faas_CaM_tup;       i=5, j=1, title="Scheme 5", 
    color=our_col,  f=f_eqs, new_axis=true, xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, faas_fixed_CaM_tup;         i=5, j=1, title="Scheme 5", 
    color=faas_col, f=f_eqs,                xlabel=nothing, ylabel=nothing, limits=lims, label="Faas et al. rates");
_, _      = plot_CaM_eqs(ca_range, pepke_m1_fixed_CaM_tup;         i=2, j=1, title="Scheme 5", 
    color=pepke_col, f=f_eqs,               xlabel=nothing, ylabel=nothing, limits=lims, label="Pepke et al. rates");
CairoMakie.scatter!(f_eqs_ax15, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
CairoMakie.Legend(f_eqs[5, 2], f_eqs_ax15, "Scheme 5", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax15, -0.3, 4.3) 



byrne_CaM_tup        = get_multiple_eqs(fpms_byrne, pract_pop, get_byrne_eqs);
byrne_fixed_CaM_tup  = get_multiple_eqs(fpms_byrne_fixed, pract_pop, get_byrne_eqs);

f_eqs, f_eqs_ax16 = plot_CaM_eqs(ca_range, byrne_CaM_tup;       i=6, j=1, title="Scheme 6", 
    color=our_col,  f=f_eqs, new_axis=true, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, faas_fixed_CaM_tup;          i=6, j=1, title="Scheme 6", 
    color=byrne_col, f=f_eqs,                ylabel=nothing, limits=lims, label="Byrne et al. rates");
CairoMakie.scatter!(f_eqs_ax16, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.3));
CairoMakie.Legend(f_eqs[6, 2], f_eqs_ax16, "Scheme 6", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax16, -0.3, 4.3) 


###CairoMakie.Label(f_eqs[:,0], L"Bound $\textrm{Ca}^{2+}$ per CaM", rotation=pi/2)
CairoMakie.Label(f_eqs[:,0], "# Bound Ca²⁺ per CaM", rotation=pi/2)
CairoMakie.colgap!(f_eqs.layout, 1, 10)
##CairoMakie.save("Plots/Shifman_joint.svg", f_eqs, pt_per_unit=1, px_per_unit = eq_px_per_unit, size=eq_size_units)
CairoMakie.save("Plots/Shifman_joint.png", f_eqs, pt_per_unit=1, px_per_unit = eq_px_per_unit, size=eq_size_units)

#############################
### PRACTICAL EQ PLOT END ###
#############################
################################
### PRACTICAL DYN PLOT START ###
################################
###dyn_pop   = gen_dyn_data()
######dyn_pop   = gen_dyn_data_2Hz()
###t = dyn_pop[1].time;
###
###bwell_CaM_dyn_tups = get_multiple_dyns(fpms_bwell, dyn_pop, get_bwell_dyns);
###bwell_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_bwell_fixed, dyn_pop, get_bwell_dyns);
###
###bhalla_CaM_dyn_tups = get_multiple_dyns(fpms_bhalla, dyn_pop, get_bhalla_dyns);
###bhalla_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_bhalla_fixed, dyn_pop, get_bhalla_dyns);
###
###shifman_CaM_dyn_tups = get_multiple_dyns(fpms_shifman, dyn_pop, get_shifman_dyns);
###shifman_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_shifman_fixed, dyn_pop, get_shifman_dyns);
###
###pepke_m2_CaM_dyn_tups = get_multiple_dyns(fpms_pepke_m2, dyn_pop, get_bwell_CN_dyns);
###pepke_m2_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_pepke_m2_fixed, dyn_pop, get_bwell_CN_dyns);
###
###faas_CaM_dyn_tups = get_multiple_dyns(fpms_faas, dyn_pop, get_faas_dyns);
###faas_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_faas_fixed, dyn_pop, get_faas_dyns);
###pepke_m1_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_pepke_m1_fixed, dyn_pop, get_pepke_m1_dyns);
###
###byrne_CaM_dyn_tups = get_multiple_dyns(fpms_byrne, dyn_pop, get_byrne_dyns);
###byrne_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_byrne_fixed, dyn_pop, get_byrne_dyns);
###
###
###dyn_size_cm     = (17, 10)
###dyn_size_px     =  dyn_size_cm .* dpc 
###dyn_size_units  =  dyn_size_cm .* pt_per_cm
###dyn_px_per_unit = (dyn_size_cm .* dpc ./ (dyn_size_cm .* pt_per_cm))[1]
######lims = (nothing, nothing, nothing, nothing)
###lims = [(-100, 1550, -1, 21);
###        (-100, 1550, -1, 21);
###        (-100, 1550, -1, 3)]
###
###lw = 1
###
###f_dyn = plot_CaM_dyns(t, bwell_CaM_dyn_tups;
###                j=1, title="Our rates", color=(our_col, 0.3),
###                f = nothing,
###                figsize=dyn_size_units, lw=lw, limits=lims);
###f_dyn = plot_CaM_dyns(t, bwell_fixed_CaM_dyn_tups;
###                j=2, title="Kim et al.", color=(kim_col, 0.3),
###                f = f_dyn,
###                lw=1, limits=lims,
###                ylabs=false, yticks=false);
###elem_our = LineElement(color = our_col, linestyle = nothing);
###elem_kim = LineElement(color = kim_col, linestyle = nothing);
###CairoMakie.Legend(f_dyn[1, 1:2], [elem_our, elem_kim], 
###    ["Our \n rates", "Kim \n et al."], "Scheme 1", framevisible=false, labelsize=8, nbanks=1, 
###    colors=[(our_col, 1.0), byrne_col], labeljustification=:center)
###
###
###f_dyn = plot_CaM_dyns(t, bhalla_CaM_dyn_tups;
###                j=3, title="Our rates", color=(our_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###f_dyn = plot_CaM_dyns(t, bhalla_fixed_CaM_dyn_tups;
###                j=4, title="Pepke et al.", color=(bhalla_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###elem_bhalla = LineElement(color = bhalla_col, linestyle = nothing);
###CairoMakie.Legend(f_dyn[1, 3:4], [elem_our, elem_bhalla],# [f_dyn.content[8].scene.plots[end], f_dyn.content[11].scene.plots[end]], 
###    ["Our \n rates", "Hayer \n and Bhalla"], "Scheme 2", framevisible=false, labelsize=8, nbanks=1, 
###    colors=[our_col, byrne_col], labeljustification=:center)
###
###
###f_dyn = plot_CaM_dyns(t, shifman_CaM_dyn_tups;
###                j=5, title="Our rates", color=(our_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###f_dyn = plot_CaM_dyns(t, shifman_fixed_CaM_dyn_tups;
###                j=6, title="Pepke et al.", color=(shifman_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###elem_shifman = LineElement(color = shifman_col, linestyle = nothing);
###CairoMakie.Legend(f_dyn[1, 5:6], [elem_our, elem_shifman], # [f_dyn.content[15].scene.plots[end], f_dyn.content[18].scene.plots[end]], 
###    ["Our \n rates", "Shifman \n et al."], "Scheme 3", framevisible=false, labelsize=8, nbanks=1, 
###    colors=[our_col, byrne_col], labeljustification=:center)
###
###
###f_dyn = plot_CaM_dyns(t, pepke_m2_CaM_dyn_tups;
###                j=7, title="Our rates", color=(our_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###f_dyn = plot_CaM_dyns(t, pepke_m2_fixed_CaM_dyn_tups;
###                j=8, title="Pepke et al.", color=(pepke_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###elem_pepke = LineElement(color = pepke_col, linestyle = nothing);
###CairoMakie.Legend(f_dyn[1, 7:8], [elem_our, elem_pepke], # [f_dyn.content[22].scene.plots[end], f_dyn.content[25].scene.plots[end]], 
###    ["Our \n rates", "Pepke \n et al."], "Scheme 4", framevisible=false, labelsize=8, nbanks=1, 
###    colors=[our_col, byrne_col], labeljustification=:center)
###
###
###f_dyn = plot_CaM_dyns(t, faas_CaM_dyn_tups;
###                j=9, title="Our rates", color=(our_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###f_dyn = plot_CaM_dyns(t, faas_fixed_CaM_dyn_tups;
###                j=10, title="Faas et al.", color=(faas_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###f_dyn = plot_CaM_dyns(t, pepke_m1_fixed_CaM_dyn_tups;
###                j=11, title="Pepke et al.", color=(pepke_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###elem_faas = LineElement(color = faas_col, linestyle = nothing);
###CairoMakie.Legend(f_dyn[1, 9:11], [elem_our, elem_faas, elem_pepke], #[f_dyn.content[29].scene.plots[end], f_dyn.content[32].scene.plots[end], f_dyn.content[35].scene.plots[end]], 
###    ["Our \n rates", "Faas \n et al.", "Pepke \n et al."], "Scheme 5", framevisible=false, labelsize=8, nbanks=2, 
###    colors=[our_col, byrne_col], labeljustification=:center)
###
###
###f_dyn = plot_CaM_dyns(t, byrne_CaM_dyn_tups;
###                j=12, title="Our rates", color=(our_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###f_dyn = plot_CaM_dyns(t, byrne_fixed_CaM_dyn_tups;
###                j=13, title="Byrne et al.", color=(byrne_col, 0.3),
###                f = f_dyn,
###                lw=lw, limits=lims,
###                ylabs=false, yticks=false);
###elem_byrne = LineElement(color = byrne_col, linestyle = nothing);
###CairoMakie.Legend(f_dyn[1, 12:13], [elem_our, elem_byrne], # [f_dyn.content[39].scene.plots[end], f_dyn.content[42].scene.plots[end]], 
###    ["Our \n rates", "Byrne \n et al."], "Scheme 6", framevisible=false, labelsize=8, nbanks=1, 
###    colors=[our_col, byrne_col], labeljustification=:center)
###
######CairoMakie.Label(f_dyn[1,0], L"$\textrm{CaM}_0$", rotation=pi/2)
######CairoMakie.Label(f_dyn[2,0], L"$\textrm{CaM}_p$", rotation=pi/2)
######CairoMakie.Label(f_dyn[3,0], L"$\textrm{CaM}_f$", rotation=pi/2)
###
###CairoMakie.Label(f_dyn[2,0], "CaM0Ca (μM)", rotation=pi/2)
###CairoMakie.Label(f_dyn[3,0], "CaM1-3Ca (μM)", rotation=pi/2)
###CairoMakie.Label(f_dyn[4,0], "CaM4Ca (μM)", rotation=pi/2)
###
###CairoMakie.Label(f_dyn[5,:], "Time (ms)")
###
###for i in 1:4
###    rowsize!(f_dyn.layout, i, Relative(1/3.8))
###end
###for i in 1:13
###    colsize!(f_dyn.layout, i, Relative(1/13))
###end
###
###CairoMakie.colgap!(f_dyn.layout, 5)
###
###colgap!(f_dyn.layout, 5)
###rowgap!(f_dyn.layout, 10)
###
######CairoMakie.save("Plots/dyn_joint.svg", f_dyn, pt_per_unit=1, px_per_unit = dyn_px_per_unit, size=dyn_size_units)
###CairoMakie.save("Plots/dyn_joint.png", f_dyn, pt_per_unit=1, px_per_unit = dyn_px_per_unit, size=dyn_size_units)
###
###
##############################
### PRACTICAL DYN PLOT END ###
##############################

##############################
### INTEGRATION PLOT START ###
##############################
###f_range = 2:10:102;
###ca_amp  = 12e-6;
f_range = 2:20:502
ca_amp = 0.7e-6

bwell_lines       = calc_integration(fpms_bwell, f_range, ca_amp, get_bwell_dyns);
bwell_fixed_lines = calc_integration(fpms_bwell_fixed, f_range, ca_amp, get_bwell_dyns);

bhalla_lines       = calc_integration(fpms_bhalla, f_range, ca_amp, get_bhalla_dyns);
bhalla_fixed_lines = calc_integration(fpms_bhalla_fixed, f_range, ca_amp, get_bhalla_dyns);

shifman_lines       = calc_integration(fpms_shifman, f_range, ca_amp, get_shifman_dyns);
shifman_fixed_lines = calc_integration(fpms_shifman_fixed, f_range, ca_amp, get_shifman_dyns);

pepke_m2_lines       = calc_integration(fpms_pepke_m2, f_range, ca_amp, get_bwell_CN_dyns);
pepke_m2_fixed_lines = calc_integration(fpms_pepke_m2_fixed, f_range, ca_amp, get_bwell_CN_dyns);

faas_lines       = calc_integration(fpms_faas, f_range, ca_amp, get_faas_dyns);
faas_fixed_lines = calc_integration(fpms_faas_fixed, f_range, ca_amp, get_faas_dyns);
pepke_m1_fixed_lines = calc_integration(fpms_pepke_m1_fixed, f_range, ca_amp, get_pepke_m1_dyns);

byrne_lines       = calc_integration(fpms_byrne, f_range, ca_amp, get_byrne_dyns);
byrne_fixed_lines = calc_integration(fpms_byrne_fixed, f_range, ca_amp, get_byrne_dyns);


int_size_cm     = (17, 10)
int_size_px     =  int_size_cm .* dpc 
int_size_units  =  int_size_cm .* pt_per_cm
int_px_per_unit = (int_size_cm .* dpc ./ (int_size_cm .* pt_per_cm))[1]

lw = 1
line_alpha = 0.6
ribbon_alpha = 0.2
###int_limits = [(nothing, nothing, -1e-3, 0.016); (nothing, nothing, -1e-4, 0.0015)]
###int_limits = [(nothing, nothing, -1e-4, 0.0015); (nothing, nothing, -1e-4, 0.0015)]
int_limits = [(nothing, nothing, nothing, nothing); (nothing, nothing, nothing, nothing)]

f_int = CairoMakie.Figure(size=int_size_units, fontsize=8)

int_ax_11 = plot_int_lines(f_int, collect(f_range), bwell_lines, 1, 2, 
    our_col, line_alpha, ribbon_alpha, lw, (false, true), (false, true), int_limits);
plot_int_lines(int_ax_11, collect(f_range), bwell_fixed_lines,   1, 2, kim_col, line_alpha, ribbon_alpha, lw)
int_ax_21 = plot_int_lines(f_int, collect(f_range), bwell_lines, 2, 2, 
    our_col, line_alpha, ribbon_alpha, lw, (true, true), (true, true), int_limits);
plot_int_lines(int_ax_21, collect(f_range), bwell_fixed_lines,   2, 2, kim_col, line_alpha, ribbon_alpha, lw);
elem_our = LineElement(color = our_col, linestyle = nothing);
elem_kim = LineElement(color = kim_col, linestyle = nothing);
CairoMakie.Legend(f_int[0, 2], [elem_our, elem_kim], 
    ["Our \n rates", "Kim \n et al."], "Scheme 1", framevisible=false, labelsize=8, nbanks=1, 
    labeljustification=:center, valign=:top)


int_ax_12 = plot_int_lines(f_int, collect(f_range), bhalla_lines, 1, 3, 
    our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
plot_int_lines(int_ax_12, collect(f_range), bhalla_fixed_lines,   1, 3, bhalla_col, line_alpha, ribbon_alpha, lw);
int_ax_22 = plot_int_lines(f_int, collect(f_range), bhalla_lines, 2, 3, 
    our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
plot_int_lines(int_ax_22, collect(f_range), bhalla_fixed_lines,   2, 3, bhalla_col, line_alpha, ribbon_alpha, lw);
elem_bhalla = LineElement(color = bhalla_col, linestyle = nothing);
CairoMakie.Legend(f_int[0, 3], [elem_our, elem_bhalla], 
    ["Our \n rates", "Hayer \n and Bhalla"], "Scheme 2", framevisible=false, labelsize=8, nbanks=1, 
    labeljustification=:center, valign=:top)


int_ax_13 = plot_int_lines(f_int, collect(f_range), shifman_lines, 1, 4, 
    our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
plot_int_lines(int_ax_13, collect(f_range), shifman_fixed_lines,   1, 4, shifman_col, line_alpha, ribbon_alpha, lw);
int_ax_23 = plot_int_lines(f_int, collect(f_range), shifman_lines, 2, 4, 
    our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
plot_int_lines(int_ax_23, collect(f_range), shifman_fixed_lines,   2, 4, shifman_col, line_alpha, ribbon_alpha, lw);
elem_shifman = LineElement(color = shifman_col, linestyle = nothing);
CairoMakie.Legend(f_int[0, 4], [elem_our, elem_shifman], 
    ["Our \n rates", "Shifman \n et al."], "Scheme 3", framevisible=false, labelsize=8, nbanks=1, 
    labeljustification=:center, valign=:top)



int_ax_14 = plot_int_lines(f_int, collect(f_range), pepke_m2_lines, 1, 5, 
    our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
plot_int_lines(int_ax_14, collect(f_range), pepke_m2_fixed_lines,   1, 5, pepke_col, line_alpha, ribbon_alpha, lw);
int_ax_24 = plot_int_lines(f_int, collect(f_range), pepke_m2_lines, 2, 5, 
    our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
plot_int_lines(int_ax_24, collect(f_range), pepke_m2_fixed_lines,   2, 5, pepke_col, line_alpha, ribbon_alpha, lw);
elem_pepke = LineElement(color = pepke_col, linestyle = nothing);
CairoMakie.Legend(f_int[0, 5], [elem_our, elem_pepke], 
    ["Our \n rates", "Pepke \n et al."], "Scheme 4", framevisible=false, labelsize=8, nbanks=1, 
    labeljustification=:center, valign=:top)


int_ax_15 = plot_int_lines(f_int, collect(f_range), faas_lines,     1, 6, 
    our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
plot_int_lines(int_ax_15, collect(f_range), faas_fixed_lines,       1, 6, faas_col, line_alpha, ribbon_alpha, lw);
plot_int_lines(int_ax_15, collect(f_range), pepke_m1_fixed_lines,   1, 6, pepke_col, line_alpha, ribbon_alpha, lw);
int_ax_25 = plot_int_lines(f_int, collect(f_range), faas_lines,     2, 6, 
    our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
plot_int_lines(int_ax_25, collect(f_range), faas_fixed_lines,       2, 6, faas_col, line_alpha, ribbon_alpha, lw);
plot_int_lines(int_ax_25, collect(f_range), pepke_m1_fixed_lines,   2, 6, pepke_col, line_alpha, ribbon_alpha, lw);
elem_faas = LineElement(color = faas_col, linestyle = nothing);
CairoMakie.Legend(f_int[0, 6], [elem_our, elem_faas, elem_pepke], 
    ["Our rates", "Faas et al.", "Pepke et al."], "Scheme 5", framevisible=false, labelsize=8, nbanks=1, 
    labeljustification=:center, valign=:top)
    

int_ax_16 = plot_int_lines(f_int, collect(f_range), byrne_lines,     1, 7, 
    our_col, line_alpha, ribbon_alpha, lw, (false, false), (false, false), int_limits);
plot_int_lines(int_ax_16, collect(f_range), byrne_fixed_lines,       1, 7, byrne_col, line_alpha, ribbon_alpha, lw);
int_ax_26 = plot_int_lines(f_int, collect(f_range), byrne_lines,     2, 7, 
    our_col, line_alpha, ribbon_alpha, lw, (true, false), (true, false), int_limits);
plot_int_lines(int_ax_26, collect(f_range), byrne_fixed_lines,       2, 7, byrne_col, line_alpha, ribbon_alpha, lw);
elem_byrne = LineElement(color = byrne_col, linestyle = nothing);
CairoMakie.Legend(f_int[0, 7], [elem_our, elem_byrne], 
    ["Our \n rates", "Byrne \n et al."], "Scheme 6", framevisible=false, labelsize=8, nbanks=1, 
    labeljustification=:center, valign=:top)

CairoMakie.Label(f_int[1,1], "Partially bound CaM", rotation=pi/2)
CairoMakie.Label(f_int[2,1], "Fully bound CaM", rotation=pi/2)
CairoMakie.Label(f_int[1:2,0], "AUC", rotation=pi/2)
CairoMakie.Label(f_int[3,:], "Frequency (Hz)")

for i in 1:2
    rowsize!(f_int.layout, i, Relative(1/3))
end
rowsize!(f_int.layout, 0, Relative(1/2.5))
for i in 2:7
    colsize!(f_int.layout, i, Relative(1/6))
end

CairoMakie.colgap!(f_int.layout, 5)
CairoMakie.rowgap!(f_int.layout, 10)

###CairoMakie.save("Plots/dyn_joint.svg", f_dyn, pt_per_unit=1, px_per_unit = dyn_px_per_unit, size=dyn_size_units)
CairoMakie.save("Plots/integration2.png", f_int, pt_per_unit=1, px_per_unit = int_px_per_unit, size=int_size_units)


##############################
### PRACTICAL DYN PLOT END ###
##############################
#########################
### PCORR PLOTS START ###
#########################


coef_bwell = get_coef_lists(fpms_bwell, 
    (:tv_kon_1, :tv_kon_2,
     :tv_Kd_1,  :tv_Kd_2);
);
dd_bwell = Dict(:tv_kon_1 => "k₁",   :tv_kon_2 => "k₃",
                :tv_Kd_1  => "K_D₁", :tv_Kd_2  => "K_D₂")


coef_bhalla = get_coef_lists(fpms_bhalla, 
    (:tv_kon_1, :tv_kon_2, :tv_kon_3,
     :tv_Kd_1,  :tv_Kd_2,  :tv_Kd_3);
);
dd_bhalla = Dict(:tv_kon_1 => "k₁",   :tv_kon_2 => "k₃",   :tv_kon_3 => "k₅",
                 :tv_Kd_1  => "K_D₁", :tv_Kd_2  => "K_D₂", :tv_Kd_3  => "K_D₃")


coef_shifman = get_coef_lists(fpms_shifman, 
    (:tv_kon_1, :tv_kon_2, :tv_kon_3, :tv_kon_4,
     :tv_Kd_1,  :tv_Kd_2,  :tv_Kd_3,  :tv_Kd_4);
);
dd_shifman = Dict(:tv_kon_1 => "k₁",   :tv_kon_2 => "k₃",   :tv_kon_3 => "k₅",   :tv_kon_4 => "k₇",
                  :tv_Kd_1  => "K_D₁", :tv_Kd_2  => "K_D₂", :tv_Kd_3  => "K_D₃", :tv_Kd_4  => "K_D₄")


pepke_m2_coef_counternames = zip(
    (:tv_kon_TC, :tv_kon_RC, :tv_Kd_TC, :tv_Kd_RC),
    (:tv_kon_TN, :tv_kon_RN, :tv_Kd_TN, :tv_Kd_RN),
)
pepke_m2_keys_ordered = (:tv_kon_TC, :tv_kon_RC, :tv_kon_TN, :tv_kon_RN,
                         :tv_Kd_TC,  :tv_Kd_RC,  :tv_Kd_TN,  :tv_Kd_RN)
pepke_m2_labs_ordered = ("k₁", "k₃", "k₅", "k₇",
                         "K_D₁", "K_D₂", "K_D₃", "K_D₄")
coef_pepke_m2 = get_coef_lists(fpms_pepke_m2, 
    (:tv_kon_RC, :tv_kon_RN, :tv_kon_TC, :tv_kon_TN,
     :tv_Kd_RC,  :tv_Kd_RN,  :tv_Kd_TC,  :tv_Kd_TN);
     sort = true,
     kd_C_name=:tv_Kd_TC,
     kd_N_name=:tv_Kd_TN,
     counternames = pepke_m2_coef_counternames
);
coef_pepke_m2 = NamedTuple{pepke_m2_keys_ordered}(coef_pepke_m2)
dd_pepke_m2 = Dict(zip(pepke_m2_keys_ordered, pepke_m2_labs_ordered))


faas_coef_counternames = zip(
    (:tv_kon_TC, :tv_kon_RC, :tv_Kd_TC, :tv_Kd_RC),
    (:tv_kon_TN, :tv_kon_RN, :tv_Kd_TN, :tv_Kd_RN),
)
faas_keys_ordered = (:tv_kon_TC, :tv_kon_RC, :tv_kon_TN, :tv_kon_RN,
                     :tv_Kd_TC,  :tv_Kd_RC,  :tv_Kd_TN,  :tv_Kd_RN)
faas_labs_ordered = ("k₁", "k₃", "k₅", "k₇",
                     "K_D₁", "K_D₂", "K_D₃", "K_D₄")
coef_faas = get_coef_lists(fpms_faas, 
    (:tv_kon_RC, :tv_kon_RN, :tv_kon_TC, :tv_kon_TN,
     :tv_Kd_RC,  :tv_Kd_RN,  :tv_Kd_TC,  :tv_Kd_TN);
     sort = true,
     kd_C_name=:tv_Kd_TC,
     kd_N_name=:tv_Kd_TN,
     counternames = faas_coef_counternames
);
coef_faas = NamedTuple{faas_keys_ordered}(coef_faas)
dd_faas = Dict(zip(faas_keys_ordered, faas_labs_ordered))


byrne_coef_counternames = zip(
    (:tv_k01_N, :tv_k02_N, :tv_k13_N, :tv_k23_N, :tv_K01d_N, :tv_K13d_N, :tv_K02d_N),
    (:tv_k01_C, :tv_k02_C, :tv_k13_C, :tv_k23_C, :tv_K01d_C, :tv_K13d_C, :tv_K02d_C)
)
byrne_keys_ordered = (:tv_k01_N,  :tv_k02_N,  :tv_k13_N, :tv_k23_N,
                      :tv_k01_C,  :tv_k02_C,  :tv_k13_C, :tv_k23_C,
                      :tv_K01d_N, :tv_K13d_N, :tv_K02d_N, 
                      :tv_K01d_C, :tv_K13d_C, :tv_K02d_C
) 
byrne_labs_ordered = ("kⁿ₀₁", "kⁿ₀₂", "kⁿ₁₃", "kⁿ₂₃",
                      "kᶜ₀₁", "kᶜ₀₂", "kᶜ₁₃", "kᶜ₂₃",
                      "K_Dⁿ₀₁", "K_Dⁿ₁₃", "K_Dⁿ₀₂", 
                      "K_Dᶜ₀₁", "K_Dᶜ₁₃", "K_Dᶜ₀₂"
)
coef_byrne = get_coef_lists(fpms_byrne, 
    (:tv_k01_N,     :tv_k02_N,   :tv_k13_N, :tv_k23_N,
     :tv_K01d_N,    :tv_K13d_N,  :tv_K02d_N, 
     :tv_k01_C,     :tv_k02_C,   :tv_k13_C, :tv_k23_C,
     :tv_K01d_C,    :tv_K13d_C,  :tv_K02d_C);
     sort = true,
     kd_C_name=:tv_K01d_C,
     kd_N_name=:tv_K01d_N,
     counternames = byrne_coef_counternames
);
coef_byrne = NamedTuple{byrne_keys_ordered}(coef_byrne)
dd_byrne = Dict(:tv_k01_N  => "kⁿ₀₁",   :tv_k02_N  => "kⁿ₀₂",   :tv_k13_N  => "kⁿ₁₃", :tv_k23_N => "kⁿ₂₃",
            :tv_K01d_N => "K_Dⁿ₀₁", :tv_K13d_N => "K_Dⁿ₁₃", :tv_K02d_N => "K_Dⁿ₀₂", 
            :tv_k01_C  => "kᶜ₀₁",   :tv_k02_C  => "kᶜ₀₂",   :tv_k13_C  => "kᶜ₁₃", :tv_k23_C => "kᶜ₂₃",
            :tv_K01d_C => "K_Dᶜ₀₁", :tv_K13d_C => "K_Dᶜ₁₃", :tv_K02d_C => "K_Dᶜ₀₂")



hm_size_cm     = (17, 22)
hm_size_px     =  hm_size_cm .* dpc 
hm_size_units  =  hm_size_cm .* pt_per_cm
hm_px_per_unit = (hm_size_cm .* dpc ./ (hm_size_cm .* pt_per_cm))[1]

hm = CairoMakie.Figure(size=hm_size_units, fontsize=08)
add_heatmap(hm, coef_bwell,   1, 1       , dd_bwell,    "Scheme 1", 05)
add_heatmap(hm, coef_bhalla,  1, 2       , dd_bhalla,   "Scheme 2", 05)
add_heatmap(hm, coef_shifman, 2, 1       , dd_shifman,  "Scheme 3", 05)
add_heatmap(hm, coef_pepke_m2,2, 2       , dd_pepke_m2, "Scheme 4", 05)
add_heatmap(hm, coef_faas,   3:4, 1:2    , dd_faas,     "Scheme 5", 05)
hm_cbar = add_heatmap(hm, coef_byrne,  5:6, 1:2    , dd_byrne,    "Scheme 6", 05)
CairoMakie.Colorbar(hm[:, 3], hm_cbar, label="Partial Correlation Coefficient", labelsize=10)

CairoMakie.rowgap!(hm.layout, 10)
CairoMakie.colgap!(hm.layout, 10)

CairoMakie.save("Plots/correlations.png", hm, pt_per_unit=1, px_per_unit = hm_px_per_unit, size=hm_size_units)

#######################
### PCORR PLOTS END ###
#######################
#########################
### PARAM PLOTS START ###
#########################

par_bwell         = get_param_lists(fpms_bwell, (:kon_1, :koff_1, :kon_2, :koff_2));
par_bwell_fixed   = get_param_lists(fpms_bwell_fixed, (:kon_1, :koff_1, :kon_2, :koff_2));

par_bhalla        = get_param_lists(fpms_bhalla,       (:kon_1, :koff_1, :kon_2, :koff_2, :kon_3, :koff_3));
par_bhalla_fixed  = get_param_lists(fpms_bhalla_fixed, (:kon_1, :koff_1, :kon_2, :koff_2, :kon_3, :koff_3));

par_shifman        = get_param_lists(fpms_shifman,       (:kon_1, :koff_1, :kon_2, :koff_2, :kon_3, :koff_3, :kon_4, :koff_4));
par_shifman_fixed  = get_param_lists(fpms_shifman_fixed, (:kon_1, :koff_1, :kon_2, :koff_2, :kon_3, :koff_3, :kon_4, :koff_4));

pepke_m2_c_names = (on = :kon_TC, off = :koff_TC)
pepke_m2_n_names = (on = :kon_TN, off = :koff_TN)
pepke_m2_counternames = zip(
    (:kon_TN, :koff_TN, :kon_RN, :koff_RN),
    (:kon_TC, :koff_TC, :kon_RC, :koff_RC))
par_pepke_m2      = get_param_lists(fpms_pepke_m2, (:kon_TN, :koff_TN, 
                                                  :kon_RN, :koff_RN,
                                                  :kon_TC, :koff_TC,
                                                  :kon_RC, :koff_RC);
                                                    sort = true,
                                                    sort_C_names = pepke_m2_c_names,
                                                    sort_N_names = pepke_m2_n_names,
                                                    counternames = pepke_m2_counternames
                                                    );
par_pepke_m2_fixed= get_param_lists(fpms_pepke_m2_fixed, (:kon_TN, :koff_TN, 
                                                  :kon_RN, :koff_RN,
                                                  :kon_TC, :koff_TC,
                                                  :kon_RC, :koff_RC));
faas_c_names = (on = :kon_TC, off = :koff_TC)
faas_n_names = (on = :kon_TN, off = :koff_TN)
faas_counternames = zip(
    (:kon_TN, :koff_TN, :kon_RN, :koff_RN),
    (:kon_TC, :koff_TC, :kon_RC, :koff_RC))
par_faas          = get_param_lists(fpms_faas, (:kon_TN, :koff_TN, 
                                                  :kon_RN, :koff_RN,
                                                  :kon_TC, :koff_TC,
                                                  :kon_RC, :koff_RC);
                                                  sort = true,
                                                  sort_C_names = faas_c_names,
                                                  sort_N_names = faas_n_names,
                                                  counternames = faas_counternames
                                                  );
par_faas_fixed    = get_param_lists(fpms_faas_fixed, (:kon_TN, :koff_TN, 
                                                  :kon_RN, :koff_RN,
                                                  :kon_TC, :koff_TC,
                                                  :kon_RC, :koff_RC));
par_pepke_m1_fixed= get_param_lists(fpms_pepke_m1_fixed, (:kon_TN, :koff_TN, 
                                                  :kon_RN, :koff_RN,
                                                  :kon_TC, :koff_TC,
                                                  :kon_RC, :koff_RC));



byrne_c_names = (on = :k01_C, off = :k10_C)
byrne_n_names = (on = :k01_N, off = :k10_N)
byrne_counternames = zip(
    (:k01_N, :k10_N, :k02_N, :k20_N, :k13_N, :k31_N, :k23_N, :k32_N),
    (:k01_C, :k10_C, :k02_C, :k20_C, :k13_C, :k31_C, :k23_C, :k32_C))
par_byrne         = get_param_lists(fpms_byrne, (:k01_N, :k10_N, 
                                                     :k02_N, :k20_N,
                                                     :k13_N, :k31_N,
                                                     :k23_N, :k32_N,
                                                     :k01_C, :k10_C, 
                                                     :k02_C, :k20_C,
                                                     :k13_C, :k31_C,
                                                     :k23_C, :k32_C);
                                                     sort = true,
                                                     sort_C_names = byrne_c_names,
                                                     sort_N_names = byrne_n_names,
                                                     counternames = byrne_counternames
                                                     );                                                 
par_byrne_fixed   = get_param_lists(fpms_byrne_fixed, (:k01_N, :k10_N, 
                                                     :k02_N, :k20_N,
                                                     :k13_N, :k31_N,
                                                     :k23_N, :k32_N,
                                                     :k01_C, :k10_C, 
                                                     :k02_C, :k20_C,
                                                     :k13_C, :k31_C,
                                                     :k23_C, :k32_C));


#############################################
### ===================================== ###
#############################################


bw_pp_size_cm     = (17, 17)
bw_pp_size_px     =  bw_pp_size_cm .* dpc 
bw_pp_size_units  =  bw_pp_size_cm .* pt_per_cm
bw_pp_px_per_unit = (bw_pp_size_cm .* dpc ./ (bw_pp_size_cm .* pt_per_cm))[1]

bwell_lims = (10^-14, 10^14, 10^-14, 10^14)
###bwell_lims = (10^-6, 10^12, 10^-6, 10^12)
bwell_pairplot = pair_plots(
    [par_bwell, par_bwell_fixed],
###    (; kon_1=L"log$_{10}(k_1)$", koff_1=L"log$_{10}(k_2)$", 
###       kon_2=L"log$_{10}(k_3)$", koff_2=L"log$_{10}(k_4)$"),
    (; kon_1="log₁₀(k₁)", koff_1="log₁₀(k₂)", 
       kon_2="log₁₀(k₃)", koff_2="log₁₀(k₄)"),
    [(our_col, 0.5), (kim_col, 1.0)],
    ["Our rates", "Kim et al."],
    [:circle, :x],
    [0, 0.5];
    lims=bwell_lims,
    fig_w = bw_pp_size_units[1],
    fig_h = bw_pp_size_units[2],
    yvis=true,
    ms=6
)
CairoMakie.Legend(bwell_pairplot[0, :], bwell_pairplot.content[1], "Scheme 1",
    unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
    nbanks=3)
for i in 1:3
    rowsize!(bwell_pairplot.layout, i, Relative(1/3.3))
    colsize!(bwell_pairplot.layout, i, Relative(1/2.9))
end
colgap!(bwell_pairplot.layout, 10)
rowgap!(bwell_pairplot.layout, 1, 5)
rowgap!(bwell_pairplot.layout, 2, -22.5)
rowgap!(bwell_pairplot.layout, 3, -22.5)

###CairoMakie.save("Plots/bwell_pplot.svg", bwell_pairplot, pt_per_unit=1, px_per_unit=bw_pp_px_per_unit,
###    size=bw_pp_size_units)
CairoMakie.save("Plots/bwell_pplot.png", bwell_pairplot, pt_per_unit=1, px_per_unit=bw_pp_px_per_unit,
    size=bw_pp_size_units)

#############################################
### ===================================== ###
#############################################

bh_pp_size_cm     = (17, 22)
bh_pp_size_px     =  bh_pp_size_cm .* dpc 
bh_pp_size_units  =  bh_pp_size_cm .* pt_per_cm
bh_pp_px_per_unit = (bh_pp_size_cm .* dpc ./ (bh_pp_size_cm .* pt_per_cm))[1]

bhalla_lims = (10^-11, 10^11, 10^-11, 10^11)

bhalla_pairplot = pair_plots(
    [par_bhalla, par_bhalla_fixed],
    (; kon_1="log₁₀(k₁)", koff_1="log₁₀(k₂)", 
       kon_2="log₁₀(k₃)", koff_2="log₁₀(k₄)",
       kon_3="log₁₀(k₅)", koff_3="log₁₀(k₆)"),
    [(our_col, 0.5), (bhalla_col, 1.0)],
    ["Our rates", "Bhalla and Iyengar."],
    [:circle, :x],
    [0, 0.5];
    lims=bhalla_lims,
    fig_w = bh_pp_size_units[1],
    fig_h = bh_pp_size_units[2],
    yvis=true,
    ms=6,
    rot_xlab = pi/2
)
CairoMakie.Legend(bhalla_pairplot[0, :], bhalla_pairplot.content[1], "Scheme 2",
    unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
    nbanks=2)
for i in 1:5
    rowsize!(bhalla_pairplot.layout, i, Relative(1/5.3))
    colsize!(bhalla_pairplot.layout, i, Relative(1/4.9))
end
colgap!(bhalla_pairplot.layout, 4)
rowgap!(bhalla_pairplot.layout, -32.5)
rowgap!(bhalla_pairplot.layout, 1, 5)

###CairoMakie.save("Plots/bwell_pplot.svg", bwell_pairplot, pt_per_unit=1, px_per_unit=bw_pp_px_per_unit,
###    size=bw_pp_size_units)
CairoMakie.save("Plots/bhalla_pplot.png", bhalla_pairplot, pt_per_unit=1, px_per_unit=bh_pp_px_per_unit,
    size=bh_pp_size_units)


#############################################
### ===================================== ###
#############################################


sh_pp_size_cm     = (17, 22)
sh_pp_size_px     =  sh_pp_size_cm .* dpc 
sh_pp_size_units  =  sh_pp_size_cm .* pt_per_cm
sh_pp_px_per_unit = (sh_pp_size_cm .* dpc ./ (sh_pp_size_cm .* pt_per_cm))[1]

shifman_lims = (10^-11, 10^11, 10^-11, 10^11)

shifman_pairplot = pair_plots(
    [par_shifman, par_shifman_fixed],
    (; kon_1="log₁₀(k₁)", koff_1="log₁₀(k₂)", 
       kon_2="log₁₀(k₃)", koff_2="log₁₀(k₄)",
       kon_3="log₁₀(k₅)", koff_3="log₁₀(k₆)",
       kon_4="log₁₀(k₇)", koff_4="log₁₀(k₈)"),
    [(our_col, 0.5), (shifman_col, 1.0)],
    ["Our rates", "Shifman et al."],
    [:circle, :x],
    [0, 0.5];
    lims=shifman_lims,
    fig_w = bh_pp_size_units[1],
    fig_h = bh_pp_size_units[2],
    yvis=true,
    ms=6,
    rot_xlab = pi/2
)
CairoMakie.Legend(shifman_pairplot[0, :], shifman_pairplot.content[1], "Scheme 3",
    unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
    nbanks=2)
for i in 1:7
    rowsize!(shifman_pairplot.layout, i, Relative(1/7.3))
    colsize!(shifman_pairplot.layout, i, Relative(1/6.9))
end
colgap!(shifman_pairplot.layout, 4)
rowgap!(shifman_pairplot.layout, -32.5)
rowgap!(shifman_pairplot.layout, 1, 5)

###CairoMakie.save("Plots/bwell_pplot.svg", bwell_pairplot, pt_per_unit=1, px_per_unit=bw_pp_px_per_unit,
###    size=bw_pp_size_units)
CairoMakie.save("Plots/shifman_pplot.png", shifman_pairplot, pt_per_unit=1, px_per_unit=sh_pp_px_per_unit,
    size=sh_pp_size_units)

#############################################
### ===================================== ###
#############################################


ppm2_pp_size_cm     = (17, 22.23)
ppm2_pp_size_px     =  ppm2_pp_size_cm .* dpc 
ppm2_pp_size_units  =  ppm2_pp_size_cm .* pt_per_cm
ppm2_pp_px_per_unit = (ppm2_pp_size_cm .* dpc ./ (ppm2_pp_size_cm .* pt_per_cm))[1]

pepke_lims = (10^-5, 10^10, 10^-5, 10^10)
pepke_m2_pairplot = pair_plots(
    [par_pepke_m2, par_pepke_m2_fixed],
###    (; kon_2C=L"log$_{10}(k_1)$", koff_2C=L"log$_{10}(k_2)$", 
###       kon_2N=L"log$_{10}(k_3)$", koff_2N=L"log$_{10}(k_4)$"),
    (; kon_TC="log₁₀(k₊ᵗᶜ)", koff_TC="log₁₀(k₋ᵗᶜ)", 
       kon_RC="log₁₀(k₊ʳᶜ)", koff_RC="log₁₀(k₋ʳᶜ)", 
       kon_TN="log₁₀(k₊ᵗⁿ)", koff_TN="log₁₀(k₋ᵗⁿ)", 
       kon_RN="log₁₀(k₊ʳⁿ)", koff_RN="log₁₀(k₋ʳⁿ)"),
    [(our_col, 0.5), (pepke_col, 1.0)],
    ["Our rates", "Pepke et al."],
    [:circle, :x],
    [0, 0.5];
    lims = pepke_lims,
    fig_w = ppm2_pp_size_units[1],
    fig_h = ppm2_pp_size_units[2],
    yvis=true,
    ms = 6,
    rot_xlab = pi/2
)
CairoMakie.Legend(pepke_m2_pairplot[0, :], pepke_m2_pairplot.content[1], "Scheme 4",
    unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
    nbanks=3)
for i in 1:7
    rowsize!(pepke_m2_pairplot.layout, i, Relative(1/7.6))
    colsize!(pepke_m2_pairplot.layout, i, Relative(1/6.5))
end
colgap!(pepke_m2_pairplot.layout, 10)
rowgap!(pepke_m2_pairplot.layout, 1, 5)

for i in 2:7
    rowgap!(pepke_m2_pairplot.layout, i, -20)
end

###CairoMakie.save("Plots/pepke_m2_pplot.svg", pepke_m2_pairplot, pt_per_unit=1, px_per_unit=ppm2_pp_px_per_unit,
###    size=ppm2_pp_size_units)
CairoMakie.save("Plots/pepke_m2_pplot.png", pepke_m2_pairplot, pt_per_unit=1, px_per_unit=ppm2_pp_px_per_unit,
    size=ppm2_pp_size_units)


#############################################
### ===================================== ###
#############################################


faas_pp_size_cm     = (17, 22.23)
faas_pp_size_px     =  faas_pp_size_cm .* dpc 
faas_pp_size_units  =  faas_pp_size_cm .* pt_per_cm
faas_pp_px_per_unit = (faas_pp_size_cm .* dpc ./ (faas_pp_size_cm .* pt_per_cm))[1]

faas_lims = (10^-5, 10^10, 10^-5, 10^10)
faas_pairplot = pair_plots(
    [par_faas, par_faas_fixed, par_pepke_m1_fixed],
###    (; kon_TC=L"log$_{10}(k_{on}^{TC})$", koff_TC=L"log$_{10}(k_{off}^{TC})$", 
###       kon_RC=L"log$_{10}(k_{on}^{RC})$", koff_RC=L"log$_{10}(k_{off}^{RC})$", 
###       kon_TN=L"log$_{10}(k_{on}^{TN})$", koff_TN=L"log$_{10}(k_{off}^{TN})$", 
###       kon_RN=L"log$_{10}(k_{on}^{RN})$", koff_RN=L"log$_{10}(k_{off}^{RN})$"),
    (; kon_TC="log₁₀(k₊ᵗᶜ)", koff_TC="log₁₀(k₋ᵗᶜ)", 
       kon_RC="log₁₀(k₊ʳᶜ)", koff_RC="log₁₀(k₋ʳᶜ)", 
       kon_TN="log₁₀(k₊ᵗⁿ)", koff_TN="log₁₀(k₋ᵗⁿ)", 
       kon_RN="log₁₀(k₊ʳⁿ)", koff_RN="log₁₀(k₋ʳⁿ)"),
    [(our_col, 0.5), (faas_col, 1.0), (pepke_col, 1.0)],
    ["Our rates", "Faas et al.", "Pepke et al."],
    [:circle, :x, :utriangle],
    [0, 0.5, 0.5];
    lims = faas_lims,
    fig_w = faas_pp_size_units[1],
    fig_h = faas_pp_size_units[2],
    ms=5,
    rot_xlab = pi/2
)
CairoMakie.Legend(faas_pairplot[0, :], faas_pairplot.content[1], "Scheme 5",
    unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
    nbanks=3)
for i in 1:7
    rowsize!(faas_pairplot.layout, i, Relative(1/7.6))
    colsize!(faas_pairplot.layout, i, Relative(1/6.5))
end
colgap!(faas_pairplot.layout, 7)
rowgap!(faas_pairplot.layout, 1, 5)

for i in 2:7
    rowgap!(faas_pairplot.layout, i, -25)
end

###CairoMakie.save("Plots/faas_pplot.svg", faas_pairplot, pt_per_unit=1, px_per_unit=faas_pp_px_per_unit,
###    size=faas_pp_size_units)
CairoMakie.save("Plots/faas_pplot.png", faas_pairplot, pt_per_unit=1, px_per_unit=faas_pp_px_per_unit,
    size=faas_pp_size_units)


#############################################
### ===================================== ###
#############################################


byrne_pp_size_cm     = (17, 22.23)
byrne_pp_size_px     =  byrne_pp_size_cm .* dpc 
byrne_pp_size_units  =  byrne_pp_size_cm .* pt_per_cm
byrne_pp_px_per_unit = (byrne_pp_size_cm .* dpc ./ (byrne_pp_size_cm .* pt_per_cm))[1]

byrne_lims = (10^-6, 10^11, 10^-6, 10^11)
###byrne_lims = (nothing, nothing, nothing, nothing)
byrne_pairplot = pair_plots(
    [par_byrne, par_byrne_fixed],
###    (; k01_N=L"log$_{10}(k_{01}^{N})$", k10_N=L"log$_{10}(k_{10}^{N})$", 
###       k02_N=L"log$_{10}(k_{02}^{N})$", k20_N=L"log$_{10}(k_{20}^{N})$",
###       k13_N=L"log$_{10}(k_{13}^{N})$", k31_N=L"log$_{10}(k_{31}^{N})$",
###       k23_N=L"log$_{10}(k_{32}^{N})$", k32_N=L"log$_{10}(k_{32}^{N})$",
###       k01_C=L"log$_{10}(k_{01}^{C})$", k10_C=L"log$_{10}(k_{10}^{C})$", 
###       k02_C=L"log$_{10}(k_{02}^{C})$", k20_C=L"log$_{10}(k_{20}^{C})$",
###       k13_C=L"log$_{10}(k_{13}^{C})$", k31_C=L"log$_{10}(k_{31}^{C})$",
###       k23_C=L"log$_{10}(k_{32}^{C})$", k32_C=L"log$_{10}(k_{32}^{C})$",
###    ),
    (; k01_N="k₀₁ⁿ", k10_N="k₁₀ⁿ", 
       k02_N="k₀₂ⁿ", k20_N="k₂₀ⁿ",
       k13_N="k₁₃ⁿ", k31_N="k₃₁ⁿ",
       k23_N="k₂₃ⁿ", k32_N="k₃₂ⁿ",
       k01_C="k₀₁ᶜ", k10_C="k₁₀ᶜ", 
       k02_C="k₀₂ᶜ", k20_C="k₂₀ᶜ",
       k13_C="k₁₃ᶜ", k31_C="k₃₁ᶜ",
       k23_C="k₂₃ᶜ", k32_C="k₃₂ᶜ",
    ),
    [(our_col, 0.5), (byrne_col, 1.0)],
    ["Our rates", "Byrne et al."],
    [:circle, :x],
    [0, 0.5];
    lims = byrne_lims,
    fig_w = byrne_pp_size_units[1],
    fig_h = byrne_pp_size_units[2],
    ms=3,
    rot_xlab = pi/2
)
CairoMakie.Legend(byrne_pairplot[0, :], byrne_pairplot.content[1], "Scheme 6",
    unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
    nbanks=3)

for i in 1:15
    rowsize!(byrne_pairplot.layout, i, Relative(1/15.5))
    colsize!(byrne_pairplot.layout, i, Relative(1/14))
end

colgap!(byrne_pairplot.layout, 5)
rowgap!(byrne_pairplot.layout, 1, 7.5)

for i in 2:15
    rowgap!(byrne_pairplot.layout, i, -32.5)
end

###CairoMakie.save("Plots/byrne_pplot.svg", byrne_pairplot, pt_per_unit=1, px_per_unit=byrne_pp_px_per_unit,
###    size=byrne_pp_size_units)
CairoMakie.save("Plots/byrne_pplot.png", byrne_pairplot, pt_per_unit=1, px_per_unit=byrne_pp_px_per_unit,
    size=byrne_pp_size_units)
########################
### PARAM PLOT START ###
########################
####################################
### TOTAL CALCIUM EQ PLOTS START ###
####################################
###tot_ca_range    = 10 .^ LinRange(-8, 0, 100);
###dyn_eq_pop   = gen_dyn_eq_data(tot_ca_range);
###
###f_crouch = XLSX.readxlsx("Crouch_Klee_1980.xlsx");
###X_crouch = f_crouch["Crouch_Klee_1980"]["A"][1:end];
###Y_crouch = f_crouch["Crouch_Klee_1980"]["B"][1:end];
###
###bwell_CaM_dyn_eqs_tups = get_multiple_dyns(fpms_bwell, dyn_eq_pop, get_bwell_dyns_last);
###bwell_fixed_CaM_dyn_eqs_tups = get_multiple_dyns(fpms_bwell_fixed, dyn_eq_pop, get_bwell_dyns_last);
###bwell_CN_CaM_dyn_eqs_tups = get_multiple_dyns(fpms_bwell_CN, dyn_eq_pop, get_bwell_CN_dyns_last);
###
###p_eqs_scheme1 = plot_CaM_eqs(tot_ca_range, bwell_CaM_dyn_eqs_tups; title="Scheme 1", color="red", label="Ours", xlabel="[Ca²⁺ total]",
###    new_plot=true, legloc=:topleft, alpha=0.1, xlim=[1e-8, 1e-3], xscale=:log10)
###p_eqs_scheme1 = plot_CaM_eqs(tot_ca_range, bwell_CN_CaM_dyn_eqs_tups; title="Scheme 1", color="green", label="Ours + lobes", 
###    new_plot=false, legloc=:bottomright, prev_plot=p_eqs_scheme1, alpha=0.1, xlim=[1e-8, 1e-3], xscale=:log10)
###p_eqs_scheme1 = plot_CaM_eqs(tot_ca_range, bwell_fixed_CaM_dyn_eqs_tups; title="Scheme 1", color="blue", label="Published", 
###    new_plot=false, legloc=:bottomright, prev_plot=p_eqs_scheme1, alpha=0.1, xlim=[1e-8, 1e-3], xscale=:log10)
###
###p_eqs_scheme1 = Plots.scatter!(10 .^ X_crouch, Y_crouch, marker=:v, label="Crouch and Klee (1980) data", color="black")
###
###
###
###faas_CaM_dyn_eqs_tups = get_multiple_dyns(fpms_faas, dyn_eq_pop, get_faas_dyns_last);
###faas_fixed_CaM_dyn_eqs_tups = get_multiple_dyns(fpms_faas_fixed, dyn_eq_pop, get_faas_dyns_last);
###
###p_eqs_scheme2 = plot_CaM_eqs(tot_ca_range, faas_CaM_dyn_eqs_tups; title="Scheme 2", color="red", label="Ours", xlabel="[Ca²⁺ total]",
###    new_plot=true, legloc=:topleft, alpha=0.1, xlim=[1e-8, 1e-3], xscale=:log10)
###p_eqs_scheme2 = plot_CaM_eqs(tot_ca_range, faas_fixed_CaM_dyn_eqs_tups; title="Scheme 2", color="blue", label="Published", 
###    new_plot=false, legloc=:bottomright, prev_plot=p_eqs_scheme2, alpha=0.1, xlim=[1e-8, 1e-3], xscale=:log10)
###
###p_eqs_scheme2 = Plots.scatter!(10 .^ X_crouch, Y_crouch, marker=:v, label="Crouch and Klee (1980) data", color="black")
###
###
###
###byrne_CaM_dyn_eqs_tups = get_multiple_dyns(fpms_byrne, dyn_eq_pop, get_byrne_dyns_last);
###byrne_fixed_CaM_dyn_eqs_tups = get_multiple_dyns(fpms_byrne_fixed, dyn_eq_pop, get_byrne_dyns_last);
###
###p_eqs_scheme3 = plot_CaM_eqs(tot_ca_range, byrne_CaM_dyn_eqs_tups; title="Scheme 3", color="red", label="Ours", xlabel="[Ca²⁺ total]",
###    new_plot=true, legloc=:topleft, alpha=0.1, xlim=[1e-8, 1e-3], xscale=:log10)
###p_eqs_scheme2 = plot_CaM_eqs(tot_ca_range, byrne_fixed_CaM_dyn_eqs_tups; title="Scheme 3", color="blue", label="Published", 
###    new_plot=false, legloc=:bottomright, prev_plot=p_eqs_scheme3, alpha=0.1, xlim=[1e-8, 1e-3], xscale=:log10)
###
###p_eqs_scheme3 = Plots.scatter!(10 .^ X_crouch, Y_crouch, marker=:v, label="Crouch and Klee (1980) data", color="black")
###
###
###
########################
### DATA PLOTS START ###
########################
pexp_pp_size_cm     = (17, 7)
pexp_pp_size_px     =  pexp_pp_size_cm .* dpc 
pexp_pp_size_units  =  pexp_pp_size_cm .* pt_per_cm
pexp_pp_px_per_unit = (pexp_pp_size_cm .* dpc ./ (pexp_pp_size_cm .* pt_per_cm))[1]

p_data = plot_exp_data([A,B,G]; figsize=pexp_pp_size_units, lw=1, ylim=[0, 22.5])
##CairoMakie.save("Plots/data_plots.svg", p_data, pt_per_unit=1, px_per_unit=pexp_pp_px_per_unit,
##    size=pexp_pp_size_units)
CairoMakie.save("Plots/data_plots.png", p_data, pt_per_unit=1, px_per_unit=pexp_pp_px_per_unit,
size=pexp_pp_size_units)


psub_pp_size_cm     = (17, 22)
psub_pp_size_px     =  psub_pp_size_cm .* dpc 
psub_pp_size_units  =  psub_pp_size_cm .* pt_per_cm
psub_pp_px_per_unit = (psub_pp_size_cm .* dpc ./ (psub_pp_size_cm .* pt_per_cm))[1]

groups = [A, B, C, D, E, F, G];
sub_groups = [subsample_start(i, sub_range) for i in groups];
p_sub_data = plot_subsampled_data(groups, sub_groups; figsize=psub_pp_size_units, lw=1, ylim=[0, 22.5])
##CairoMakie.save("Plots/data_plots.svg", p_data, pt_per_unit=1, px_per_unit=pexp_pp_px_per_unit,
##    size=pexp_pp_size_units)
CairoMakie.save("Plots/subsample_plot.png", p_sub_data, pt_per_unit=1, px_per_unit=psub_pp_px_per_unit,
size=psub_pp_size_units)


pshif_pp_size_cm     = (17, 7)
pshif_pp_size_px     =  pshif_pp_size_cm .* dpc 
pshif_pp_size_units  =  pshif_pp_size_cm .* pt_per_cm
pshif_pp_px_per_unit = (pshif_pp_size_cm .* dpc ./ (pshif_pp_size_cm .* pt_per_cm))[1]

p_shif = CairoMakie.Figure(size=pshif_pp_size_units, fontsize=8)
ax_shif = CairoMakie.Axis(p_shif[1,1], 
    ylabel="# Bound Ca²⁺ per CaM", 
    xlabel="Free Ca²⁺ (M)", 
 ###   limits = (-3, 43, ylim[1], ylim[2]), 
    yticklabelsize=8, 
    xticklabelsize=8,
)
CairoMakie.scatter!(ax_shif, X_shifman, Y_shifman, color=(:black, 0.3), marker=:utriangle)

###CairoMakie.save("Plots/data_shifman.svg", p_shif, pt_per_unit=1, px_per_unit=pshif_pp_px_per_unit,
###    size=pshif_pp_size_units)
CairoMakie.save("Plots/data_shifman.png", p_shif, pt_per_unit=1, px_per_unit=pshif_pp_px_per_unit,
    size=pshif_pp_size_units)

###############################
### BYRNE C/N LOBE ANALYSES ###
###############################

c_mat, t_mat, catot = byrne_N_C_lobe_binding_comp(fpms_byrne_fixed);

bb_size_cm     = (17, 17)
bb_size_px     =  bb_size_cm .* dpc 
bb_size_units  =  bb_size_cm .* pt_per_cm
bb_px_per_unit = (bb_size_cm .* dpc ./ (bb_size_cm .* pt_per_cm))[1]


linecol = (:black, 0.025)
bb = CairoMakie.Figure(size=bb_size_units, fontsize=8)

bb_ax11 = CairoMakie.Axis(bb[1,1], #ylabel="N0 / CaM_total", 
    title = "k₀₁ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₁ᴰ = 33×10⁻⁶ M \n
k₀₂ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₂ᴰ = 229.88×10⁻⁶ M",
    xscale=log10,
    limits = (nothing, nothing, 0, 1.1),
    xticks = [0.16, 1, 5, 35],
    xticksvisible = false,
    xticklabelsvisible = false
    )
for (x, y) in zip(t_mat.n0, c_mat.n0)  
    CairoMakie.lines!(bb_ax11,
    x, y, color=linecol )
end

bb_ax12 = CairoMakie.Axis(bb[1,2], #ylabel="N1 / CaM_total", 
    title = "k₀₁ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₁ᴰ = 33×10⁻⁶ M \n
k₁₃ = 5.072×10⁵ M⁻¹ms⁻¹ \n K₁₃ᴰ = 3.45×10⁻⁶ M",
    xscale=log10,
    limits = (nothing, nothing, 0, 1.1),
    xticks = [0.16, 1, 5, 35],
    xticksvisible = false,
    xticklabelsvisible = false,
    yticksvisible = false,
    yticklabelsvisible = false
    )
for (x, y) in zip(t_mat.n1, c_mat.n1)  
    CairoMakie.lines!(bb_ax12,
    x, y, color=linecol )
end

bb_ax13 = CairoMakie.Axis(bb[1,3], #ylabel="N2 / CaM_total", 
    title = "k₀₂ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₂ᴰ = 229.88×10⁻⁶ M \n
k₂₃ = 5×10⁵ M⁻¹ms⁻¹ \n K₂₃ᴰ = 0.5×10⁻⁶ M",
    xscale=log10,
    limits = (nothing, nothing, 0, 1.1),
    xticks = [0.16, 1, 5, 35],
    xticksvisible = false,
    xticklabelsvisible = false,
    yticksvisible = false,
    yticklabelsvisible = false
    )
for (x, y) in zip(t_mat.n2, c_mat.n2)  
    CairoMakie.lines!(bb_ax13,
    x, y, color=linecol)
end

bb_ax14 = CairoMakie.Axis(bb[1,4], #ylabel="N3 / CaM_total", 
    title = "k₁₃ = 5.072×10⁵ M⁻¹ms⁻¹ \n K₁₃ᴰ=3.45×10⁻⁶ M \n
k₂₃=5×10⁵ M⁻¹ms⁻¹ \n K₂₃ᴰ=0.5×10⁻⁶ M",
    xscale=log10,
    limits = (nothing, nothing, 0, 1.1),
    xticks = [0.16, 1, 5, 35],
    xticksvisible = false,
    xticklabelsvisible = false,
    yticksvisible = false,
    yticklabelsvisible = false
    )
for (x, y) in zip(t_mat.n3, c_mat.n3)  
    CairoMakie.lines!(bb_ax14,
    x, y, color=linecol)
end



bb_ax21 = CairoMakie.Axis(bb[2,1], #ylabel="C0 / CaM_total", 
    title = "k₀₁ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₁ᴰ = 18.5×10⁻⁶ M \n
k₀₂ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₂ᴰ = 116×10⁻⁶ M",
    xscale=log10,
    limits = (nothing, nothing, 0, 1.1),
    xlabel="log₁₀(t (ms))",
    xticks = [0.16, 1, 5, 35]
    )
for (x, y) in zip(t_mat.c0, c_mat.c0)  
    CairoMakie.lines!(bb_ax21,
    x, y, color=linecol )
end

bb_ax22 = CairoMakie.Axis(bb[2,2], #ylabel="C1 / CaM_total", 
    title = "k₀₁ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₁ᴰ = 18.5×10⁻⁶ M \n
k₁₃ = 3.71×10³ M⁻¹ms⁻¹ \n K₁₃ᴰ = 0.38×10⁻⁶ M",
    xscale=log10,
    limits = (nothing, nothing, 0, 1.1),
    xlabel="log₁₀(t (ms))",
    xticks = [0.16, 1, 5, 35],
    yticksvisible = false,
    yticklabelsvisible = false
)
for (x, y) in zip(t_mat.c1, c_mat.c1)  
    CairoMakie.lines!(bb_ax22,
    x, y, color=linecol)
end

bb_ax23 = CairoMakie.Axis(bb[2,3], #ylabel="C2 / CaM_total", 
    title = "k₀₂ = 2.75×10⁵ M⁻¹ms⁻¹ \n K₀₂ᴰ = 116×10⁻⁶ M \n
k₂₃ = 1.18×10⁵ M⁻¹ms⁻¹ \n K₂₃ᴰ = 0.06×10⁻⁶ M",
    xscale=log10,
    limits = (nothing, nothing, 0, 1.1),
    xticks = [0.16, 1, 5, 35],
    yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="log₁₀(t (ms))")
for (x, y) in zip(t_mat.c2, c_mat.c2)  
    CairoMakie.lines!(bb_ax23,
    x, y, color=linecol)
end

bb_ax24 = CairoMakie.Axis(bb[2,4], #ylabel="C3 / CaM_total", 
    title = "k₁₃ = 3.71×10³ M⁻¹ms⁻¹ \n K₁₃ᴰ = 0.38×10⁻⁶ M \n
k₂₃ = 1.18×10⁵ M⁻¹ms⁻¹ \n K₂₃ᴰ = 0.06×10⁻⁶ M",
    xscale=log10,
    limits = (nothing, nothing, 0, 1.1),
    xticks = [0.16, 1, 5, 35],
    yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="log₁₀(t (ms))")
for (x, y) in zip(t_mat.c3, c_mat.c3)  
    CairoMakie.lines!(bb_ax24,
    x, y, color=linecol)
end

CairoMakie.Label(bb[1,0], "N lobe (fraction)", rotation=pi/2)
CairoMakie.Label(bb[2,0], "C lobe (fraction)", rotation=pi/2)

CairoMakie.Label(bb[0,1], "CaM + 0Ca")
CairoMakie.Label(bb[0,2], "CaM + 1Ca \n First site")
CairoMakie.Label(bb[0,3], "CaM + 1Ca \n Second site")
CairoMakie.Label(bb[0,4], "CaM + 2Ca")

for i in 1:2
    rowsize!(bb.layout, i, Relative(1/2))
end
for i in 1:4
    colsize!(bb.layout, i, Relative(1/4))
end
colgap!(bb.layout, 1, 5)
rowgap!(bb.layout, 1, 5)
###CairoMakie.save("Plots/byrne_lobe_analysis.svg", bb, pt_per_unit=1, px_per_unit=bb_px_per_unit,
###    size=bb_size_units)
CairoMakie.save("Plots/byrne_lobe_analysis.png", bb, pt_per_unit=1, px_per_unit=bb_px_per_unit,
    size=bb_size_units)

###############################
### BYRNE C/N LOBE ANALYSES ###
###############################

############################
### PEPKE M1/M2 ANALYSES ###
############################


###plot_training_by_condition([fpm_12, fpm_22, fpm_dif_rates_2], 
###                           ["red", "blue", "green", "orange"]; 
###                     names=["Scheme", "Faas et al", "Byrne et al."], add_legend=true)
###plot_unseen_by_condition([fpms_bwell[1], fpms_bwell_fixed[1]],
###                          val_pop,
###                         ["red", "blue"]; 
###                   names=["Scheme 1", "Scheme 1 fixed rates"], add_legend=true)
###
###
###plot_fens([fpms_bwell, fpms_bwell_fixed, fpms_faas, fpms_faas_fixed, fpms_byrne, fpms_byrne_fixed], 
###    E[[1, 4, 7, 10, 14]], 
###    ["red", "blue"]; 
###    add_legend=true
###)












###all_pre = [fpm.model.pre(coef(fpm), empirical_bayes(fpm)[i], pop[i]).ret for i in 1:length(fpm.data)];
######
######
######preds = predict(fpm);
###sims = simobs(fpm_dif_rates);

i_p = (; Ω=1, σ=10)
model = Byrne_scheme
sims = simobs(model, train_pop, merge(i_p, Blackwell_params));
sims = simobs(model, full_pop[1], init_params(model));
sims = simobs(model, full_pop[1], sample_params(model));
sims = simobs(model, train_pop[1], coef(fpm_Pepke_TR));

######
######
######Plots.plot()
######for p in preds
######    Plots.plot!(p.subject.time, p.pred.F_F0, label=false, ls=:dash, color="red" )
######    Plots.plot!(p.subject.time, p.subject.observations.F_F0, label=false, color="black")
######end
######Plots.plot!()
######
######
######Plots.plot()
######for p in sims 
######    Plots.plot!(p.subject.time, p.observations.F_F0, label=false, ls=:dash, color="red" )
######    Plots.plot!(p.subject.time, p.subject.observations.F_F0, label=false, color="black")
######end
######Plots.plot!()
######
###
aa = sims
aa = ss[3]
p1 = Plots.plot()
Plots.plot!(aa.time, aa.dynamics.DMn_s,     label="DMn_s")
Plots.plot!(aa.time, aa.dynamics.CaDMn_s,   label="CaDMn_s")

p2 = Plots.plot()
Plots.plot!(aa.time, aa.dynamics.DMn_f,     label="DMn_f")
Plots.plot!(aa.time, aa.dynamics.CaDMn_f,   label="CaDMn_f")

p3 = Plots.plot()
Plots.plot!(aa.time, aa.dynamics.PP,        label="PP")
Plots.plot!(aa.time, aa.dynamics.CaPP,      label="CaPP")

p4 = Plots.plot()
Plots.plot!(aa.time, aa.dynamics.OGB5,      label="OGB5")
Plots.plot!(aa.time, aa.dynamics.CaOGB5,    label="CaOGB5")

p5 = Plots.plot()
Plots.plot!(aa.time, aa.dynamics.N0,   label="N0")
Plots.plot!(aa.time, aa.dynamics.N1,   label="N1")
Plots.plot!(aa.time, aa.dynamics.N2,   label="N2")
Plots.plot!(aa.time, aa.dynamics.N3,   label="N3")
Plots.plot!(aa.time, aa.dynamics.C0,   label="C0")
Plots.plot!(aa.time, aa.dynamics.C1,   label="C1")
Plots.plot!(aa.time, aa.dynamics.C2,   label="C2")
Plots.plot!(aa.time, aa.dynamics.C3,   label="C3")

###p5 = Plots.plot()
###Plots.plot!(aa.time, aa.dynamics.T0,    label="T0")
###Plots.plot!(aa.time, aa.dynamics.T1C,   label="T1C")
###Plots.plot!(aa.time, aa.dynamics.T2C,   label="T2C")
###Plots.plot!(aa.time, aa.dynamics.T1N,   label="T1N")
###Plots.plot!(aa.time, aa.dynamics.T2N,   label="T2N")
###Plots.plot!(aa.time, aa.dynamics.T1N1C, label="T1N1C")
###Plots.plot!(aa.time, aa.dynamics.T1N2C, label="T1N2C")
###Plots.plot!(aa.time, aa.dynamics.T2N1C, label="T2N1C")
###Plots.plot!(aa.time, aa.dynamics.T4,    label="T4")
###
###
###p6 = Plots.plot()
###Plots.plot!(aa.time, aa.dynamics.R0,    label="R0")
###Plots.plot!(aa.time, aa.dynamics.R1C,   label="R1C")
###Plots.plot!(aa.time, aa.dynamics.R2C,   label="R2C")
###Plots.plot!(aa.time, aa.dynamics.R1N,   label="R1N")
###Plots.plot!(aa.time, aa.dynamics.R2N,   label="R2N")
###Plots.plot!(aa.time, aa.dynamics.R1N1C, label="R1N1C")
###Plots.plot!(aa.time, aa.dynamics.R1N2C, label="R1N2C")
###Plots.plot!(aa.time, aa.dynamics.R2N1C, label="R2N1C")
###Plots.plot!(aa.time, aa.dynamics.R4,    label="R4")
###
###
###p8 = Plots.plot()
###tot_CaM = aa.dynamics.T0 + aa.dynamics.T1C + aa.dynamics.T1N + 
###          aa.dynamics.T2C + aa.dynamics.T1N1C + aa.dynamics.T2N +
###          aa.dynamics.T1N2C + aa.dynamics.T2N1C + aa.dynamics.T4 +
###          aa.dynamics.R0 + aa.dynamics.R1C + aa.dynamics.R1N + 
###          aa.dynamics.R2C + aa.dynamics.R1N1C + aa.dynamics.R2N +
###          aa.dynamics.R1N2C + aa.dynamics.R2N1C + aa.dynamics.R4
###Plots.plot!(aa.time, tot_CaM, label="CaM_T")


###p5 = Plots.plot()
###Plots.plot!(aa.time, aa.dynamics.T0,   label="T0")
###Plots.plot!(aa.time, aa.dynamics.T_A,   label="T_A")
###Plots.plot!(aa.time, aa.dynamics.T_B,   label="T_B")
###Plots.plot!(aa.time, aa.dynamics.T_C,   label="T_C")
###Plots.plot!(aa.time, aa.dynamics.T_D,   label="T_D")
###
###p6 = Plots.plot()
###Plots.plot!(aa.time, aa.dynamics.T_AB,    label="T_AB")
###Plots.plot!(aa.time, aa.dynamics.T_AC,   label="T_AC")
###Plots.plot!(aa.time, aa.dynamics.T_AD,   label="T_AD")
###Plots.plot!(aa.time, aa.dynamics.T_BC,   label="T_BC")
###Plots.plot!(aa.time, aa.dynamics.T_BD,   label="T_BD")
###Plots.plot!(aa.time, aa.dynamics.T_CD,   label="T_CD")
###
###p7 = Plots.plot()
###Plots.plot!(aa.time, aa.dynamics.T_ABC,    label="T_ABC")
###Plots.plot!(aa.time, aa.dynamics.T_ABD,   label="T_ABD")
###Plots.plot!(aa.time, aa.dynamics.T_ACD,   label="T_ACD")
###Plots.plot!(aa.time, aa.dynamics.T_BCD,   label="T_BCD")
###Plots.plot!(aa.time, aa.dynamics.T_ABCD,   label="T_ABCD")
###
###p8 = Plots.plot()
###tot_CaM = aa.dynamics.T0 + aa.dynamics.T_A +aa.dynamics.T_B +aa.dynamics.T_C +aa.dynamics.T_D +aa.dynamics.T_AB +aa.dynamics.T_AC +aa.dynamics.T_AD +aa.dynamics.T_BC +aa.dynamics.T_BD +aa.dynamics.T_CD +aa.dynamics.T_ABC +aa.dynamics.T_ABD +aa.dynamics.T_ACD +aa.dynamics.T_BCD +aa.dynamics.T_ABCD
###Plots.plot!(aa.time, tot_CaM, label="CaM_T")
###p5 = Plots.plot()
###Plots.plot!(aa.time, aa.dynamics.CaM,      label="CaM")
###Plots.plot!(aa.time, aa.dynamics.CaM2Ca,   label="CaM2Ca")
###Plots.plot!(aa.time, aa.dynamics.CaM4Ca,   label="CaM4Ca")

###p5 = Plots.plot()
###Plots.plot!(aa.time, aa.dynamics.NtNt,      label="NtNt")
###Plots.plot!(aa.time, aa.dynamics.NtNr,      label="NtNr")
###Plots.plot!(aa.time, aa.dynamics.NrNr,      label="NrNr")
###
###p6 = Plots.plot()
###Plots.plot!(aa.time, aa.dynamics.CtCt,      label="CtCt")
###Plots.plot!(aa.time, aa.dynamics.CtCr,      label="CtCr")
###Plots.plot!(aa.time, aa.dynamics.CrCr,      label="CrCr")

p7 = Plots.plot()
Plots.plot!(aa.time, aa.dynamics.Ca,        label="Ca")

###Plots.plot(p1, p2, p3, p4, p5, p7)
Plots.plot(p1, p2, p3, p4, p5, p6, p7)