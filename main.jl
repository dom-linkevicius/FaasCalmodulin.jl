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



include("models.jl");
include("SymbolicsForwardDiffExt.jl")
include("process_data.jl");
include("convenience.jl");
include("symbolic_manipulations.jl");
include("models_with_NN.jl");
include("eq_comparisons.jl");
include("dyn_comparisons.jl");


f_xlsx = XLSX.readxlsx("dataset.xlsx");
full_pop  = fetch_data(f_xlsx);
###og_pop    = fetch_data(f_xlsx);
###og_pop    = subsample_start(og_pop, 1:10:230)

A = full_pop[1:13];
B = full_pop[14:24];
C = full_pop[25:35];
D = full_pop[36:51];
E = full_pop[52:65];
F = full_pop[66:80];
G = full_pop[81:94];


#####train_pop = [A[1:7  ]; B[1:7  ]; C[1:7  ]; D[1:7  ]; E[1:7  ]; F[1:7  ]; G[1:7]];
#####val_pop   = [A[8:10 ]; B[8:9  ]; C[8:9  ]; D[8:12 ]; E[8:11 ]; F[8:11 ]; G[8:11]];
#####test_pop  = [A[11:13]; B[10:11]; C[10:11]; D[13:16]; E[12:14]; F[12:15]; G[12:14]];
###train_pop = [A[1:7  ]; B[1:7  ]; D[1:7  ]; E[1:7  ]; F[1:7  ]; G[1:7]];
###val_pop   = [A[8:10 ]; B[8:9  ]; D[8:12 ]; E[8:11 ]; F[8:11 ]; G[8:11]];
###test_pop  = [A[11:13]; B[10:11]; D[13:16]; E[12:14]; F[12:15]; G[12:14]];


n_models = 20
seeds = 1:n_models
sub_range = nothing  #1:30:231

ss = subsample_start(full_pop, sub_range)
A_ss = ss[1:13];
B_ss = ss[14:24];
C_ss = ss[25:35];
D_ss = ss[36:51];
E_ss = ss[52:65];
F_ss = ss[66:80];
G_ss = ss[81:94];

fpms_bwell, p_bwell_train, p_bwell_val, val_idxs_bwell, test_idxs_bwell = train_plot("Blackwell", 
    full_pop, 
    n_models,
    seeds;
    sub_range = sub_range,
);
###serialize("Data/fpms_bwell.jls",fpms_bwell)
###serialize("Data/val_idxs_bwell.jls",val_idxs_bwell)
###serialize("Data/test_idxs_bwell.jls",test_idxs_bwell)

fpms_bwell = deserialize("Data/fpms_bwell.jls")
val_idxs_bwell = deserialize("Data/val_idxs_bwell.jls")
test_idxs_bwell = deserialize("Data/test_idxs_bwell.jls")


fpms_bwell_fixed, p_bwell_train_fixed, p_bwell_val_fixed, val_idxs_bwell_fixed, test_idxs_bwell_fixed = train_plot("Blackwell_fixed", 
    full_pop, 
    n_models,
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_bwell_fixed.jls",fpms_bwell_fixed)
fpms_bwell_fixed = deserialize("Data/fpms_bwell_fixed.jls")


fpms_pepke_m2, p_pepke_m2_train, p_pepke_m2_val, val_idxs_pepke_m2, test_idxs_pepke_m2 = train_plot("Pepke_m2", 
    full_pop, 
    n_models, 
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_pepke_m2.jls",fpms_pepke_m2)
fpms_pepke_m2 = deserialize("Data/fpms_pepke_m2.jls")


fpms_pepke_m2_fixed, p_pepke_m2_train_fixed, p_pepke_m2_val_fixed, val_idxs_pepke_m2_fixed, test_idxs_pepke_m2_fixed = train_plot("Pepke_m2_fixed", 
    full_pop, 
    n_models, 
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_pepke_m2_fixed.jls",fpms_pepke_m2_fixed)
fpms_pepke_m2_fixed = deserialize("Data/fpms_pepke_m2_fixed.jls")



fpms_faas, p_faas_train, p_faas_val, val_idxs_faas, test_idxs_faas = train_plot("Faas", 
    full_pop, 
    n_models,
    seeds;
    sub_range = sub_range
);
###serialize("Data/fpms_faas.jls",fpms_faas)
fpms_faas = deserialize("Data/fpms_faas.jls")

fpms_faas_fixed, p_faas_train_fixed, p_faas_val_fixed, val_idxs_faas_fixed, test_idxs_faas_fixed = train_plot("Faas_fixed", 
    full_pop, 
    n_models,
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_faas_fixed.jls",fpms_faas_fixed)
fpms_faas_fixed = deserialize("Data/fpms_faas_fixed.jls")

fpms_pepke_m1_fixed, p_pepke_m1_train_fixed, p_pepke_m1_val_fixed, val_idxs_pepke_m1_fixed, test_idxs_pepke_m1_fixed = train_plot("Pepke_m1_fixed", 
    full_pop, 
    n_models,
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_pepke_m1_fixed.jls",fpms_pepke_m1_fixed)
fpms_pepke_m1_fixed = deserialize("Data/fpms_pepke_m1_fixed.jls")


fpms_byrne, p_byrne_train, p_byrne_val, val_idxs_byrne, test_idxs_byrne = train_plot("Byrne", 
    full_pop,
    n_models,
    seeds; 
    sub_range = sub_range,
);
###serialize("Data/fpms_byrne.jls",fpms_byrne)
fpms_byrne = deserialize("Data/fpms_byrne.jls")

fpms_byrne_fixed, p_byrne_train_fixed, p_byrne_val_fixed, val_idxs_byrne_fixed, test_idxs_byrne_fixed = train_plot("Byrne_fixed", 
    full_pop, 
    n_models,
    seeds; 
    sub_range = sub_range
);
###serialize("Data/fpms_byrne_fixed.jls",fpms_byrne_fixed)
fpms_byrne_fixed = deserialize("Data/fpms_byrne_fixed.jls")


###full_pop[1].observations.F_F0
###
###fpm_ms = msfit(
###    Pepke_m2_scheme,
###    fpms_bwell[1].data,
###    init_params(Pepke_m2_scheme),
###    MAP(FOCE()),
###    100,
###    5e-16, #or standard deviation of shooting penalty hyper prior if below argument is not nothing
###    :L1; #:L1 or :L2, :L1 should be highly favoured because it performed much better in experiments
###    hyper_prior_type = nothing, #nothing/:L1/:L2,
###    return_with_shooting = false, #true/false, if false, does a 0 iteration fit and returns a FittePumasModel with the original model
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
kim_col     = RGB(247/255,072/255,165/255);
pepke_col   = RGB(053/255,155/255,115/255);
faas_col    = RGB(213/255,094/255,000/255);
byrne_col   = RGB(240/255,228/255,066/255);


###########################
### AIC BOX PLOTS START ###
###########################

all_fpms = [fpms_bwell;    fpms_bwell_fixed;
            fpms_pepke_m2; fpms_pepke_m2_fixed;
            fpms_faas;     fpms_faas_fixed; fpms_pepke_m1_fixed;
            fpms_byrne;    fpms_byrne_fixed];


val_pops = [idx_to_val_vec(full_pop, idxs) for idxs in val_idxs_bwell];
all_val_pops = repeat(val_pops, 9);


all_names = [repeat(["Scheme 1 + our rates"], length(val_pops)); 
             repeat(["Scheme 1 + Kim et al."], length(val_pops));
             repeat(["Scheme 2 + our rates"], length(val_pops));
             repeat(["Scheme 2 + Pepke et al."], length(val_pops));
             repeat(["Scheme 3 + our rates"], length(val_pops));
             repeat(["Scheme 3 + Faas et al."], length(val_pops));
             repeat(["Scheme 3 + Pepke et al."], length(val_pops));
             repeat(["Scheme 4 + our rates"], length(val_pops));
             repeat(["Scheme 4 + Byrne et al."], length(val_pops))];

aic_pairs = zip(all_fpms, all_val_pops, all_names);

x, y = aic_boxplots(aic_pairs);
x_makie = sort(repeat(1:length(unique(x)), Int(length(x)/length(unique(x)))));
y_makie = [i for i in y];
boxplot_cols = [repeat([our_col], 20); 
        repeat([kim_col], 20);
        repeat([our_col], 20);
        repeat([pepke_col], 20);
        repeat([our_col], 20);
        repeat([faas_col], 20);
        repeat([pepke_col], 20);
        repeat([our_col], 20);
        repeat([byrne_col], 20)];

boxplot_size_cm     = (17, 15)
boxplot_size_px     =  boxplot_size_cm .* dpc 
boxplot_size_units  =  boxplot_size_cm .* pt_per_cm
boxplot_px_per_unit = (boxplot_size_cm .* dpc ./ (boxplot_size_cm .* pt_per_cm))[1]

boxplot = CairoMakie.Figure(size=boxplot_size_units, fontsize=8)
ax_boxplot = CairoMakie.Axis(boxplot[1,1], 
    xticks=(unique(x_makie), unique(x)),
    xlabel="Scheme + parameters",
    ylabel="AIC",
    xticklabelrotation=pi/4,
    xticklabelpad=10,
    xlabelpadding=30,
    xticklabelsize=8,
    yticklabelsize=8,
    xgridvisible=false,
    )
CairoMakie.boxplot!(ax_boxplot, x_makie, y_makie, show_notch=false, color=boxplot_cols)
CairoMakie.save("Plots/AIC_boxplot.svg", boxplot, pt_per_unit=1, 
    px_per_unit = boxplot_px_per_unit, size=boxplot_size_units)
    

#########################
### AIC BOX PLOTS END ###
#########################
#############################
### MSE HISTS PLOTS START ###
#############################
bwell_mses           = [mean(unseen_rmse(i, full_pop[j])) for (i,j) in zip(fpms_bwell, val_idxs_bwell)]
bwell_fixed_mses     = [mean(unseen_rmse(i, full_pop[j])) for (i,j) in zip(fpms_bwell_fixed, val_idxs_bwell)]

pepke_m2_mses        = [mean(unseen_rmse(i, full_pop[j])) for (i,j) in zip(fpms_pepke_m2, val_idxs_bwell)]
pepke_m2_fixed_mses  = [mean(unseen_rmse(i, full_pop[j])) for (i,j) in zip(fpms_pepke_m2_fixed, val_idxs_bwell)]

faas_mses            = [mean(unseen_rmse(i, full_pop[j])) for (i,j) in zip(fpms_faas, val_idxs_bwell)]
faas_fixed_mses      = [mean(unseen_rmse(i, full_pop[j])) for (i,j) in zip(fpms_faas_fixed, val_idxs_bwell)]
pepke_m1_fixed_mses  = [mean(unseen_rmse(i, full_pop[j])) for (i,j) in zip(fpms_pepke_m1_fixed, val_idxs_bwell)]

byrne_mses           = [mean(unseen_rmse(i, full_pop[j])) for (i,j) in zip(fpms_byrne, val_idxs_bwell)]
byrne_fixed_mses     = [mean(unseen_rmse(i, full_pop[j])) for (i,j) in zip(fpms_byrne_fixed, val_idxs_bwell)]



hist_size_cm     = (17, 10)
hist_size_px     =  hist_size_cm .* dpc 
hist_size_units  =  hist_size_cm .* pt_per_cm
hist_px_per_unit = (hist_size_cm .* dpc ./ (hist_size_cm .* pt_per_cm))[1]


hist = CairoMakie.Figure(size=hist_size_units, fontsize=8)


ax11 = CairoMakie.Axis(hist[2,1])
h_111 = CairoMakie.hist!(ax11, bwell_mses,       bins=15, color=(our_col, 0.5), strokewidth=1, strokecolor=:black, label="Our rates")
h_112 = CairoMakie.hist!(ax11, bwell_fixed_mses, bins=5, color=(kim_col, 0.5), strokewidth=1, strokecolor=:black, label="Kim et al. rates")
CairoMakie.Legend(hist[1, 1], [h_111, h_112], ["Our rates", "Kim et al. rates"], "Scheme 1", 
    titleposition=:left, framevisible=false, labelsize=8, nbanks=1)


ax12 = CairoMakie.Axis(hist[2,2])
h_121 = CairoMakie.hist!(ax12, pepke_m2_mses,       bins=15, color=(our_col, 0.5), strokewidth=1, strokecolor=:black, label="Our rates")
h_122 = CairoMakie.hist!(ax12, pepke_m2_fixed_mses, bins=6, color=(pepke_col, 0.5), strokewidth=1, strokecolor=:black, label="Pepke et al. rates")
CairoMakie.Legend(hist[1, 2], [h_121, h_122], ["Our rates", "Pepke et al. rates"], "Scheme 2", 
    titleposition=:left, framevisible=false, labelsize=8, nbanks=1)


ax21 = CairoMakie.Axis(hist[4,1])
h_211 = CairoMakie.hist!(ax21, faas_mses,           bins=15, color=(our_col, 0.5),   strokewidth=1, strokecolor=:black, label="Our rates")
h_212 = CairoMakie.hist!(ax21, faas_fixed_mses,     bins=3, color=(faas_col, 0.5),  strokewidth=1, strokecolor=:black, label="Faas et al. rates")
h_213 = CairoMakie.hist!(ax21, pepke_m1_fixed_mses, bins=4, color=(pepke_col, 0.5), strokewidth=1, strokecolor=:black, label="Pepke et al. rates")
CairoMakie.Legend(hist[3, 1], [h_211, h_212, h_213], ["Our rates", "Faas et al. rates", "Pepke et al. rates"], "Scheme 3", 
    titleposition=:left, framevisible=false, labelsize=8, nbanks=2)


ax22 = CairoMakie.Axis(hist[4,2])
h_221 = CairoMakie.hist!(ax22, byrne_mses,       bins=15, color=(our_col, 0.5),    strokewidth=1, strokecolor=:black, label="Our rates")
h_222 = CairoMakie.hist!(ax22, byrne_fixed_mses, bins=5, color=(byrne_col, 0.5),  strokewidth=1, strokecolor=:black, label="Byrne et al. rates")
CairoMakie.Legend(hist[3, 2], [h_221, h_222], ["Our rates", "Byrne et al. rates"], "Scheme 4", 
    titleposition=:left, framevisible=false, labelsize=8, nbanks=1)

CairoMakie.Label(hist[:,0], "Count", rotation=pi/2)
CairoMakie.Label(hist[5,:], "RMSE")

CairoMakie.rowsize!(hist.layout, 1, Relative(0.125))
CairoMakie.rowsize!(hist.layout, 3, Relative(0.125))
CairoMakie.rowsize!(hist.layout, 5, Relative(0.00))

CairoMakie.colsize!(hist.layout, 1, Relative(0.5))
CairoMakie.colsize!(hist.layout, 2, Relative(0.5))

CairoMakie.rowgap!(hist.layout, 4, 5)
CairoMakie.colgap!(hist.layout, 1, 10)


CairoMakie.save("Plots/hist.svg", hist, pt_per_unit=1, 
    px_per_unit=hist_px_per_unit, size=hist_size_units)
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
bp_colors = [our_col; kim_col; our_col; pepke_col; our_col; faas_col; pepke_col; our_col; byrne_col];
bp_rowlist = [2;3;5;6;8;9;10;12;13]
bp = bigplot([fpms_bwell[1];     fpms_bwell_fixed[1];
         fpms_pepke_m2[1];  fpms_pepke_m2_fixed[2];
         fpms_faas[1];  fpms_faas_fixed[1]; fpms_pepke_m1_fixed[1];
         fpms_byrne[1];  fpms_byrne_fixed[1]],
         val_pops_bigplot,
         bp_colors,
         bp_rowlist,
         names = ["Our rates"; "Kim et al.";
                  "Our rates"; "Pepke et al.";
                  "Our rates"; "Faas et al."; "Pepke et al.";
                  "Our rates"; "Byrne et al."],
        figsize=bp_size_units)

CairoMakie.Legend(bp[1, :], [bp.content[1].scene.plots[6],  bp.content[7].scene.plots[end]], 
    ["Our rates", "Kim et al. rates"], "Scheme 1", titleposition=:left, framevisible=false, labelsize=8, nbanks=3, 
    colors=[our_col, kim_col])
CairoMakie.Legend(bp[4, :], [bp.content[13].scene.plots[end], bp.content[19].scene.plots[end]], 
    ["Our rates", "Pepke et al. rates"], "Scheme 2", titleposition=:left, framevisible=false, labelsize=8, nbanks=3, 
    colors=[our_col, pepke_col])
CairoMakie.Legend(bp[7, :], [bp.content[25].scene.plots[end], bp.content[31].scene.plots[end], bp.content[37].scene.plots[end]], 
    ["Our rates", "Faas et al. rates", "Pepke et al."], "Scheme 3", titleposition=:left, framevisible=false, labelsize=8, nbanks=3, 
    colors=[our_col, faas_col, pepke_col])
CairoMakie.Legend(bp[11, :], [bp.content[43].scene.plots[end], bp.content[49].scene.plots[end]], 
    ["Our rates", "Byrne et al. rates"], "Scheme 4", titleposition=:left, framevisible=false, labelsize=8, nbanks=3, 
    colors=[our_col, byrne_col])

for i in 1:6
    colsize!(bp.layout, i, Relative(1/6))
end

for i in [2, 3, 5, 6, 8, 9, 10, 12, 13]
    CairoMakie.rowsize!(bp.layout, i, Relative(1/9))
end

CairoMakie.rowgap!(bp.layout, 10)
CairoMakie.rowgap!(bp.layout, 1, 15)
CairoMakie.rowgap!(bp.layout, 14, 5)

CairoMakie.colgap!(bp.layout, 0)
CairoMakie.colgap!(bp.layout, 1, 0)

CairoMakie.save("Plots/bigplot.svg", bp, pt_per_unit=1, 
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



bwell_CaM_tup = get_multiple_eqs(fpms_bwell, pract_pop, get_bwell_eqs);
bwell_fixed_CaM_tup = get_multiple_eqs(fpms_bwell_fixed, pract_pop, get_bwell_eqs);


eq_size_cm     = (17, 22.23)
eq_size_px     =  eq_size_cm .* dpc 
eq_size_units  =  eq_size_cm .* pt_per_cm
eq_px_per_unit = (eq_size_cm .* dpc ./ (eq_size_cm .* pt_per_cm))[1]

lims = (nothing, nothing, nothing, nothing)
f_eqs, f_eqs_ax11 = plot_CaM_eqs(ca_range, bwell_CaM_tup;       i=1, j=1, title="Scheme 1", 
    color=(our_col, 0.3),  f=nothing, xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates", figsize=eq_size_units);
_, _      = plot_CaM_eqs(ca_range, bwell_fixed_CaM_tup; i=1, j=1, title="Scheme 1", 
    color=(kim_col, 0.3), f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Kim et al. rates");
CairoMakie.scatter!(f_eqs_ax11, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.5));
CairoMakie.Legend(f_eqs[1, 2], f_eqs_ax11, "Scheme 1", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax11, -0.3, 4.3) 


pepke_m2_CaM_tup        = get_multiple_eqs(fpms_pepke_m2, pract_pop, get_pepke_m2_eqs);
pepke_m2_fixed_CaM_tup  = get_multiple_eqs(fpms_pepke_m2_fixed, pract_pop, get_pepke_m2_eqs);

f_eqs, f_eqs_ax12 = plot_CaM_eqs(ca_range, pepke_m2_CaM_tup;       i=2, j=1, title="Scheme 2", 
    color=(our_col, 0.3),  f=f_eqs, new_axis=true,   xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, pepke_m2_fixed_CaM_tup;         i=2, j=1, title="Scheme 2", 
    color=(pepke_col, 0.3), f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Pepke et al. rates");
CairoMakie.scatter!(f_eqs_ax12, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.5));
CairoMakie.Legend(f_eqs[2, 2], f_eqs_ax12, "Scheme 2", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax12, -0.3, 4.3) 




faas_CaM_tup        = get_multiple_eqs(fpms_faas, pract_pop, get_faas_eqs);
faas_fixed_CaM_tup  = get_multiple_eqs(fpms_faas_fixed, pract_pop, get_faas_eqs);
pepke_m1_fixed_CaM_tup  = get_multiple_eqs(fpms_pepke_m1_fixed, pract_pop, get_faas_eqs);

f_eqs, f_eqs_ax21 = plot_CaM_eqs(ca_range, faas_CaM_tup;       i=3, j=1, title="Scheme 3", 
    color=(our_col, 0.3),  f=f_eqs, new_axis=true, xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, faas_fixed_CaM_tup;         i=3, j=1, title="Scheme 3", 
    color=(faas_col, 0.3), f=f_eqs,                xlabel=nothing, ylabel=nothing, limits=lims, label="Faas et al. rates");
_, _      = plot_CaM_eqs(ca_range, pepke_m1_fixed_CaM_tup;         i=2, j=1, title="Scheme 3", 
    color=(pepke_col, 0.3), f=f_eqs,               xlabel=nothing, ylabel=nothing, limits=lims, label="Pepke et al. rates");
CairoMakie.scatter!(f_eqs_ax21, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.5));
CairoMakie.Legend(f_eqs[3, 2], f_eqs_ax21, "Scheme 3", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax21, -0.3, 4.3) 



byrne_CaM_tup        = get_multiple_eqs(fpms_byrne, pract_pop, get_byrne_eqs);
byrne_fixed_CaM_tup  = get_multiple_eqs(fpms_byrne_fixed, pract_pop, get_byrne_eqs);

f_eqs, f_eqs_ax22 = plot_CaM_eqs(ca_range, byrne_CaM_tup;       i=4, j=1, title="Scheme 4", 
    color=(our_col, 0.3),  f=f_eqs, new_axis=true, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, faas_fixed_CaM_tup;          i=4, j=1, title="Scheme 4", 
    color=(byrne_col, 0.3), f=f_eqs,                ylabel=nothing, limits=lims, label="Byrne et al. rates");
CairoMakie.scatter!(f_eqs_ax22, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. data", color=(:black, 0.5));
CairoMakie.Legend(f_eqs[4, 2], f_eqs_ax22, "Scheme 4", unique=true, valign=:center, framevisible=false, labelsize=8)
CairoMakie.ylims!(f_eqs_ax22, -0.3, 4.3) 


###CairoMakie.Label(f_eqs[:,0], L"Bound $\textrm{Ca}^{2+}$ per CaM", rotation=pi/2)
CairoMakie.Label(f_eqs[:,0], "# Bound Ca²⁺ per CaM", rotation=pi/2)
CairoMakie.colgap!(f_eqs.layout, 1, 10)
CairoMakie.save("Plots/Shifman_joint.svg", f_eqs, pt_per_unit=1, px_per_unit = eq_px_per_unit, size=eq_size_units)

#############################
### PRACTICAL EQ PLOT END ###
#############################
################################
### PRACTICAL DYN PLOT START ###
################################
dyn_pop   = gen_dyn_data()
t = dyn_pop[1].time;

bwell_CaM_dyn_tups = get_multiple_dyns(fpms_bwell, dyn_pop, get_bwell_dyns);
bwell_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_bwell_fixed, dyn_pop, get_bwell_dyns);

pepke_m2_CaM_dyn_tups = get_multiple_dyns(fpms_pepke_m2, dyn_pop, get_bwell_CN_dyns);
pepke_m2_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_pepke_m2_fixed, dyn_pop, get_bwell_CN_dyns);

faas_CaM_dyn_tups = get_multiple_dyns(fpms_faas, dyn_pop, get_faas_dyns);
faas_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_faas_fixed, dyn_pop, get_faas_dyns);
pepke_m1_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_pepke_m1_fixed, dyn_pop, get_faas_dyns);

byrne_CaM_dyn_tups = get_multiple_dyns(fpms_byrne, dyn_pop, get_byrne_dyns);
byrne_fixed_CaM_dyn_tups = get_multiple_dyns(fpms_byrne_fixed, dyn_pop, get_byrne_dyns);


dyn_size_cm     = (17, 12)
dyn_size_px     =  dyn_size_cm .* dpc 
dyn_size_units  =  dyn_size_cm .* pt_per_cm
dyn_px_per_unit = (dyn_size_cm .* dpc ./ (dyn_size_cm .* pt_per_cm))[1]
###lims = (nothing, nothing, nothing, nothing)
lims = (-50, 1550, -1, 21)
lw = 1

f_dyn = plot_CaM_dyns(t, bwell_CaM_dyn_tups;
                j=1, title="Our rates", color=(our_col, 0.3),
                f = nothing,
                figsize=dyn_size_units, lw=lw, limits=lims);
f_dyn = plot_CaM_dyns(t, bwell_fixed_CaM_dyn_tups;
                j=2, title="Kim et al.", color=(kim_col, 0.3),
                f = f_dyn,
                lw=1, limits=lims,
                ylabs=false, yticks=false);
CairoMakie.Legend(f_dyn[1, 1:2], [f_dyn.content[1].scene.plots[end], f_dyn.content[4].scene.plots[end]], 
    ["Our rates", "Kim et al."], "Scheme 1", framevisible=false, labelsize=8, nbanks=1, 
    colors=[our_col, byrne_col])


f_dyn = plot_CaM_dyns(t, pepke_m2_CaM_dyn_tups;
                j=3, title="Our rates", color=(our_col, 0.3),
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);
f_dyn = plot_CaM_dyns(t, pepke_m2_fixed_CaM_dyn_tups;
                j=4, title="Pepke et al.", color=(pepke_col, 0.3),
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);
CairoMakie.Legend(f_dyn[1, 3:4], [f_dyn.content[8].scene.plots[end], f_dyn.content[11].scene.plots[end]], 
    ["Our rates", "Pepke et al."], "Scheme 2", framevisible=false, labelsize=8, nbanks=1, 
    colors=[our_col, byrne_col])


f_dyn = plot_CaM_dyns(t, faas_CaM_dyn_tups;
                j=5, title="Our rates", color=(our_col, 0.3),
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);
f_dyn = plot_CaM_dyns(t, faas_fixed_CaM_dyn_tups;
                j=6, title="Faas et al.", color=(faas_col, 0.3),
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);
f_dyn = plot_CaM_dyns(t, pepke_m1_fixed_CaM_dyn_tups;
                j=7, title="Pepke et al.", color=(pepke_col, 0.3),
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);
CairoMakie.Legend(f_dyn[1, 5:7], [f_dyn.content[15].scene.plots[end], f_dyn.content[18].scene.plots[end], f_dyn.content[21].scene.plots[end]], 
    ["Our rates", "Faas et al.", "Pepke et al."], "Scheme 3", framevisible=false, labelsize=8, nbanks=2, 
    colors=[our_col, byrne_col])


f_dyn = plot_CaM_dyns(t, byrne_CaM_dyn_tups;
                j=8, title="Our rates", color=(our_col, 0.3),
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);
f_dyn = plot_CaM_dyns(t, byrne_fixed_CaM_dyn_tups;
                j=9, title="Byrne et al.", color=(byrne_col, 0.3),
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);
CairoMakie.Legend(f_dyn[1, 8:9], [f_dyn.content[25].scene.plots[end], f_dyn.content[28].scene.plots[end]], 
    ["Our rates", "Byrne et al."], "Scheme 4", framevisible=false, labelsize=8, nbanks=1, 
    colors=[our_col, byrne_col])

###CairoMakie.Label(f_dyn[1,0], L"$\textrm{CaM}_0$", rotation=pi/2)
###CairoMakie.Label(f_dyn[2,0], L"$\textrm{CaM}_p$", rotation=pi/2)
###CairoMakie.Label(f_dyn[3,0], L"$\textrm{CaM}_f$", rotation=pi/2)

CairoMakie.Label(f_dyn[2,0], "0Ca²⁺-CaM (μM)", rotation=pi/2)
CairoMakie.Label(f_dyn[3,0], "1-3Ca²⁺-CaM (μM)", rotation=pi/2)
CairoMakie.Label(f_dyn[4,0], "4Ca²⁺-CaM (μM)", rotation=pi/2)

CairoMakie.Label(f_dyn[5,:], "Time (ms)")

for i in 2:4
    rowsize!(f_dyn.layout, i, Relative(1/3.5))
end
for i in 0:9
    if i == 0
        colsize!(f_dyn.layout, i, Relative(0))
    else
        colsize!(f_dyn.layout, i, Relative(1/9))
    end
end

CairoMakie.colgap!(f_dyn.layout, 5)

colgap!(f_dyn.layout, 5)
colgap!(f_dyn.layout, 1, 10)
rowgap!(f_dyn.layout, 10)
rowgap!(f_dyn.layout, 1, 20)
rowgap!(f_dyn.layout, 4, 5)


CairoMakie.save("Plots/dyn_joint.svg", f_dyn, pt_per_unit=1, px_per_unit = dyn_px_per_unit, size=dyn_size_units)


##############################
### PRACTICAL DYN PLOT END ###
##############################
#########################
### PARAM PLOTS START ###
#########################
par_bwell         = get_param_lists(fpms_bwell, (:kon_1, :koff_1, :kon_2, :koff_2));
par_bwell_fixed   = get_param_lists(fpms_bwell_fixed, (:kon_1, :koff_1, :kon_2, :koff_2));

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


### ===================================== ###
bw_pp_size_cm     = (17, 17)
bw_pp_size_px     =  bw_pp_size_cm .* dpc 
bw_pp_size_units  =  bw_pp_size_cm .* pt_per_cm
bw_pp_px_per_unit = (bw_pp_size_cm .* dpc ./ (bw_pp_size_cm .* pt_per_cm))[1]

bwell_lims = (10^-12, 10^14, 10^-7, 10^14)
bwell_pairplot = pair_plots(
    [par_bwell, par_bwell_fixed],
###    (; kon_1=L"log$_{10}(k_1)$", koff_1=L"log$_{10}(k_2)$", 
###       kon_2=L"log$_{10}(k_3)$", koff_2=L"log$_{10}(k_4)$"),
    (; kon_1="log₁₀(k₁)", koff_1="log₁₀(k₂)", 
       kon_2="log₁₀(k₃)", koff_2="log₁₀(k₄)"),
    [(our_col, 0.5), (kim_col, 1.0)],
    ["Our rates", "Kim et al."];
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
colgap!(bwell_pairplot.layout, 5)
rowgap!(bwell_pairplot.layout, 1, 5)
rowgap!(bwell_pairplot.layout, 2, -22.5)
rowgap!(bwell_pairplot.layout, 3, -22.5)

CairoMakie.save("Plots/bwell_pplot.svg", bwell_pairplot, pt_per_unit=1, px_per_unit=bw_pp_px_per_unit,
    size=bw_pp_size_units)

### ===================================== ###
ppm2_pp_size_cm     = (17, 22.23)
ppm2_pp_size_px     =  ppm2_pp_size_cm .* dpc 
ppm2_pp_size_units  =  ppm2_pp_size_cm .* pt_per_cm
ppm2_pp_px_per_unit = (ppm2_pp_size_cm .* dpc ./ (ppm2_pp_size_cm .* pt_per_cm))[1]

pepke_lims = (10^-5, 10^10, 10^-5, 10^10)
###pepke_lims = (nothing, nothing, nothing, nothing)
pepke_m2_pairplot = pair_plots(
    [par_pepke_m2, par_pepke_m2_fixed],
###    (; kon_2C=L"log$_{10}(k_1)$", koff_2C=L"log$_{10}(k_2)$", 
###       kon_2N=L"log$_{10}(k_3)$", koff_2N=L"log$_{10}(k_4)$"),
    (; kon_TC="log₁₀(k₊ᵗᶜ)", koff_TC="log₁₀(k₋ᵗᶜ)", 
       kon_RC="log₁₀(k₊ʳᶜ)", koff_RC="log₁₀(k₋ʳᶜ)", 
       kon_TN="log₁₀(k₊ᵗⁿ)", koff_TN="log₁₀(k₋ᵗⁿ)", 
       kon_RN="log₁₀(k₊ʳⁿ)", koff_RN="log₁₀(k₋ʳⁿ)"),
    [(our_col, 0.5), (pepke_col, 1.0)],
    ["Our rates", "Pepke et al."];
    lims = pepke_lims,
    fig_w = ppm2_pp_size_units[1],
    fig_h = ppm2_pp_size_units[2],
    yvis=true,
    ms = 6,
    rot_xlab = pi/2
)
CairoMakie.Legend(pepke_m2_pairplot[0, :], pepke_m2_pairplot.content[1], "Scheme 2",
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

CairoMakie.save("Plots/pepke_m2_pplot.svg", pepke_m2_pairplot, pt_per_unit=1, px_per_unit=ppm2_pp_px_per_unit,
    size=ppm2_pp_size_units)

### ===================================== ###

faas_pp_size_cm     = (17, 22.23)
faas_pp_size_px     =  faas_pp_size_cm .* dpc 
faas_pp_size_units  =  faas_pp_size_cm .* pt_per_cm
faas_pp_px_per_unit = (faas_pp_size_cm .* dpc ./ (faas_pp_size_cm .* pt_per_cm))[1]

faas_lims = (10^-13, 10^13, 10^-13, 10^13)
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
    ["Our rates", "Faas et al.", "Pepke et al."];
    lims = faas_lims,
    fig_w = faas_pp_size_units[1],
    fig_h = faas_pp_size_units[2],
    ms=5,
    rot_xlab = pi/2
)
CairoMakie.Legend(faas_pairplot[0, :], faas_pairplot.content[1], "Scheme 3",
    unique=true, valign=:center, framevisible=false, labelsize=8, titlesize=8,
    nbanks=3)
for i in 1:7
    rowsize!(faas_pairplot.layout, i, Relative(1/7.6))
    colsize!(faas_pairplot.layout, i, Relative(1/6.5))
end
colgap!(faas_pairplot.layout, 3)
rowgap!(faas_pairplot.layout, 1, 5)

for i in 2:7
    rowgap!(faas_pairplot.layout, i, -25)
end

CairoMakie.save("Plots/faas_pplot.svg", faas_pairplot, pt_per_unit=1, px_per_unit=faas_pp_px_per_unit,
    size=faas_pp_size_units)

### ===================================== ###
byrne_pp_size_cm     = (17, 22.23)
byrne_pp_size_px     =  byrne_pp_size_cm .* dpc 
byrne_pp_size_units  =  byrne_pp_size_cm .* pt_per_cm
byrne_pp_px_per_unit = (byrne_pp_size_cm .* dpc ./ (byrne_pp_size_cm .* pt_per_cm))[1]

byrne_lims = (10^-18, 10^18, 10^-15, 10^18)
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
    ["Our rates", "Byrne et al."];
    lims = byrne_lims,
    fig_w = byrne_pp_size_units[1],
    fig_h = byrne_pp_size_units[2],
    ms=3,
    rot_xlab = pi/2
)
CairoMakie.Legend(byrne_pairplot[0, :], byrne_pairplot.content[1], "Scheme 4",
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

CairoMakie.save("Plots/byrne_pplot.svg", byrne_pairplot, pt_per_unit=1, px_per_unit=byrne_pp_px_per_unit,
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
CairoMakie.save("Plots/data_plots.svg", p_data, pt_per_unit=1, px_per_unit=pexp_pp_px_per_unit,
    size=pexp_pp_size_units)




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