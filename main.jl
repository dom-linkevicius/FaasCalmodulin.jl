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
##usingSerialization



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

###A = full_pop[1:13];
###B = full_pop[14:24];
###C = full_pop[25:35];
###D = full_pop[36:51];
###E = full_pop[52:65];
###F = full_pop[66:80];
###G = full_pop[81:94];

#####train_pop = [A[1:7  ]; B[1:7  ]; C[1:7  ]; D[1:7  ]; E[1:7  ]; F[1:7  ]; G[1:7]];
#####val_pop   = [A[8:10 ]; B[8:9  ]; C[8:9  ]; D[8:12 ]; E[8:11 ]; F[8:11 ]; G[8:11]];
#####test_pop  = [A[11:13]; B[10:11]; C[10:11]; D[13:16]; E[12:14]; F[12:15]; G[12:14]];
###train_pop = [A[1:7  ]; B[1:7  ]; D[1:7  ]; E[1:7  ]; F[1:7  ]; G[1:7]];
###val_pop   = [A[8:10 ]; B[8:9  ]; D[8:12 ]; E[8:11 ]; F[8:11 ]; G[8:11]];
###test_pop  = [A[11:13]; B[10:11]; D[13:16]; E[12:14]; F[12:15]; G[12:14]];


n_models = 20
seeds = 1:n_models
sub_range = nothing ###1:10:230



fpms_bwell, p_bwell_train, p_bwell_val, val_idxs_bwell, test_idxs_bwell = train_plot("Blackwell", 
    full_pop, 
    n_models,
    seeds;
    sub_range = sub_range,
);
serialize("Data/fpms_bwell.jls",fpms_bwell)
serialize("Data/val_idxs_bwell.jls",val_idxs_bwell)
serialize("Data/test_idxs_bwell.jls",test_idxs_bwell)

fpms_bwell_fixed, p_bwell_train_fixed, p_bwell_val_fixed, val_idxs_bwell_fixed, test_idxs_bwell_fixed = train_plot("Blackwell_fixed", 
    full_pop, 
    n_models,
    seeds; 
    sub_range = sub_range
);
serialize("Data/fpms_bwell_fixed.jls",fpms_bwell_fixed)



fpms_pepke_m2, p_pepke_m2_train, p_pepke_m2_val, val_idxs_pepke_m2, test_idxs_pepke_m2 = train_plot("Pepke_m2", 
    full_pop, 
    n_models, 
    seeds; 
    sub_range = sub_range
);
serialize("Data/fpms_pepke_m2.jls",fpms_pepke_m2)

fpms_pepke_m2_fixed, p_pepke_m2_train_fixed, p_pepke_m2_val_fixed, val_idxs_pepke_m2_fixed, test_idxs_pepke_m2_fixed = train_plot("Pepke_m2_fixed", 
    full_pop, 
    n_models, 
    seeds; 
    sub_range = sub_range
);
serialize("Data/fpms_pepke_m2_fixed.jls",fpms_pepke_m2_fixed)



fpms_faas, p_faas_train, p_faas_val, val_idxs_faas, test_idxs_faas = train_plot("Faas", 
    full_pop, 
    n_models,
    seeds;
    sub_range = sub_range
);
serialize("Data/fpms_faas.jls",fpms_faas)

fpms_faas_fixed, p_faas_train_fixed, p_faas_val_fixed, val_idxs_faas_fixed, test_idxs_faas_fixed = train_plot("Faas_fixed", 
    full_pop, 
    n_models,
    seeds; 
    sub_range = sub_range
);
serialize("Data/fpms_faas_fixed.jls",fpms_faas_fixed)

fpms_pepke_m1_fixed, p_pepke_m1_train_fixed, p_pepke_m1_val_fixed, val_idxs_pepke_m1_fixed, test_idxs_pepke_m1_fixed = train_plot("Pepke_m1_fixed", 
    full_pop, 
    n_models,
    seeds; 
    sub_range = sub_range
);
serialize("Data/fpms_pepke_m1_fixed.jls",fpms_pepke_m1_fixed)



fpms_byrne, p_byrne_train, p_byrne_val, val_idxs_byrne, test_idxs_byrne = train_plot("Byrne", 
    full_pop,
    n_models,
    seeds; 
    sub_range = sub_range,
);
serialize("Data/fpms_byrne.jls",fpms_byrne)

fpms_byrne_fixed, p_byrne_train_fixed, p_byrne_val_fixed, val_idxs_byrne_fixed, test_idxs_byrne_fixed = train_plot("Byrne_fixed", 
    full_pop, 
    n_models,
    seeds; 
    sub_range = sub_range
);
serialize("Data/fpms_byrne_fixed.jls",fpms_byrne_fixed)



###fpm_ms = msfit(
###    Blackwell_scheme,
###    fpms_bwell[1].data,
###    init_params(Blackwell_scheme),
###    JointMAP(),
###    3,
###    1, #or standard deviation of shooting penalty hyper prior if below argument is not nothing
###    :L1; #:L1 or :L2, :L1 should be highly favoured because it performed much better in experiments
###    hyper_prior_type = nothing, #nothing/:L1/:L2,
###    return_with_shooting = false, #true/false, if false, does a 0 iteration fit and returns a FittePumasModel with the original model
###    #regular Pumas.fit kwargs, e.g. optim_options, diffeq_options, constantcoef, etc.
###    optim_options=(; iterations=2000, store_trace=true),
###    diffeq_options=(;alg=Rodas5P(), abstol=1e-16)
###)


###########################
### AIC BOX PLOTS START ###
###########################

all_fpms = [fpms_bwell;    fpms_bwell_fixed;
            fpms_pepke_m2; fpms_pepke_m2_fixed;
            fpms_faas;     fpms_faas_fixed; fpms_pepke_m1_fixed;
            fpms_byrne;    fpms_byrne_fixed];


val_pops = [idx_to_val_vec(full_pop, idxs) for idxs in val_idxs_bwell];
all_val_pops = repeat(val_pops, 9)


all_names = [repeat(["Scheme 1 + our rates"], length(val_pops)); 
             repeat(["Scheme 1 + Kim et al. (2010)"], length(val_pops));
             repeat(["Scheme 2 + our rates"], length(val_pops));
             repeat(["Scheme 2 + Pepke et al. (2010)"], length(val_pops));
             repeat(["Scheme 3 + our rates"], length(val_pops));
             repeat(["Scheme 3 + Faas et al. (2011)"], length(val_pops));
             repeat(["Scheme 3 + Pepke et al. (2010)"], length(val_pops));
             repeat(["Scheme 4 + our rates"], length(val_pops));
             repeat(["Scheme 4 + Byrne et al. (2009)"], length(val_pops))];

aic_pairs = zip(all_fpms, all_val_pops, all_names);

x, y = aic_boxplots(aic_pairs);
x_makie = sort(repeat(1:length(unique(x)), Int(length(x)/length(unique(x)))))
y_makie = [i for i in y]

boxplot_w = 17*28.3465;
boxplot_h = 15*28.3465;
boxplot = CairoMakie.Figure(size=(boxplot_w, boxplot_h), fontsize=12)
ax_boxplot = CairoMakie.Axis(boxplot[1,1], 
    xticks=(unique(x_makie), unique(x)),
    xlabel="Scheme + parameters",
    ylabel="AIC",
    xticklabelrotation=pi/4,
    xticklabelpad=10,
    xlabelpadding=30,
    xticklabelsize=10,
    yticklabelsize=10
    )
CairoMakie.boxplot!(ax_boxplot, x_makie, y_makie, show_notch=false)
CairoMakie.save("Plots/AIC_boxplot.svg", boxplot, pt_per_unit=1, resolution=(boxplot_w, boxplot_h))



###Plots.box_makieplot(x, y, legend=false, ylabel="AIC", xlabel="", xrotation=35, ylim=[3e3, 2.25e4], 
###    bottom_margin=8mm, left_margin=6mm)
#########################
### AIC BOX PLOTS END ###
#########################

#######################
### VAL PLOTS START ###
#######################
val_vec = idx_to_val_vec(full_pop, val_idxs_byrne_fixed[1]);

im_bwell  = load("Plots/Simple_scheme.png");
p_m_std_bwell = plot_mean_std([fpms_bwell[1:1], fpms_bwell_fixed[1:1]],
           val_vec,
           ["red", "blue"]; im=im_bwell, alpha_val=0.2, legendfontsize=8, ylim=[0, 9])

im_bwell_CN  = load("Plots/Simple_scheme_CN.png");
p_m_std_bwell_CN = plot_mean_std([fpms_bwell_CN[1:1], fpms_bwell_fixed], 
    val_vec, 
    ["red", "blue"]; im=im_bwell_CN, alpha_val=0.2, legendfontsize=8)

im_faas  = load("Plots/Faas_scheme.png");
p_m_std_faas = plot_mean_std([fpms_faas, fpms_faas_fixed], 
    val_vec, 
    ["red", "blue"]; im=im_faas, alpha_val=0.2, legendfontsize=8)

im_byrne  = load("Plots/Byrne_scheme.png");
p_m_std_byrne = plot_mean_std([fpms_byrne, fpms_byrne_fixed], 
    val_vec, 
    ["red", "blue"]; im=im_byrne, alpha_val=0.2, legendfontsize=8)


val_pops_bigplot = idx_to_val_vec(full_pop, val_idxs_bwell[1]);
bp_w = 17*28.3465;
bp_h = 22*28.3465;
bp = bigplot([fpms_bwell[1];     fpms_bwell_fixed[1];
         fpms_pepke_m2[1];  fpms_pepke_m2_fixed[2];
         fpms_faas[1];  fpms_faas_fixed[1]; fpms_pepke_m1_fixed[1];
         fpms_byrne[1];  fpms_byrne_fixed[1]],
         val_pops_bigplot,
         names = ["Scheme 1 \n our rates"; "Scheme 1 \n Kim et al. (2010)";
                  "Scheme 2 \n our rates"; "Scheme 2 \n Pepke et al. (2010)";
                  "Scheme 3 \n our rates"; "Scheme 3 \n Faas et al. (2011) rates"; "Scheme 3 \n Pepke et al. (2010) rates";
                  "Scheme 4 \n our rates"; "Scheme 4 \n Byrne et al. (2009) rates"],
        figsize=(bp_w, bp_h))

CairoMakie.save("Plots/bigplot.svg", bp, pt_per_unit=1, resolution=(bp_w, bp_h))

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

ca_range    = 10 .^ LinRange(-9, 0, 100);
pract_pop   = gen_dummy_data(ca_range);



bwell_CaM_tup = get_multiple_eqs(fpms_bwell, pract_pop, get_bwell_eqs);
bwell_fixed_CaM_tup = get_multiple_eqs(fpms_bwell_fixed, pract_pop, get_bwell_eqs);

lims = (0, 7e-5, -0.3, 4.3);
eq_w = 17*28.3465;
eq_h = 22*28.3465;
f_eqs, f_eqs_ax11 = plot_CaM_eqs(ca_range, bwell_CaM_tup;       i=1, j=1, title="Scheme 1", 
    color=(:red, 0.5),  f=nothing, xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates", figsize=(eq_h, eq_h));
_, _      = plot_CaM_eqs(ca_range, bwell_fixed_CaM_tup; i=1, j=1, title="Scheme 1", 
    color=(:blue, 0.5), f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Kim et al.(2010) rates");
CairoMakie.scatter!(f_eqs_ax11, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. (2006) data", color=(:black, 0.5));
CairoMakie.Legend(f_eqs[1, 2], f_eqs_ax11, "Scheme 1", unique=true, valign=:center, framevisible=false, labelsize=8)
  


pepke_m2_CaM_tup        = get_multiple_eqs(fpms_pepke_m2, pract_pop, get_pepke_m2_eqs);
pepke_m2_fixed_CaM_tup  = get_multiple_eqs(fpms_pepke_m2_fixed, pract_pop, get_pepke_m2_eqs);

f_eqs, f_eqs_ax12 = plot_CaM_eqs(ca_range, pepke_m2_CaM_tup;       i=2, j=1, title="Scheme 2", 
    color=(:red, 0.5),  f=f_eqs, new_axis=true,   xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, pepke_m2_fixed_CaM_tup;         i=2, j=1, title="Scheme 2", 
    color=(:blue, 0.5), f=f_eqs,   xlabel=nothing, ylabel=nothing, limits=lims, label="Pepke et al.(2010) rates");
CairoMakie.scatter!(f_eqs_ax12, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. (2006) data", color=(:black, 0.5));
CairoMakie.Legend(f_eqs[2, 2], f_eqs_ax11, "Scheme 2", unique=true, valign=:center, framevisible=false, labelsize=8)




faas_CaM_tup        = get_multiple_eqs(fpms_faas, pract_pop, get_faas_eqs);
faas_fixed_CaM_tup  = get_multiple_eqs(fpms_faas_fixed, pract_pop, get_faas_eqs);
pepke_m1_fixed_CaM_tup  = get_multiple_eqs(fpms_pepke_m1_fixed, pract_pop, get_faas_eqs);

f_eqs, f_eqs_ax21 = plot_CaM_eqs(ca_range, faas_CaM_tup;       i=3, j=1, title="Scheme 3", 
    color=(:red, 0.5),  f=f_eqs, new_axis=true, xlabel=nothing, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, faas_fixed_CaM_tup;         i=3, j=1, title="Scheme 3", 
    color=(:blue, 0.5), f=f_eqs,                xlabel=nothing, ylabel=nothing, limits=lims, label="Faas et al.(2011) rates");
_, _      = plot_CaM_eqs(ca_range, pepke_m1_fixed_CaM_tup;         i=2, j=1, title="Scheme 3", 
    color=(:green, 0.5), f=f_eqs,               xlabel=nothing, ylabel=nothing, limits=lims, label="Pepke et al.(2010) rates");
CairoMakie.scatter!(f_eqs_ax21, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. (2006) data", color=(:black, 0.5));
CairoMakie.Legend(f_eqs[3, 2], f_eqs_ax11, "Scheme 3", unique=true, valign=:center, framevisible=false, labelsize=8)



byrne_CaM_tup        = get_multiple_eqs(fpms_byrne, pract_pop, get_byrne_eqs);
byrne_fixed_CaM_tup  = get_multiple_eqs(fpms_byrne_fixed, pract_pop, get_byrne_eqs);

f_eqs, f_eqs_ax22 = plot_CaM_eqs(ca_range, byrne_CaM_tup;       i=4, j=1, title="Scheme 4", 
    color=(:red, 0.5),  f=f_eqs, new_axis=true, ylabel=nothing, limits=lims, label="Our rates");
_, _      = plot_CaM_eqs(ca_range, faas_fixed_CaM_tup;          i=4, j=1, title="Scheme 4", 
    color=(:blue, 0.5), f=f_eqs,                ylabel=nothing, limits=lims, label="Byrne et al.(2009) rates");
CairoMakie.scatter!(f_eqs_ax22, X_shifman, Y_shifman, 
    marker=:utriangle, label="Shifman et al. (2006) data", color=(:black, 0.5));
CairoMakie.Legend(f_eqs[4, 2], f_eqs_ax22, "Scheme 4", unique=true, valign=:center, framevisible=false, labelsize=8)

CairoMakie.Label(f_eqs[:,0], L"Bound $\textrm{Ca}^{2+}$ per CaM", rotation=pi/2)
CairoMakie.save("Plots/Shifman_joint2.svg", f_eqs, pt_per_unit=1, resolution=(eq_w, eq_h))

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


col = (:red, 0.5)

dyn_w = 17*28.3465;
dyn_h = 22*28.3465;

###lims = (nothing, nothing, nothing, nothing)
lims = (-50, 1550, -1, 21)
lw = 1

f_dyn = plot_CaM_dyns(t, bwell_CaM_dyn_tups;
                j=1, title="Scheme 1 \n Our rates", color=col,
                f = nothing,
                figsize=(dyn_w, dyn_h), lw=lw, limits=lims);

f_dyn = plot_CaM_dyns(t, bwell_fixed_CaM_dyn_tups;
                j=2, title="Scheme 1 \n Kim et al.", color=col,
                f = f_dyn,
                lw=1, limits=lims,
                ylabs=false, yticks=false);

f_dyn = plot_CaM_dyns(t, pepke_m2_CaM_dyn_tups;
                j=3, title="Scheme 2 \n Our rates", color=col,
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);

f_dyn = plot_CaM_dyns(t, pepke_m2_fixed_CaM_dyn_tups;
                j=4, title="Scheme 2 \n Pepke et al.", color=col,
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);

f_dyn = plot_CaM_dyns(t, faas_CaM_dyn_tups;
                j=5, title="Scheme 3 \n Our rates", color=col,
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);

f_dyn = plot_CaM_dyns(t, faas_fixed_CaM_dyn_tups;
                j=6, title="Scheme 3 \n Faas et al.", color=col,
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);

f_dyn = plot_CaM_dyns(t, pepke_m1_fixed_CaM_dyn_tups;
                j=7, title="Scheme 3 \n Pepke et al.", color=col,
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);

f_dyn = plot_CaM_dyns(t, byrne_CaM_dyn_tups;
                j=8, title="Scheme 4 \n Our rates", color=col,
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);

f_dyn = plot_CaM_dyns(t, byrne_fixed_CaM_dyn_tups;
                j=9, title="Scheme 4 \n Byrne et al.", color=col,
                f = f_dyn,
                lw=lw, limits=lims,
                ylabs=false, yticks=false);

CairoMakie.Label(f_dyn[1,0], L"$\textrm{CaM}_0$", rotation=pi/2)
CairoMakie.Label(f_dyn[2,0], L"$\textrm{CaM}_p$", rotation=pi/2)
CairoMakie.Label(f_dyn[3,0], L"$\textrm{CaM}_f$", rotation=pi/2)

CairoMakie.Label(f_dyn[4,:], "Time (ms)")

for i in 1:3
    rowsize!(f_dyn.layout, i, Relative(1/3))
end
for i in 0:9
    if i == 0
        colsize!(f_dyn.layout, i, Relative(0))
    else
        colsize!(f_dyn.layout, i, Relative(1/9))
    end
end

CairoMakie.save("Plots/dyn_joint.svg", f_dyn, pt_per_unit=1, resolution=(dyn_w, dyn_h))


##############################
### PRACTICAL DYN PLOT END ###
##############################
#########################
### PARAM PLOTS START ###
#########################


par_bwell         = get_param_lists(fpms_bwell, (:kon_1, :koff_1, :kon_2, :koff_2));
par_bwell_fixed   = get_param_lists(fpms_bwell_fixed, (:kon_1, :koff_1, :kon_2, :koff_2));

par_pepke_m2      = get_param_lists(fpms_pepke_m2, (:kon_2C, :koff_2C, :kon_2N, :koff_2N));
par_pepke_m2_fixed= get_param_lists(fpms_pepke_m2_fixed, (:kon_2C, :koff_2C, :kon_2N, :koff_2N));


par_faas          = get_param_lists(fpms_faas, (:kon_TN, :koff_TN, 
                                                  :kon_RN, :koff_RN,
                                                  :kon_TC, :koff_TC,
                                                  :kon_RC, :koff_RC));
par_faas_fixed    = get_param_lists(fpms_faas_fixed, (:kon_TN, :koff_TN, 
                                                  :kon_RN, :koff_RN,
                                                  :kon_TC, :koff_TC,
                                                  :kon_RC, :koff_RC));
par_pepke_m1_fixed= get_param_lists(fpms_pepke_m1_fixed, (:kon_TN, :koff_TN, 
                                                  :kon_RN, :koff_RN,
                                                  :kon_TC, :koff_TC,
                                                  :kon_RC, :koff_RC));

par_byrne         = get_param_lists(fpms_byrne, (:k01_N, :k10_N, 
                                                     :k02_N, :k20_N,
                                                     :k13_N, :k31_N,
                                                     :k23_N, :k32_N,
                                                     :k01_C, :k10_C, 
                                                     :k02_C, :k20_C,
                                                     :k13_C, :k31_C,
                                                     :k23_C, :k32_C));                                                 

par_byrne_fixed   = get_param_lists(fpms_byrne_fixed, (:k01_N, :k10_N, 
                                                     :k02_N, :k20_N,
                                                     :k13_N, :k31_N,
                                                     :k23_N, :k32_N,
                                                     :k01_C, :k10_C, 
                                                     :k02_C, :k20_C,
                                                     :k13_C, :k31_C,
                                                     :k23_C, :k32_C));


par_w = 17*28.3465;
par_h = 23*28.3465;
msize = 15

lims = (nothing, nothing, nothing, nothing)
c1 = (:red, 0.5)
c2 = (:blue, 0.5)
c3 = (:green, 0.5)

fp = CairoMakie.Figure(size=(par_w, par_h), fontsize=6)

bwell_lims = (1e3, 1e10, 1e-6, 1e2)
plot_params_ax(fp, 1, 1:4,
    [par_bwell[:kon_1], par_bwell_fixed[:kon_1]],
    [par_bwell[:koff_1], par_bwell_fixed[:koff_1]],
    L"$k_1$", 
    L"$k_2$",
    [c1, c2, c3],
    ["Our rates", "Kim et al. (2010)"];
    lims=bwell_lims
)
plot_params_ax(fp, 1, 5:8,
    [par_bwell[:kon_2], par_bwell_fixed[:kon_2]],
    [par_bwell[:koff_2], par_bwell_fixed[:koff_2]],
    L"$k_3$", 
    L"$k_4$",
    [c1, c2, c3],
    ["Our rates", "Kim et al. (2010)"];
    lims=bwell_lims,
    ylabs=false,
    yticks=false
)


pepke_m2_lims = (1e4, 1e11, 1e-3, 1e4)
plot_params_ax(fp, 2, 1:4,
    [par_pepke_m2[:kon_2C],  par_pepke_m2_fixed[:kon_2C]],
    [par_pepke_m2[:koff_2C], par_pepke_m2_fixed[:koff_2C]],
    L"$k_1$", 
    L"$k_2$",
    [c1, c2, c3],
    ["Our rates", "Pepke et al. (2010)"],
    lims=pepke_m2_lims,
)
plot_params_ax(fp, 2, 5:8,
    [par_pepke_m2[:kon_2N],  par_pepke_m2_fixed[:kon_2N]],
    [par_pepke_m2[:koff_2N], par_pepke_m2_fixed[:koff_2N]],
    L"$k_3$", 
    L"$k_4$",
    [c1, c2, c3],
    ["Our rates", "Pepke et al. (2010)"],
    lims=pepke_m2_lims,
    ylabs=false,
    yticks=false
)


faas_lims = (1e1, 1e9, 1e-4, 1e9)
plot_params_ax(fp, 3, 1:2,
    [par_faas[:kon_TC],  par_faas_fixed[:kon_TC], par_pepke_m1_fixed[:kon_TC]],
    [par_faas[:koff_TC],  par_faas_fixed[:koff_TC], par_pepke_m1_fixed[:koff_TC]],
    L"$k_1$", 
    L"$k_2$",
    [c1, c2, c3],
    ["Our rates", "Faas et al. (2011)", "Pepke et al. (2010)"],
    lims=faas_lims,
)
plot_params_ax(fp, 3, 3:4,
    [par_faas[:kon_RC],  par_faas_fixed[:kon_RC], par_pepke_m1_fixed[:kon_RC]],
    [par_faas[:koff_RC],  par_faas_fixed[:koff_RC], par_pepke_m1_fixed[:koff_RC]],
    L"$k_3$", 
    L"$k_4$",
    [c1, c2, c3],
    ["Our rates", "Faas et al. (2011)", "Pepke et al. (2010)"],
    lims=faas_lims,
    ylabs=false,
    yticks=false
)
plot_params_ax(fp, 3, 5:6,
    [par_faas[:kon_TN],  par_faas_fixed[:kon_TN], par_pepke_m1_fixed[:kon_TN]],
    [par_faas[:koff_TN],  par_faas_fixed[:koff_TN], par_pepke_m1_fixed[:koff_TN]],
    L"$k_5$", 
    L"$k_6$",
    [c1, c2, c3],
    ["Our rates", "Faas et al. (2011)", "Pepke et al. (2010)"],
    lims=faas_lims,
    ylabs=false,
    yticks=false
)
plot_params_ax(fp, 3, 7:8,
    [par_faas[:kon_TC],  par_faas_fixed[:kon_TC], par_pepke_m1_fixed[:kon_TC]],
    [par_faas[:koff_TC],  par_faas_fixed[:koff_TC], par_pepke_m1_fixed[:koff_TC]],
    L"$k_7$", 
    L"$k_8$",
    [c1, c2, c3],
    ["Our rates", "Faas et al. (2011)", "Pepke et al. (2010)"],
    lims=faas_lims,
    ylabs=false,
    yticks=false
)



byrne_lims = (1e1, 1e10, 1e-4, 1e9)
plot_params_ax(fp, 4, 1,
    [par_byrne[:k01_C], par_byrne_fixed[:k01_C]],
    [par_byrne[:k10_C], par_byrne_fixed[:k10_C]],
    L"$k^C_{01}$", 
    L"$k^C_{10}$",
    [c1, c2, c3],
    ["Our rates", "Byrne et al. (2009)"],
    lims=byrne_lims,
)
plot_params_ax(fp, 4, 2,
    [par_byrne[:k02_C], par_byrne_fixed[:k02_C]],
    [par_byrne[:k20_C], par_byrne_fixed[:k20_C]],
    L"$k^C_{02}$", 
    L"$k^C_{20}$",
    [c1, c2, c3],
    ["Our rates", "Byrne et al. (2009)"],
    lims=byrne_lims,
    ylabs=false,
    yticks=false
)
plot_params_ax(fp, 4, 3,
    [par_byrne[:k13_C], par_byrne_fixed[:k13_C]],
    [par_byrne[:k31_C], par_byrne_fixed[:k31_C]],
    L"$k^C_{13}$", 
    L"$k^C_{31}$",
    [c1, c2, c3],
    ["Our rates", "Byrne et al. (2009)"],
    lims=byrne_lims,
    ylabs=false,
    yticks=false
)
plot_params_ax(fp, 4, 4,
    [par_byrne[:k23_C], par_byrne_fixed[:k23_C]],
    [par_byrne[:k32_C], par_byrne_fixed[:k32_C]],
    L"$k^C_{23}$", 
    L"$k^C_{32}$",
    [c1, c2, c3],
    ["Our rates", "Byrne et al. (2009)"],
    lims=byrne_lims,
    ylabs=false,
    yticks=false
)
plot_params_ax(fp, 4, 5,
    [par_byrne[:k01_N], par_byrne_fixed[:k01_N]],
    [par_byrne[:k10_N], par_byrne_fixed[:k10_N]],
    L"$k^N_{01}$", 
    L"$k^N_{10}$",
    [c1, c2, c3],
    ["Our rates", "Byrne et al. (2009)"],
    lims=byrne_lims,
    ylabs=false,
    yticks=false
)
plot_params_ax(fp, 4, 6,
    [par_byrne[:k02_N], par_byrne_fixed[:k02_N]],
    [par_byrne[:k20_N], par_byrne_fixed[:k20_N]],
    L"$k^N_{02}$", 
    L"$k^N_{20}$",
    [c1, c2, c3],
    ["Our rates", "Byrne et al. (2009)"],
    lims=byrne_lims,
    ylabs=false,
    yticks=false
)
plot_params_ax(fp, 4, 7,
    [par_byrne[:k13_N], par_byrne_fixed[:k13_N]],
    [par_byrne[:k31_N], par_byrne_fixed[:k31_N]],
    L"$k^N_{13}$", 
    L"$k^N_{31}$",
    [c1, c2, c3],
    ["Our rates", "Byrne et al. (2009)"],
    lims=byrne_lims,
    ylabs=false,
    yticks=false
)
plot_params_ax(fp, 4, 8,
    [par_byrne[:k23_N], par_byrne_fixed[:k23_N]],
    [par_byrne[:k32_N], par_byrne_fixed[:k32_N]],
    L"$k^N_{23}$", 
    L"$k^N_{32}$",
    [c1, c2, c3],
    ["Our rates", "Byrne et al. (2009)"],
    lims=byrne_lims,
    ylabs=false,
    yticks=false
)

CairoMakie.Label(fp[1,:], "Scheme 1", valign=:top)
CairoMakie.Label(fp[2,:], "Scheme 2", valign=:top)
CairoMakie.Label(fp[3,:], "Scheme 3", valign=:top)
CairoMakie.Label(fp[4,:], "Scheme 4", valign=:top)

for i in 1:4
    rowsize!(fp.layout, i, Relative(0.25))
end
for i in 0:9
    if i == 0
        colsize!(f_dyn.layout, i, Relative(0))
    else
        colsize!(f_dyn.layout, i, Relative(1/9))
    end
end

CairoMakie.save("Plots/par_joint.svg", fp, pt_per_unit=1, resolution=(par_w, par_h))


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
pexp_w = 17*28.3465;
pexp_h = 7*28.3465;
p_data = plot_exp_data([A,B,G]; figsize=(pexp_w, pexp_h), lw=1, ylim=[0, 22.5])
CairoMakie.save("Plots/data_plots.svg", p_data, pt_per_unit=1, resolution=(pexp_w, pexp_h))




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