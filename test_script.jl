using FaasCalmodulin
using Pumas
using DeepPumas
using CairoMakie

faas_data = fetch_faas("data/Faas_2011.xlsx");
shifman_data = fetch_shifman("data/Shifman_2006.xlsx");

###crouch_raw  = XLSX.readxlsx("data/Crouch_Klee_1980.xlsx");
###shifman_raw = XLSX.readxlsx("data/Shifman_2006.xlsx");
###
###f = Figure()
###ax = Axis(f[1,1], 
###    xlabel="Free Ca²⁺",
###    ylabel="Ca²⁺ per CaM",
###    limits=(nothing, 1e-4, nothing, nothing)
###    )
###scatter!(ax, 10 .^(crouch_raw["Crouch_Klee_1980"]["A"][:]), 
###    Float64.(crouch_raw["Crouch_Klee_1980"]["B"][:]), 
###    label="Crouch & Klee 1980 (25°C)", alpha=0.5)
###scatter!(ax, Float64.(shifman_raw["Shifman_2006.csv"]["A"][2:end]) * 1e-6, 
###    Float64.(shifman_raw["Shifman_2006.csv"]["B"][2:end]), 
###    label="Shifman et al. 2006 (35°C)", alpha=0.5)
###Legend(f[1,2], ax)

################################
### GLOBAL PLOTTING SETTINGS ###
################################

regular_font_string = "fonts/Arial Regular.ttf"
bold_font_string = "fonts/Arial Bold.ttf"
CairoMakie.set_theme!(fonts=(;regular = regular_font_string, bold= bold_font_string))

gs = (
    fontsize=8,
    dpc = 236.22, ### dots per cm, this is equal to 600dpi
    pt_per_cm = 28.3465,
    our_col     = CairoMakie.RGB( 34/255,113/255,178/255),
    bhalla_col  = CairoMakie.RGB(061/255,183/255,233/255),
    kim_col     = CairoMakie.RGB(247/255,072/255,165/255),
    pepke_col   = CairoMakie.RGB(053/255,155/255,115/255),
    faas_col    = CairoMakie.RGB(213/255,094/255,000/255),
    shifman_col = CairoMakie.RGB(230/255,159/255,000/255),
    byrne_col   = CairoMakie.RGB(240/255,228/255,066/255)
);

######################
### LOADING MODELS ###
###################### 

fpms_nt, val_idxs, test_idxs = load_fpms("trained_models/",
    (:blackwell,
    :blackwell_fixed,
    :bhalla,
    :bhalla_fixed,
    :shifman,
    :shifman_fixed,
    :pepke_m2,
    :pepke_m2_fixed,
    :faas,
    :faas_fixed,
    :pepke_m1_fixed,
    :byrne,
    :byrne_fixed,
    ),
    ".jls"
);



####################################
### TRAINING MODELS FROM SCRATCH ###
#################################### 

n_runs = 20;
seeds = 1:n_runs;
sub_idxs = 1:20:201;
n_iterations = 2000;

fpms_blackwell, _, _, val_idxs, test_idxs = train_n_models(
    :Blackwell, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=false, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_blackwell_fixed, _, _, _, _ = train_n_models(
    :Blackwell, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=true, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_bhalla, _, _, _, _ = train_n_models(
    :Bhalla, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=false, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_bhalla_fixed, _, _, _, _ = train_n_models(
    :Bhalla, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=true, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_shifman, _, _, _, _ = train_n_models(
    :Shifman, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=false, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_shifman_fixed, _, _, _, _ = train_n_models(
    :Shifman, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=true, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_pepke_m2, _, _, _, _ = train_n_models(
    :Pepke_m2, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=false, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_pepke_m2_fixed, _, _, _, _ = train_n_models(
    :Pepke_m2, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=true, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_faas, _, _, _, _ = train_n_models(
    :Faas, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=false, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_faas_fixed, _, _, _, _ = train_n_models(
    :Faas, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=true, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_pepke_m1_fixed, _, _, _, _ = train_n_models(
    :Pepke_m1, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=true, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_byrne, _, _, _, _ = train_n_models(
    :Byrne, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=false, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


fpms_byrne_fixed, _, _, _, _ = train_n_models(
    :Byrne, 
    faas_data, 
    shifman_data, 
    n_runs, 
    seeds, 
    sub_range = sub_idxs,
    fixed=true, 
    plot=false,
    alg=JointMAP(), 
    diffeq_options=(;alg=Rodas5P(), abstol=1e-16),
    optim_options=(;iterations=n_iterations, 
        store_trace=true, 
        allow_f_increases=true, 
        g_tol=1e-6
    )
);


save_fpms("trained_models/",
    (
    blackwell           = fpms_blackwell,
    blackwell_fixed     = fpms_blackwell_fixed,
    bhalla              = fpms_bhalla,
    bhalla_fixed        = fpms_bhalla_fixed,
    shifman             = fpms_shifman,
    shifman_fixed       = fpms_shifman_fixed,
    pepke_m2            = fpms_pepke_m2,
    pepke_m2_fixed      = fpms_pepke_m2_fixed,
    faas                = fpms_faas,
    faas_fixed          = fpms_faas_fixed,
    pepke_m1_fixed      = fpms_pepke_m1_fixed,
    byrne               = fpms_byrne,
    byrne_fixed         = fpms_byrne_fixed
    ),
    val_idxs,
    test_idxs,
    ".jls"
);

save_non_fpms_outputs("trained_params/",
###    (
###    blackwell           = fpms_blackwell,
###    blackwell_fixed     = fpms_blackwell_fixed,
###    bhalla              = fpms_bhalla,
###    bhalla_fixed        = fpms_bhalla_fixed,
###    shifman             = fpms_shifman,
###    shifman_fixed       = fpms_shifman_fixed,
###    pepke_m2            = fpms_pepke_m2,
###    pepke_m2_fixed      = fpms_pepke_m2_fixed,
###    faas                = fpms_faas,
###    faas_fixed          = fpms_faas_fixed,
###    pepke_m1_fixed      = fpms_pepke_m1_fixed,
###    byrne               = fpms_byrne,
###    byrne_fixed         = fpms_byrne_fixed
###    ),
    fpms_nt, 
    test_mses_nt,
    ".jls"
);


#########################
### PLOTS AND OUTPUTS ###
#########################

train_mses_nt = print_mses_training(fpms_nt);
valid_mses_nt = print_mses_validation(fpms_nt, faas_data, val_idxs);
test_mses_nt  = print_mses_testing(fpms_nt, faas_data, test_idxs);

faas_figure        = faas_data_plots(faas_data, gs; save=true) 
shifman_figure     = shifman_data_plot("data/Shifman_2006.xlsx", gs; save=true)
violin_figure      = violin_plot(test_mses_nt, gs; save=true)
big_figure         = big_plot(fpms_nt, faas_data, val_idxs, gs; save=true)
equilibrium_figure = equilibrium_plot(fpms_nt, gs; save=true)
aic_figure         = aic_plot(fpms_nt, faas_data, test_idxs, gs; save=true)
integration_figure = integration_plot(fpms_nt, gs; save=true)
correlation_figure = correlation_plot(fpms_nt, gs; save=true)
byrne_c_n_figure   = byrne_c_n_plot(fpms_nt[:byrne_fixed], gs; save=true)

############################
### SUPPLEMENTAL FIGURES ###
############################

blackwell_pairplot_figure = blackwell_parameter_pairplot(
    fpms_nt[:blackwell], fpms_nt[:blackwell_fixed], gs; save=true)
bhalla_pairplot_figure = bhalla_parameter_pairplot(
    fpms_nt[:bhalla], fpms_nt[:bhalla_fixed], gs; save=true)
shifman_pairplot_figure = shifman_parameter_pairplot(
    fpms_nt[:shifman], fpms_nt[:shifman_fixed], gs; save=true)
pepke_m2_pairplot_figure = pepke_m2_parameter_pairplot(
    fpms_nt[:pepke_m2], fpms_nt[:pepke_m2_fixed], gs; save=true)
faas_pairplot_figure = faas_parameter_pairplot(
    fpms_nt[:faas], fpms_nt[:faas_fixed], fpms_nt[:pepke_m1_fixed], gs; save=true)
byrne_pairplot_figure = byrne_parameter_pairplot(
    fpms_nt[:byrne], fpms_nt[:byrne_fixed], gs; save=true)


###################################
### FENS 2024 PLOTS AND OUTPUTS ###
###################################

faas_figure_fens    = data_plots_fens(faas_data, "data/Shifman_2006.xlsx", gs; save=true, w=30, h=13) 