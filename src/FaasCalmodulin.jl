module FaasCalmodulin


using Serialization
using Random
using XLSX
using NumericalIntegration
using HypothesisTests

using Symbolics
using Pumas
using DeepPumas
import NNlib.sigmoid
using CairoMakie
using FileIO


include("symbolic_manipulations.jl")
include("process_data.jl")
export fetch_faas, fetch_shifman
include("models.jl")
include("model_training.jl")
export train_n_models, plot_by_condition, load_fpms, save_fpms, save_non_fpms_outputs
include("metrics.jl")
export rmse,
    shifman_error,
    my_aic,
    variance,
    cohen_d,
    print_mses_training,
    print_mses_validation,
    print_mses_testing
include("plot_recipes.jl")
export plot_by_condition, faas_data_plots, shifman_data_plot, violin_plot, big_plot
export equilibrium_plot, aic_plot, integration_plot, correlation_plot, byrne_c_n_plot
export blackwell_parameter_pairplot,
    bhalla_parameter_pairplot, shifman_parameter_pairplot, pepke_m2_parameter_pairplot
export faas_parameter_pairplot, byrne_parameter_pairplot

end
