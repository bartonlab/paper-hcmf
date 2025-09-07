using DelimitedFiles
using Distributed
using Random
using Distributions
using LinearAlgebra
using Printf
rng = Random.MersenneTwister(1234);

include("./systematic_performance_evaluation.jl");
include("./basic_inference.jl")

L = parse(Int, ARGS[1])
time_upto = parse(Int, ARGS[2])
t_interval = parse(Int, ARGS[3])
id_ensemble_begin = parse(Int, ARGS[4])
id_ensemble_max = parse(Int, ARGS[5])
flag_resample = parse(Bool, ARGS[6]) # true
n_resample = parse(Int, ARGS[7]) # 30 

dir_in = ARGS[8] # "./data_2024/WF/temp/"
file_key1 = ARGS[9] # "allele-traject-time-10"
file_key2 = ARGS[10] # "_N-1000_simple_wt_freq-selection.txt"
dir_out = ARGS[11]
file_key = ARGS[12]

mu = 10^-3
rec = 10^-4;
#density_coupling = 0.5*density_field
LLhalf = Int( 0.5 * L * (L-1) );
rank_x = L+LLhalf;
# ------------------------------------------------------------ #

for id_ensemble in id_ensemble_begin:id_ensemble_max
    # ------------------ process the data ------------------------ #
    (data, time_list) = get_data_time(file_key1, file_key2, dir_in, id_ensemble, time_upto);
    if(flag_resample) # Resampling
        data = resampling_with_N(data, n_resample)
        time_list = data[:, 1]
    end;
    time_unique = unique(time_list);

    index_to_look = data[:, 1] .%t_interval .== 0 
    data = copy(data[index_to_look, :])
    time_list = copy(time_list[time_list .%t_interval .== 0]) # The entries of time steps are only unique ones. 
    time_list = sort(unique(time_list)) ;
    # ------------------------------------------------------------ #
    (iCov, d_mu_tot, d_rec_tot, Δx_tot) = get_iCov_drift(rank_x, L, time_list, data);

    fout_iCov = open( @sprintf("%siCov_%s_%d.txt", dir_out, file_key, id_ensemble) , "w")
    fout_drift = open( @sprintf("%sidrift_%s_%d.txt", dir_out, file_key, id_ensemble) , "w")
    for n in 1:size(iCov, 1)
        str_out = join([@sprintf("%.8e", x) for x in iCov[n,:]], " ") 
        println(fout_iCov, str_out)
    end
    for n in 1:size(d_mu_tot, 1)
        str_out = @sprintf("%.8e %.8e %.8e %.8e", Δx_tot[n] - mu*d_mu_tot[n] + rec*d_rec_tot[n], d_mu_tot[n], d_rec_tot[n], Δx_tot[n])
        println(fout_drift, str_out)
    end
    close(fout_iCov); close(fout_drift)
    # ------------------------ Output ot computational time ------------------------- #
end;
