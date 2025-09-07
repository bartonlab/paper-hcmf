
include("./systematic_performance_evaluation.jl");

function get_Δt_set(t_set) 
    Δt_set = zeros(length(t_set))
    Δt_set[1] = 0.5 * (t_set[2] - t_set[1])
    Δt_set[end] = 0.5 * (t_set[end] - t_set[end-1])
    for i in 2:(length(t_set)-1)
        Δt_set[i] = 0.5 * (t_set[i+1] - t_set[i-1]) 
    end
    return Δt_set
end

function get_iCov_drift(rank_x, L, time_list, data)
    n_time_list_max = length(time_list);
    Δt_set = get_Δt_set(time_list);
    # ----------- For MPL based inference ---------- #
    d_mu_tot, d_rec_tot, Δx_tot = zeros(rank_x), zeros(rank_x), zeros(rank_x)
    #Ξ = zeros(rank_x, size(data,1) + length(time_list)-1);
    iCov = zeros(rank_x, rank_x)
    #-------------
    n_count = 1
    for id_t in 1:n_time_list_max
        t = time_list[id_t]
        Δt = Δt_set[id_t]
        
        (n_t, sample_t) = get_sample_at_t(data, t)
        (x_1w2_set, x_1w2, d_mu, d_rec) = get_x_1w2_set_and_d_mu_random_compression(sample_t, n_t, L, rank_x);
        if(id_t==1) 
            Δx_tot -= x_1w2
        end
        if(id_t==n_time_list_max) 
            Δx_tot += x_1w2 
        end
    
        d_mu_tot += d_mu * Δt
        d_rec_tot += d_rec * Δt
        for n in 1:size(x_1w2_set,1)
            iCov += (Δt * x_1w2_set[n, :]) * x_1w2_set[n, :]'
        end
    end
    return (iCov, d_mu_tot, d_rec_tot, Δx_tot)
end;

function resampling_with_N(data, n_resample)
    time_list_unique = sort(unique(data[:, 1]));
    data_out = zeros(size(data));
    n_walker = 0;
    for i_t in 1:length(time_list_unique)
        t = time_list_unique[i_t]
        id_t = data[:, 1] .== t;
        n_sample_list = data[id_t, 2];
        data_t = copy(data[id_t, :]);
        prob_resample = n_sample_list / sum(n_sample_list);
        dist = Multinomial(n_resample, prob_resample)
        n_new_sample_list = rand(dist);
        data_out_temp = copy(data_t[n_new_sample_list .> 0, :])
        data_out_temp[:, 2] = copy(n_new_sample_list[n_new_sample_list .>0])
        n_added = size(data_out_temp, 1)
        data_out[(n_walker+1):(n_walker+n_added), :] = copy(data_out_temp)
        n_walker += n_added
    end
    return data_out = copy(data_out[1:n_walker, :]);
end;
