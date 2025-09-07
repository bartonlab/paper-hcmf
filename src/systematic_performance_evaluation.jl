get_time_list(data) = Int.(unique(data[:,1]));

function get_data_time(file_key1, file_key2, data_dir, id_ensemble, time_upto=300)
    fname_in = data_dir * file_key1 *"_id-"*string(id_ensemble) * file_key2
    data = readdlm(fname_in);
    read_upto = count(data[:,1] .<= time_upto)
    data = copy(data[1:read_upto, :])
    time_list = get_time_list(data);
    return (data, time_list)
end;

function get_sample_at_t(data, t_get)
    n_t = []
    sample_t = []
    for n in 1:size(data,1)
        if(Int(data[n, 1]) == t_get)
            #n_t = push!(n_t, data[n,2]) 
            #sample_t = push!(sample_t, data[n,4:end]) 
            if(length(sample_t)>0)
                sample_t = hcat(sample_t, data[n,4:end]) 
                n_t = vcat(n_t, data[n,2])
            end
            if(length(sample_t)==0)
                sample_t = copy(data[n,4:end])
                n_t = data[n,2]
            end
        end

        if(data[n, 1] > t_get)
            break
        end
    end
    return (n_t, sample_t)
end;


function get_x1_x2(L, n_t, sample_t)
    Neff = Int(sum(n_t))
    x1 = zeros(L)
    x2 = zeros(L,L);
    n_species = size(sample_t,2);
    for m in 1:n_species
        x1 += n_t[m] / Neff * sample_t[1:L,m]
        x2 += n_t[m] / Neff * sample_t[1:L,m] * sample_t[1:L,m]'
    end

    return (x1, x2)
end;

function get_indices_i_and_j_matrix_to_ij_vector(mat_size)
    mat_temp = zeros(mat_size, mat_size)
    n_ij = 1
    index_upper_triangle = []
    for i in 1:mat_size
        for j in 1:mat_size
            mat_temp[i,j] = n_ij
            n_ij += 1
        end
    end;
    
    for i in 1:mat_size
        for j in (i+1):mat_size
            push!(index_upper_triangle, Int(mat_temp[j,i]))
        end
    end
    
    # mat_temp[index_upper_triangle] # gives upper triangle. 
    return index_upper_triangle
end;

function get_4th_Cov_mut(L, n_t, sample_t, index_upper_triangle)
    Neff = Int(sum(n_t))
    x1 = zeros(L)
    x2 = zeros(L,L);
    LLhalf = Int(L*(L-1)/2)
    x3 = zeros(L, LLhalf);
    x4 = zeros(LLhalf, LLhalf);
    n_species = size(sample_t,2);
    for m in 1:n_species
        n_t_scale = n_t[m] / Neff
        x1 += n_t_scale * sample_t[1:L,m]
        outer_prod_sample = copy(sample_t[1:L,m] * sample_t[1:L,m]')
        x2 += n_t_scale * outer_prod_sample
        x3 += n_t_scale * sample_t[1:L,m] * outer_prod_sample[index_upper_triangle]'
        x4 += n_t_scale * outer_prod_sample[index_upper_triangle] * outer_prod_sample[index_upper_triangle]' 
    end
    x_1w2 = [x1; x2[index_upper_triangle]]
    
    # The diagonal parts becom automatically beeng the correct ones. # eg. "diagonal of C_sinal" = x1 .* (1 .- x1)
    C_single = x2 - x1 * x1'
    C_triangle = x3 - x1 * x2[index_upper_triangle]'
    C_epistasis = x4 - x2[index_upper_triangle] * x2[index_upper_triangle]'
    
    d_mu = zeros(L+LLhalf)
    d_mu[1:L] = 1 .- 2 * x1;
    
    mat_temp = repeat(copy(x1), 1,L);
    mat_temp = mat_temp + mat_temp'
    d_mu[(L+1):end] = mat_temp[index_upper_triangle] - 4 * x2[index_upper_triangle]
    return (C_single, C_triangle, C_epistasis, d_mu, x_1w2);
end;

function flat_J_to_mat_J(s_epis, L)
    J=zeros(L,L)
    n = 1    
    for i in 1:L
        for j in (i+1):L
            J[i,j] = s_epis[L+n]
            J[j,i] = s_epis[L+n]
            n += 1
        end
    end
    return J
end;

function get_x_1w2(seq_in, L, LLhalf)
    x1 = seq_in[1:L];
    x2_vec = zeros(LLhalf);
    for i in 1:L
        for j in (i+1):L
            ξ = G(i,j,L)
            x2_vec[ξ] = seq_in[i]*seq_in[j]
        end
    end
    return ([x1;x2_vec])
end;

function get_4th_Cov_mut_simple(L, n_t, sample_t)
    Neff = Int(sum(n_t))
    LLhalf = Int(L*(L-1)/2)
    n_species = size(sample_t,2);

    x_1w2 = zeros(L+LLhalf) 
    x_3w4 = zeros(L+LLhalf, L+LLhalf)
    for m in 1:n_species
        n_t_scale = n_t[m] / Neff

        #should use L and LLhalf = L*(L-1)/2
        x_1w2_m = get_x_1w2(sample_t[1:L, m], L, LLhalf)
        x_1w2 += n_t_scale * x_1w2_m
        x_3w4 += n_t_scale * x_1w2_m * x_1w2_m'
    end
    C_tot = x_3w4 - x_1w2 * x_1w2'
    d_mu = get_d_mu(x_1w2, L, LLhalf);
    d_rec = get_d_rec(x_1w2, L, LLhalf);
    return (C_tot, x_1w2, d_mu, d_rec)
end;

function re_shape_additive_selection_epistasis(s_A, s_E, L, LLhalf)
    # s_A:  additive selection
    # s_E:  epistasis selection
    s_summarized = zeros(L+LLhalf)
    s_summarized[1:L] = copy(s_A)
    for i in 1:L
        for j in (i+1):L
            ξ = G(i,j,L) + L
            s_summarized[ξ] = s_E[i,j]
        end
    end
    return s_summarized
end;

function get_x_1w2_set_and_d_mu_random_compression(sample_t, n_t, L, rank_x)
    Neff = Int(sum(n_t))
    LLhalf = Int(L*(L-1)/2)
    n_species = size(sample_t,2);
    x_1w2 = zeros(L+LLhalf) 
    
    for m in 1:n_species
        n_t_scale = n_t[m] / Neff
        x_1w2_m = get_x_1w2(sample_t[1:L, m], L, LLhalf) 
        x_1w2 += n_t_scale * x_1w2_m
    end
    d_mu = get_d_mu(x_1w2, L, LLhalf);
    d_rec = get_d_rec(x_1w2, L, LLhalf);
    # make a table keeping x_1w2
    x_1w2_set = zeros(n_species, rank_x)
    
    for m in 1:n_species
        n_t_scale = n_t[m] / Neff
        # centerize and normalize
        x_1w2_m = sqrt(n_t[m] / Neff) * ( get_x_1w2(sample_t[1:L, m], L, LLhalf) - x_1w2 )
        x_1w2_set[m, :] = copy(x_1w2_m)
    end

    return (x_1w2_set, x_1w2, d_mu, d_rec)
end;
    
    

function get_d_mu(x_1w2, L, LLhalf)
    d_mu_temp = zeros(L+LLhalf)
    for i in 1:L
        d_mu_temp[i] = 1 - 2 * x_1w2[i]
        for j in (i+1):L
            ξ = G(i,j,L) + L
            d_mu_temp[ξ] = x_1w2[i] + x_1w2[j] - 4 * x_1w2[ ξ]
        end
    end
    return d_mu_temp
end;

function get_d_rec(x_1w2, L, LLhalf)
    d_rec_temp = zeros(L+LLhalf) 
    for i in 1:L
        for j in (i+1):L
            ξ = G(i,j,L) + L
            d_rec_temp[ξ] = (j-i) *(x_1w2[ξ] - x_1w2[i]*x_1w2[j])
        end
    end
    return d_rec_temp
end;


# ---------------- For Analysis ------------------------#
stepfunc(x) = x>0 ? 1 : 0
G(i,j,L) = stepfunc(j-i) * Int((i-1)*(L-i/2.0) + j-i);
#G(i,j,L) = (j-1) * L + i;

linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y;

get_time_list(data) = Int.(unique(data[:,1]))

# ------ this function enhacne the genetic drift ----
#N_subset is a number of the populatin size which is smaller than original population size N.
function resumpling_data(data, N, N_subset)
	data_resumpled = zeros(size(data));
	n_shift = 0
	for t_sampling in unique(data[:, 1])
	    #----- extract a specific time dependent subsamples
	    data_temp = data[ collect(1:size(data,1))[data[:,1] .== t_sampling],  : ]
	    #----- extract the number and id of sequences
	    temp_index = []
	    for n in 1:size(data_temp,1)
		temp_index = vcat(temp_index, [n for _ in 1:Int(data_temp[n,2]) ])
	    end
	    #------ Selecting the subsamples
	    index_subsamples = temp_index[shuffle(collect(1:N))[1:N_subset]];
	    #------ records the sequence id in the t-dependent subset
	    index_subsamples_unique = unique(index_subsamples)
	    #------ records the number of the corresponding sequences
	    n_count_subsamples = [count(index_subsamples .== n) for n in index_subsamples_unique];

	    #------ slect subsamples from the original time depending subsamples.
	    data_temp2 = data_temp[index_subsamples_unique, :]
	    #------ filling the number of the sequences
	    for n in 1:size(n_count_subsamples,1)
		data_temp2[n, 2] = n_count_subsamples[n]
	    end

	    #------ Insert resampled sequences to the output matrx.
	    data_resumpled[(n_shift+1):(n_shift+size(data_temp2,1)), :] = copy(data_temp2)
	    n_shift += size(data_temp2,1)
	end

	return data_resumpled = copy(data_resumpled[1:n_shift, :]);
end;


function get_epis_scatter_plot(epis_true, epis_est, ylabel_in; xlabel_in="True")
    th_top_10_percent = sort(vec(abs.(epis_est)), rev=true)[1000]
    index_to_look = abs.(epis_est) .> th_top_10_percent;
    x_axis = vec(epis_true[index_to_look])
    y_axis = vec(epis_est[index_to_look])
    cor_val_str = @sprintf("ρ=%.2f", cor(x_axis, y_axis))
    x_axis = vec(epis_true[index_to_look])
    y_axis = vec(epis_est[index_to_look])
    val_max = maximum([maximum(x_axis), maximum(y_axis)])
    Plots.plot(collect(-val_max:(val_max/5):val_max), collect(-val_max:(val_max/5):val_max), c="gray", label=:false)
    p1 = scatter!(x_axis, y_axis,c="blue",
    xlabel=xlabel_in,
    ylabel=ylabel_in,
    label=cor_val_str,
    legendfontsize=10,
    tickfontsize=10,
    markerstrokewidth=0,
    legend=(0.45,0.9),
    margin=4mm,
    markeralpha=0.4,
    xrotation=30,
    foreground_color_legend = nothing,
    size=(400,400))
    return p1
end;

function get_Ξ_masked(L, Ξ, my_num, NB)
    N_eff  = size(Ξ, 2)
    index_to_look = collect(1:L);
    #size_Ξ_eff = L + 2*NB*L;
    size_Ξ_eff = size(Ξ, 1)

    Ξ_masked = zeros(size_Ξ_eff, N_eff)
    my_num_masked = zeros(size_Ξ_eff)
    Ξ_masked[1:L, :] = copy(Ξ[1:L, :])
    my_num_masked[1:L] = copy(my_num[1:L])

    n_eff = L
    for i in 1:L
        for j in (i+1):L
            if(j <= i+NB)
                ξ = G(i,j,L) + L
                n_eff += 1
                Ξ_masked[n_eff, :] = copy(Ξ[ξ, :])
                my_num_masked[n_eff] = copy(my_num[ξ])
                push!(index_to_look, ξ)
            end
        end
    end
    return (index_to_look, Ξ_masked[1:n_eff, :], my_num_masked[1:n_eff])
end;

function get_S_J_Ground_Trueth(data_dir, file_key)
    # -------------------- Read True Parameters for Comparison -------------#
    @show data_dir * "log_" * file_key * ".txt"
    hyper_para = readdlm(data_dir * "log_" * file_key * ".txt");
    mu,rec, L, N, T_max = 0,0,0,0,0
    for n in 1:size(hyper_para,1)
        x = split(hyper_para[n, 1], ":")
        if(x[1]=="mu")  mu  = parse(Float64, x[2]) end
        if(x[1]=="rec") rec = parse(Float64, x[2]) end
        if(x[1]=="L" )  L   = parse(Int, x[2]) end
        if(x[1]=="N" )  N   = parse(Int, x[2]) end
        if(x[1]=="T" )  T_max  = parse(Int, x[2]) end
    end
    @printf("mu:%.1e, rec:%.1e, L:%.0e, N:%.0e, T:%.0e\n", mu, rec, L, N, T_max)

    local_field_raw = readdlm(data_dir * "local-field_" * file_key * ".txt");
    """
    S_GT = zeros(L) 
    for n in 1:size(local_field_raw,1)
        x = local_field_raw[n, :]
        i,v = Int(x[1]), x[2]
        S_GT[i] = v
    end;
    """
    S_GT = copy(local_field_raw)

    coupling_raw = readdlm(data_dir * "coupling_" * file_key * ".txt");
    """
    J_GT = zeros(L,L)
    for n in 1:size(coupling_raw,1)
        x = coupling_raw[n, :]
        i,j,v = Int(x[1]), Int(x[2]), x[3]
        J_GT[i,j] = v
        J_GT[j,i] = v
    end;
    """
    J_GT = copy(coupling_raw)
    return (mu, rec, L, N, T_max, S_GT, J_GT)
end;

function get_MinMax_additive_epistasis_GT(σ, h2)
    s = σ/5 # suppose the varinace of each Gaussiand distribution
    θ = sqrt(3/2) * (σ - s)
    J_GT_min, J_GT_max = +s-θ, -s+θ

    h = sqrt(h2)
    σ_A = h*σ / sqrt(1-h^2)
    s_A = σ_A/5
    θ_A = sqrt(3/2) * (σ_A - s_A)
    S_GT_min, S_GT_max = +s_A-θ_A, -s_A+θ_A
    return (J_GT_min, J_GT_max, S_GT_min, S_GT_max)
end;


# ---------------- Get Selection -----------------------#
function selction_coeff_LargeMat(r_vec, C_tot, Δx, Δμ)
    return s = (C_tot + diagm(0=>r_vec)) \ (Δx-Δμ)
end;

function selction_coeff_LargeMat_SL(r_vec, C_tot, Δx, Δμ)
    return s = (Δx-Δμ) ./ (C_tot[diagind(C_tot)] + r_vec)
end;


function selction_coeff_cholfact(mu, drift, C_tot, x1_init, x1_fini)
    return s = inv(cholesky(C_tot'*C_tot)) * C_tot' * (x1_fini-x1_init-mu*drift)
end;

function selction_coeff_SL(mu, drift, C_tot, x1_init, x1_fini)
    return s =  (x1_fini-x1_init-mu*drift) ./ C_tot[diagind(C_tot)]
end;


function get_x3(sample_t, n_t, L)
    Neff = Int(sum(n_t))
    LLhalf = Int(L*(L-1)/2)
    n_species = size(sample_t,2);
    x3 = zeros(L, L, L) 

    for m in 1:n_species
        n_t_scale = n_t[m] / Neff
        for i in 1:L
            x3[:, :, i] += ( sample_t[i, m] * n_t_scale * sample_t[1:L, m] ) * sample_t[1:L, m]'
        end
    end
    return x3
end;

function get_3rd_order_UFE(L, x1, x2, x3)
    UFE_3rd = zeros(L, L)
    for i in 1:L
        for j in (i+1):L
            if(x2[i,j]>1e-2) # if x2[i,j] is smaller than 0.01, we cant acuretly detect. 
                vec_UFE_3rd = zeros(L)
                for k in 1:L 
                    if(i!=k & j!=k)
                        num_num = x2[i,j] - x3[i,j,k] + 1e-3
                        num_den = 1 - x3[i,j,k] + (x2[i,j]+x2[j,k]+x2[k,i]) - (x1[i]+x1[j]+x1[k]) + 1e-3
                        log_num = log(abs(num_num / num_den))
                        #
                        den_num = (x1[j]+x3[i,j,k]-(x2[i,j]+x2[j,k])) * (x1[i]+x3[i,j,k]-(x2[i,j]+x2[i,k])) + 1e-3
                        log_den = log(abs(den_num / num_den^2))
                        if(abs(log_den)>1e-5)
                            vec_UFE_3rd[k] = log_num / log_den
                        end
                    end 
                end
                UFE_3rd[i,j] = minimum(vec_UFE_3rd)
            end
        end
    end
    return UFE_3rd
end;

function get_UFE_2nd(L, x1, x2, threshold)
    UFE = zeros(L, L)
    for i in 1:L
        for j in (i+1):L
            x11 = x2[i,j]
            x01 = x1[j] - x2[i,j]
            x10 = x1[i] - x2[i,j]
            x00 = 1 + x2[i,j] - x1[i] - x1[j]
            #@printf("i:%d, j:%d x00:%.1e x01:%.1e x10:%.1e x11:%.1e\n", i, j, x00, x01, x10, x11  )
            if(x00>threshold && x01>threshold && x10>threshold && x11>threshold)
                UFE[i,j] = 1 - log(x11 / x00) / log( (x01 * x10) / (x00 * x00))
            end
        end
    end
    return UFE 
end;


function get_averaged_Z(L, x1, x2, threshold)
    Z_value = zeros(L, L)
    Z_av_value = 0
    n_count = 0
    for i in 1:L
        for j in (i+1):L
            x11 = x2[i,j]
            x01 = x1[j] - x2[i,j]
            x10 = x1[i] - x2[i,j]
            x00 = 1 + x2[i,j] - x1[i] - x1[j]
            #@printf("i:%d, j:%d x00:%.1e x01:%.1e x10:%.1e x11:%.1e\n", i, j, x00, x01, x10, x11  )
            if(x00>threshold && x01>threshold && x10>threshold && x11>threshold)
                Z_temp = (x11 / x00) / (x01 * x10)
                Z_value[i,j] = Z_temp; Z_value[j,i] = Z_temp
                Z_av_value += Z_temp
                n_count += 1
            end
        end
    end
    if(n_count>0) Z_av_value /= n_count end 
    
    return (Z_value, Z_av_value)
end;


function get_averaged_LD(L, x1, x2, J_GT, threshold)
    LD_av_value = 0
    r2_tot, r2_neu, r2_pos, r2_neg = 0.0, 0.0, 0.0, 0.0
    n_count = 0
    n_neu, n_pos, n_neg = 0, 0, 0
    C = x2 - x1 * x1'
    LD_value = zeros(L, L)
    
    for i in 1:L
        for j in (i+1):L
            if((threshold < x1[i] < 1-threshold) && (threshold < x1[j] < 1-threshold))
                LD_temp = C[i,j] / sqrt(x1[i] * (1 - x1[i]) * x1[j] * (1 - x1[j]))  
                LD_value[i, j] = LD_temp; LD_value[j, i] = LD_temp
                LD_av_value += LD_temp
                n_count += 1
		
		r2_tot += LD_temp^2

		if(J_GT[i,j]> 0) r2_pos += LD_temp^2; n_pos += 1 end
		if(J_GT[i,j]==0) r2_neu += LD_temp^2; n_neu += 1 end
		if(J_GT[i,j]< 0) r2_neg += LD_temp^2; n_neg += 1 end
            end
        end
    end
    if(n_count>0) LD_av_value /= n_count; r2_tot /= n_count end 
    if(n_pos>0) r2_pos /= n_pos end
    if(n_neu>0) r2_neu /= n_neu end
    if(n_neg>0) r2_neg /= n_neg end
    return (LD_value, LD_av_value, r2_tot, r2_neu, r2_pos, r2_neg)
end;


""" get_PPV(mat_in, positive_epis_index, negative_epis_index, neutral_index, n_top_x_predition = 100)
"""
function get_PPV(mat_in, positive_epis_index, negative_epis_index, n_top_x_predition = 100)
    sorted_pos = sortperm(vec(mat_in), rev=true)
    sorted_neg = sortperm(vec(mat_in))
    n_TP_pos = sum(vec(positive_epis_index)[sorted_pos][1:n_top_x_predition])
    n_TP_neg = sum(vec(negative_epis_index)[sorted_neg][1:n_top_x_predition])
    PPV_pos = n_TP_pos/n_top_x_predition
    PPV_neg = n_TP_neg/n_top_x_predition
    return (PPV_pos, PPV_neg);
     
end;

""" get_PPV_Sensitivity(mat_in, n_top_types, L, α_lsit, 
    positive_epis_index, negative_epis_index)
"""
function get_PPV_Sensitivity(mat_in, n_top_types, L, α_lsit, 
    positive_epis_index, negative_epis_index)
p_pos = 1.0 / 3.0
accuracy_out = zeros(n_top_types, 2)   
sensitivity_out = zeros(n_top_types, 2)   
n_TP_pos_tot = sum(vec(positive_epis_index))
n_TP_neg_tot = sum(vec(negative_epis_index))

sorted_pos = sortperm(vec(mat_in), rev=true)
sorted_neg = sortperm(vec(mat_in))

for i_x in 1:n_top_types
    α_top_x_predition = α_lsit[i_x]

    n_top_eff = Int( floor(L^2  * (α_top_x_predition * p_pos) ) )
    n_TP_pos = sum(vec(positive_epis_index)[sorted_pos][1:n_top_eff])
    n_TP_neg = sum(vec(negative_epis_index)[sorted_neg][1:n_top_eff])

    accuracy_out[i_x, 1] = n_TP_pos/n_top_eff
    accuracy_out[i_x, 2] = n_TP_neg/n_top_eff
    #
    sensitivity_out[i_x, 1] = (n_TP_pos / n_TP_pos_tot) / (α_top_x_predition * p_pos); 
    sensitivity_out[i_x, 2] = (n_TP_neg / n_TP_neg_tot) / (α_top_x_predition * p_pos)
end
    return (accuracy_out, sensitivity_out)
end;


""" get_PPV_Sensitivity_distant_dependent(L, n_top_types, mat_in, dist_i_j, α_lsit, 
    positive_epis_index, negative_epis_index)
"""
function get_PPV_Sensitivity_distant_dependent(L, n_top_types, mat_in, dist_i_j, α_lsit, 
    positive_epis_index, negative_epis_index)
PPV_pos_dist_depend = zeros(n_top_types, L-1)
PPV_neg_dist_depend = zeros(n_top_types, L-1)
PPV_pos_dist_depend_within = zeros(n_top_types, L-1)
PPV_neg_dist_depend_within = zeros(n_top_types, L-1)
#
sensitivity_pos_dist_depend = zeros(n_top_types, L-1)
sensitivity_neg_dist_depend = zeros(n_top_types, L-1)
sensitivity_pos_dist_depend_within = zeros(n_top_types, L-1)
sensitivity_neg_dist_depend_within = zeros(n_top_types, L-1)
#
p_pos = 1.0 / 3
for i in 1:(L-1)
#        n =  i
    index_to_look = dist_i_j .== i
    index_to_look_within = dist_i_j .<= i
    sorted_pos = sortperm(vec(mat_in[index_to_look]), rev=true)
    sorted_neg = sortperm(vec(mat_in[index_to_look]));
    sorted_pos_within = sortperm(vec(mat_in[index_to_look_within]), rev=true)
    sorted_neg_within = sortperm(vec(mat_in[index_to_look_within]));

    n_TP_pos_tot = sum(vec(positive_epis_index[index_to_look]))
    n_TP_neg_tot = sum(vec(negative_epis_index[index_to_look]))
    n_TP_pos_within_tot = sum(vec(positive_epis_index[index_to_look_within]))
    n_TP_neg_within_tot = sum(vec(negative_epis_index[index_to_look_within]))

    for i_x in 1:n_top_types
        α_top_x_predition = α_lsit[i_x] #accuracy_top_x_set[i_x] / accuracy_top_x_set[end]
        n_top_eff = Int( floor(count(index_to_look)* (α_top_x_predition * p_pos) ) ) # Uniform probability of TP = 1/3
        n_top_eff_within = Int(floor(count(index_to_look_within) * (α_top_x_predition * p_pos))) # Uniform probability of TP = 1/3
        if(n_top_eff<1) n_top_eff = 1 end
        if(n_top_eff_within<1) n_top_eff_within = 1 end

        n_TP_pos = sum(vec(positive_epis_index[index_to_look])[sorted_pos][1:n_top_eff])
        n_TP_neg = sum(vec(negative_epis_index[index_to_look])[sorted_neg][1:n_top_eff])

        n_TP_pos_within = sum(vec(positive_epis_index[index_to_look_within])[sorted_pos_within][1:n_top_eff_within])
        n_TP_neg_within = sum(vec(negative_epis_index[index_to_look_within])[sorted_neg_within][1:n_top_eff_within])

        PPV_pos_dist_depend[i_x, i] = n_TP_pos / n_top_eff; 
        PPV_neg_dist_depend[i_x, i] = n_TP_neg / n_top_eff
        PPV_pos_dist_depend_within[i_x, i] = n_TP_pos_within / n_top_eff_within; 
        PPV_neg_dist_depend_within[i_x, i] = n_TP_neg_within / n_top_eff_within

        sensitivity_pos_dist_depend[i_x, i] = (n_TP_pos / n_TP_pos_tot) / (α_top_x_predition * p_pos)
        sensitivity_neg_dist_depend[i_x, i] = (n_TP_neg / n_TP_neg_tot) / (α_top_x_predition * p_pos)
        sensitivity_pos_dist_depend_within[i_x, i] = (n_TP_pos_within / n_TP_pos_within_tot) / (α_top_x_predition * p_pos)
        sensitivity_neg_dist_depend_within[i_x, i] = (n_TP_neg_within / n_TP_neg_within_tot) / (α_top_x_predition * p_pos)
    end
end

return (PPV_pos_dist_depend, PPV_neg_dist_depend, 
        PPV_pos_dist_depend_within, PPV_neg_dist_depend_within, 
        sensitivity_pos_dist_depend, sensitivity_neg_dist_depend, 
        sensitivity_pos_dist_depend_within, sensitivity_neg_dist_depend_within )
end;

function get_expected_changes_in_moment(LLhalf, L, S_GT, J_GT, x1, x2, x_1w2_set, mu, d_mu, rec, d_rec)
    S_J_true_flat = zeros(L+LLhalf)
    for i in 1:L
        S_J_true_flat[i] = S_GT[i]
        for j in (i+1):L
            ξ = G(i,j,L) + L
            S_J_true_flat[ξ] = J_GT[i,j]
        end
    end

    Cs_projected_temp = zeros(LLhalf+L)
    ΞS_projected = (x_1w2_set * S_J_true_flat)
    for m in 1:size(x_1w2_set[:, 1],1)
        Cs_projected_temp += x_1w2_set[m,:]　*　ΞS_projected[m]
    end
    expected_changes_in_1st_moment = Cs_projected_temp[1:L] - mu*d_mu[1:L] - rec*d_rec[1:L]
    expected_changes_in_2nd_moment = Cs_projected_temp[(L+1):end] - mu*d_mu[(L+1):end] - rec*d_rec[(L+1):end]
    denominator_dot_x_ij_temp = 0; numerator_dot_x_ij_temp = 0;
    denominator_dot_x_ij_scaled_temp = 0; numerator_dot_x_ij_scaled_temp = 0;
    
    for i in 1:L
        for j in (i+1):L
            ξ = G(i,j,L) 
            denominator_dot_x_ij_temp += (expected_changes_in_1st_moment[i] * expected_changes_in_1st_moment[j])^2
            numerator_dot_x_ij_temp += (expected_changes_in_2nd_moment[ξ] - (expected_changes_in_1st_moment[i] * x1[j] + expected_changes_in_1st_moment[j] * x1[i] ))^2
            
            denominator_dot_x_ij_scaled_temp += x2[i,j] * (expected_changes_in_1st_moment[i] * expected_changes_in_1st_moment[j])^2
            numerator_dot_x_ij_scaled_temp += x2[i,j] * (expected_changes_in_2nd_moment[ξ] - (expected_changes_in_1st_moment[i] * x1[j] + expected_changes_in_1st_moment[j] * x1[i] ))^2
        end
    end
    return (sqrt(denominator_dot_x_ij_temp), sqrt(numerator_dot_x_ij_temp), 
            sqrt(denominator_dot_x_ij_scaled_temp), sqrt(numerator_dot_x_ij_scaled_temp))
end;

function get_expected_changes_in_moment_discritized(L, LLhalf, x_1w2, x_1w2_old)
    denominator_dot_x_ij_temp, numerator_dot_x_ij_temp = 0, 0
    for i in 1:L
        for j in (i+1):L
            ξ = G(i,j,L) + L
            denominator_dot_x_ij_temp += ((x_1w2[i] - x_1w2_old[i]) * (x_1w2[j] - x_1w2_old[j]))^2
            numerator_dot_x_ij_temp += ((x_1w2[ξ] - x_1w2[i]*x_1w2[j]) - (x_1w2_old[ξ] - x_1w2_old[i]*x_1w2_old[j]))^2
        end
    end
    return (denominator_dot_x_ij_temp, numerator_dot_x_ij_temp)
end;


"""
    resampling_with_N(data, n_resample)
"""
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
