
function get_pair_codon_freq(x_vec, y_vec)
    pair_nuc_freq = zeros(q,q)
    for (x,y) in zip(x_vec, y_vec)
        i_a, i_b = nuc2num[x], nuc2num[y]
        pair_nuc_freq[i_a, i_b] += 1
    end;
    return pair_nuc_freq /= length(x_vec)
end;

function get_models(sel_eps, sel_add);
    q = length(NUC);
    L = maximum(sel_add.i);
    model_add = zeros(q*L)
    model_eps = zeros(q*L, q*L); # Suppose the diagonal is additive selection, and off-diagonal elements of i <j are epistatic interactions. 
    
    for n in 1:length(sel_add.i)
        x = sel_add[n, :]
        v, i, a = x.s_MPL, x.i, x.a
        model_add[km(i,a,q)] = v
    end;
    model_temp = zeros(L, L)
    for n in 1:length(sel_eps.i)
        x = sel_eps[n, :]
        v, i, j, a, b = x.s_MPL, x.i, x.j, x.a, x.b
        model_eps[km(i,a,q), km(j,b,q)] = v
        model_temp[i,j] += abs(v)
    end;
    return (model_add, model_eps, model_temp)
end;

function get_epistatic_fit(q, L, csv_index, model_eps, seq_temp, dist_cut = 50)
    
    #f_epis_temp = sum( model_eps[km.(1:L, seq_temp, q), km.(1:L, seq_temp, q)] )
    f_epis_temp = 0
    for i in 1:L
        # Maybe better to consider only short-range interactions? 
        for j in (i+1):L 
            if(abs( Pol2HXB2[i] - Pol2HXB2[j] ) < dist_cut)
                f_epis_temp += model_eps[km(i, seq_temp[i], q), km(j, seq_temp[j], q)]
            end
        end
    end
    return f_epis_temp
end;

km(i,a,q) = (i-1) * q + a;

function i_is_within_K(hxb2_i, epitope_sites, K)
    bool_out = false
    for k in 0:K
        if(hxb2_i + k ∈ epitope_sites ) bool_out = true; return bool_out; break end;
        if(hxb2_i - k ∈ epitope_sites ) bool_out = true; return bool_out; break end;
    end;
    return bool_out
end;



function gauge_TF_zero_epistatic(q, L, seq_TF_nuc, s_add_input, s_epis_input)
    qL = q*L
    qqLLhalf = Int(q^2*L*(L-1)*0.5);
    s_add_gauged = zeros(qL)
    s_epis_gauged = zeros(qL, qL)
    
    s_ij_aq_vec, s_ij_qa_vec = zeros(qL), zeros(qL)
    for i in 1:L
        a_TF = seq_TF_nuc[i]
        x_ia_TF = km(i,a_TF,q)
        for j in (i+1):L # G function should be i<j
            b_TF = seq_TF_nuc[j]
            x_jb_TF = km(j,b_TF,q)
            for a in 1:q 
                x_ia = km(i,a,q)
                s_ij_aq_vec[x_ia] += s_epis_input[x_ia, x_jb_TF] - s_epis_input[x_ia_TF, x_jb_TF]
            end
            for b in 1:q
                x_jb = km(j,b,q)
                s_ij_aq_vec[x_jb] += s_epis_input[x_ia_TF, x_jb] - s_epis_input[x_ia_TF, x_jb_TF]
            end
        end
    end
    
    for i in 1:L
        a_TF = seq_TF_nuc[i]
        x_ia_TF = km(i,a_TF,q)
        for a in 1:q
            x_ia = km(i,a,q)
            s_add_gauged[x_ia] = s_add_input[x_ia] - s_add_input[x_ia_TF] + s_ij_aq_vec[x_ia] + s_ij_qa_vec[x_ia]
            for j in (i+1):L
                b_TF = seq_TF_nuc[j]
                x_jb_TF = km(j,b_TF,q)
                for b in 1:q
                    x_jb = km(j,b,q)
                    s_iajb_temp = s_epis_input[x_ia, x_jb] - s_epis_input[x_ia, x_jb_TF] - s_epis_input[x_ia_TF, x_jb] + s_epis_input[x_ia_TF, x_jb_TF]
                    s_epis_gauged[x_ia, x_jb] = s_iajb_temp
                    s_epis_gauged[x_jb, x_ia] = s_iajb_temp
                end
            end
        end
    end
    return (s_add_gauged, s_epis_gauged)
end;

function gauge_TF_zero_additive(q, L, seq_TF_nuc, s_add_input)
    qL = q*L
    s_add_gauged = zeros(qL)
    for i in 1:L
        a_TF = seq_TF_nuc[i]
        x_ia_TF = km(i,a_TF,q)
        for a in 1:q
            x_ia = km(i,a,q)
            s_add_gauged[x_ia] = s_add_input[x_ia] - s_add_input[x_ia_TF]
        end
    end
    return s_add_gauged
end;;

# Function to extract integer from a string
function extract_integer(str)
    # Use a regular expression to find all sequences of digits
    matches = match(r"\d+", str)

    # If there are matches, convert the first match to an integer
    if matches !== nothing
        return parse(Int, matches.match)
    else
        return nothing # or some default value/error handling
    end
end;


# The following function should be in a script. 
function get_TP_FP(bool_vec, idx_sort)
    TP_pos = cumsum(bool_vec[idx_sort])
    FP_pos = cumsum(.!bool_vec[idx_sort])
    return (TP_pos/TP_pos[end], FP_pos/FP_pos[end])
end;

function get_ROC_AUC(TP, FP)
    myauc = 0.0 
    total_positive = maximum(TP)
    total_negative = maximum(FP)
    #Assume dx is zero at n=1.
    for n in 2:length(TP)
        # Trapezoid integration
        y = 0.5 * (TP[n] + (TP[n-1]) )/ total_positive
        dx = (FP[n] - FP[n-1]) / total_negative
        x = FP[n]
        myauc += y*dx
        #@printf("n=%d, x=%.2f, y=%.2f\n",  n, x/total_negative, auc_MPL_ben)
    end
    return myauc
end;

function count_files_with_key(directory::String, file_key::String)
    # Check if the directory exists
    if !isdir(directory)
        error("Directory does not exist")
    end
    
    # Create a pattern to match files containing the file key
    pattern = joinpath(directory, "*$file_key*")
    
    # Find all files that match the pattern
    files = glob(pattern)
    
    # Return the count of files
    return length(files)
end

