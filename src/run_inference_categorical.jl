using DelimitedFiles
using Distributed
using Random
using Distributions
using LinearAlgebra
using Printf
using CSV

using DataFrames

rng = Random.MersenneTwister(1234);

include("./systematic_performance_evaluation.jl");
include("./efficient_epistasis_MPL_categorical.jl");
γ_A, γ_E =  parse(Int, ARGS[1]), parse(Int, ARGS[2]) #5, 1e10 # Regularization for additive and epistatic fitness
file_key = ARGS[3] # "700010058-3"
file_dir = ARGS[4] # /net/dali/home/barton/kais/paper-MPL-inference-master/src/HIV
file_dir2 = ARGS[5] # /net/dali/home/barton/kais/paper-MPL-inference-master/data/HIV/processed
dir_out = ARGS[6] #
fname_mut = ARGS[7] #
flag_week = parse(Bool, ARGS[8]) # false
flag_reg_E_is_L = parse(Bool, ARGS[9]) # true or false
flag_mask_gaps = parse(Bool, ARGS[10]) # true or false
Lc = parse(Int, ARGS[11]) # true or false

#fname_mut = "../data_2024/Zanini-extended.dat"
rec = 1e-5
muMat = readdlm(fname_mut);

fname_seq = @sprintf("%s/%s-poly-seq2state.dat", file_dir, file_key)
data_raw = readdlm(fname_seq, Int);

# ---------------------------- Temporal treatment ---------------------------- #
# This is a temporal treatment only for CH848; taking evolution upto 4 years.
#idx = data_raw[:, 1] .< 356 * 4.0
#data_raw = data_raw[idx, :]
# ----------------------------------------------------------------------------#

if(flag_week)
    data_raw[:,1] *= 7;
end;
data_raw[:, 3:end] .+= 1;
time_list = data_raw[:, 1]
# this csv file is for computing the linear distance on HXB2 sequence. 
csv_index = CSV.read(@sprintf("%s/%s-index.csv", file_dir2, file_key), DataFrame);

#csv_index.alignment, csv_index.polymorphic
idx = csv_index.polymorphic .!= "NA"
seq_TF_nuc = [NUC2NUM[x] for x in csv_index.TF[idx]];
indicies_mapped = Dict(a => b for (a, b) in zip(parse.(Int, csv_index.polymorphic[idx]) .+ 1, extract_integer.(string.(csv_index.HXB2[idx]))) );

L = length(data_raw[1, 3:end])
Lc = minimum([L, Lc])
q = length(unique(data_raw[:, 3:end]));
qL = q*L
qqLLhalf = Int(q^2*L*(L-1)*0.5);
rank_x = qL+qqLLhalf;

time_list_unique = sort(unique(time_list))
n_time_list_max = length(time_list_unique)
time_list = data_raw[:, 1]
idx_sort = sortperm(time_list)
time_list = time_list[idx_sort]
data = copy(data_raw[idx_sort, :]);

if(flag_reg_E_is_L)
    @printf("Regularization for epistatic fitness is set to L.\n")
    γ_E = L
end;

@time (Ξ, d_mu_tot, d_rec_tot, Δx_tot, sites_to_look) = get_Ξ_frequency_change_fast(rank_x, q, L, time_list_unique, data, muMat, indicies_mapped);
#γ_vec = [γ_A * ones(count(sites_to_look[1:qL])); γ_E * ones(count(sites_to_look[(1+qL):end]))];
#reg_short, reg_long = γ_E*Lc, 5*γ_E*L # Regularization for additive and epistatic fitness
reg_short, reg_long = γ_E, 1e10 # Regularization for additive and epistatic fitness
#γ_E_block_dependent = get_strong_reg_in_block(reg_short, reg_long, Lc, q, L)
γ_E_block_dependent = get_strong_reg_in_block_HXB2(csv_index, reg_short, reg_long, Lc, q, L)
γ_vec = [γ_A * ones(count(sites_to_look[1:qL])); γ_E_block_dependent[sites_to_look[(1+qL):end]]];
if(flag_mask_gaps)
    γ_vec = [γ_A * ones(qL); γ_E_block_dependent];
    site_gaps = get_index_involved_in_gaps(q, L)
    γ_vec[Bool.(site_gaps)] .= reg_long
    γ_vec = copy(γ_vec[sites_to_look])
end

# ---------------------------- MPL -------------------------------- #
s_MPL, elapsed_time, bytes_allocated, gc_time, memory_allocated = @timed begin
    s_MPL = get_selection_efficient(Δx_tot, d_mu_tot, rec*d_rec_tot, Ξ, γ_vec)
end;
# ----------------------------- SL -------------------------------- #
s_SL = get_selection_efficient_SL(Δx_tot, d_mu_tot, rec*d_rec_tot, Ξ, γ_vec);

# ------------------------ Output ot computational time ------------------------- #
fout = open( @sprintf("%s/computational_time_efficient-EpisMPL_%s.txt", dir_out, file_key) , "w")
println(fout, @sprintf("Elapsed_time\t%.3e", elapsed_time)) # 354.607983826 seconds --> for uniform Gamma, arbitrary diagonal matrix case is 373.22659798 seconds
#println(fout, "Bytes_allocated $bytes_allocated")
println(fout, @sprintf("GC_time\t%.3e", gc_time))
total_bytes = memory_allocated.allocd
bytes_to_gib = total_bytes / 2^30;
println(fout, @sprintf("Memory_allocated\t%.3e", bytes_to_gib)); # GB
println(fout, @sprintf("Number_of_sites\t%d", count(sites_to_look)))
close(fout)
# ---------------------------------------------------------------------------------- #
s_MPL_copy = zeros(rank_x); s_MPL_copy[sites_to_look] = s_MPL
s_SL_copy = zeros(rank_x); s_SL_copy[sites_to_look] = s_SL;
Δx_tot_copy = zeros(rank_x); Δx_tot_copy[sites_to_look] = Δx_tot;
d_mu_tot_copy = zeros(rank_x); d_mu_tot_copy[sites_to_look] = d_mu_tot;
d_rec_tot_copy = zeros(rank_x); d_rec_tot_copy[sites_to_look] = d_rec_tot;

fname_out = @sprintf("%s/selection_efficient-EpisMPL_%s.txt", dir_out, file_key )
i_set, j_set, a_set, b_set, s_MPL_set_out, s_SL_set_out, Δx_tot_set_out, d_mu_tot_set_out, d_rec_tot_set_out = [], [], [], [], [], [], [], [], [];
for i in 1:L
    for a in 1:q
        u = km(i,a,q)
        push!(i_set, i); push!(j_set, i); push!(a_set, a); push!(b_set, a);
        push!(s_MPL_set_out, s_MPL_copy[u]); push!(s_SL_set_out, s_SL_copy[u]); push!(Δx_tot_set_out, Δx_tot_copy[u]); 
        push!(d_mu_tot_set_out, d_mu_tot_copy[u]); push!(d_rec_tot_set_out, d_rec_tot_copy[u]);
    end
end
for i in 1:L
    for j in (i+1):L
        for a in 1:q
            for b in 1:q
                ξ = G(i,j,a,b,q,L) + qL
                if(sites_to_look[ξ])
                    push!(i_set, i); push!(j_set, j); push!(a_set, a); push!(b_set, b);
                    push!(s_MPL_set_out, s_MPL_copy[ξ]); push!(s_SL_set_out, s_SL_copy[ξ]); push!(Δx_tot_set_out, Δx_tot_copy[ξ]);
                    push!(d_mu_tot_set_out, d_mu_tot_copy[ξ]); push!(d_rec_tot_set_out, d_rec_tot_copy[ξ]);
                end
            end
        end
    end
end
df = DataFrame(
    i = i_set, j = j_set, a = a_set, b = b_set,
    s_MPL = s_MPL_set_out, s_SL = s_SL_set_out,
    Δx_tot = Δx_tot_set_out, d_mu_tot = d_mu_tot_set_out, d_rec_tot = d_rec_tot_set_out
)
CSV.write(fname_out, df)
df = []; i_set = []; j_set = []; a_set = []; b_set = []; s_MPL_set_out = []; s_SL_set_out = []; Δx_tot_set_out = []; d_mu_tot_set_out = []; d_rec_tot_set_out = [];
# ---------------------------------------------------------------------------------- # 

# ----------------------------- Copying the selection coefficients and other aceccary vectors ---------------- #
idx = collect(1:qL)[sites_to_look[1:qL]];
Δx_tot_additive = copy(Δx_tot[1:count(sites_to_look[1:qL])])
d_mu_tot_additive = copy(d_mu_tot[1:count(sites_to_look[1:qL])])
d_rec_tot_additive = copy(d_rec_tot[1:count(sites_to_look[1:qL])])
Ξ_mini = copy(Ξ[1:count(sites_to_look[1:qL]), :]);
γ_vec_add = γ_vec[1:count(sites_to_look[1:qL])]
#---------- MPL only additive ------------#
s_MPL_additive, elapsed_time, bytes_allocated, gc_time, memory_allocated = @timed begin
    #s_MPL_additive = get_selection_efficient(Δx_tot_additive, d_mu_tot_additive, rec*d_rec_tot_additive, Ξ_mini, γ_A * ones(length(idx)) )
    s_MPL_additive = get_selection_efficient(Δx_tot_additive, d_mu_tot_additive, rec*d_rec_tot_additive, Ξ_mini, vec(γ_vec_add) )
end;
# -------- SL ----------#
s_SL_additive = get_selection_efficient_SL(Δx_tot_additive, d_mu_tot_additive, rec*d_rec_tot_additive, Ξ_mini, vec(γ_vec_add) );
s_MPL_additive_copy = zeros(qL); s_MPL_additive_copy[idx] = s_MPL_additive
s_SL_additive_copy = zeros(qL); s_SL_additive_copy[idx] = s_SL_additive;
Δx_tot_additive_copy = zeros(qL); Δx_tot_additive_copy[idx] = Δx_tot_additive;
d_mu_tot_additive_copy = zeros(qL); d_mu_tot_additive_copy[idx] = d_mu_tot_additive;
# ------------------------ Output ot computational time ------------------------- #
fout = open( @sprintf("%s/computational_time_efficient-AdditMPL_%s.txt", dir_out, file_key) , "w")
println(fout, @sprintf("Elapsed_time\t%.3e", elapsed_time)) # 354.607983826 seconds --> for uniform Gamma, arbitrary diagonal matrix case is 373.22659798 seconds
#println(fout, "Bytes_allocated $bytes_allocated")
println(fout, @sprintf("GC_time\t%.3e", gc_time))
total_bytes = memory_allocated.allocd
bytes_to_gib = total_bytes / 2^30;
println(fout, @sprintf("Memory_allocated\t%.3e", bytes_to_gib)); # GB
println(fout, @sprintf("Number_of_sites\t%d", count(sites_to_look)))
close(fout)
# ---------------------------------------------------------------------------------- #
fname_out = @sprintf("%s/selection_efficient-AdditMPL_%s.txt", dir_out, file_key) 
i_set, a_set, s_MPL_set_out, s_SL_set_out, Δx_tot_set_out, d_mu_tot_set_out = [], [], [], [], [], [];
for i in 1:L
    for a in 1:q
        u = km(i,a,q)
        push!(i_set, i); push!(a_set, a);
        push!(s_MPL_set_out, s_MPL_additive_copy[u]); push!(s_SL_set_out, s_SL_additive_copy[u]); 
        push!(Δx_tot_set_out, Δx_tot_additive_copy[u]); push!(d_mu_tot_set_out, d_mu_tot_additive_copy[u]);
    end
end
df = DataFrame(
    i = i_set, a = a_set,
    s_MPL = s_MPL_set_out, s_SL = s_SL_set_out,
    Δx_tot = Δx_tot_set_out, d_mu_tot = d_mu_tot_set_out
)
CSV.write(fname_out, df)
df = []; i_set = []; a_set = []; s_MPL_set_out = []; s_SL_set_out = []; Δx_tot_set_out = []; d_mu_tot_set_out = [];
# ---------------------------------------------------------------------------------- # 
