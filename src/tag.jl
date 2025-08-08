function sigma_tensor(ind1::Index, ind2::Index)
    return diagITensor([1.0, -1.0], ind1, ind2)
end
function make_1D_index_tags(s::String, N::Int)
    s_list = Vector{Index}(undef, m) #undefで初期化
    for i in 1:N
        s_list[i] = Index(2, "$s,m=$i")
    end
    return s_list
end
function make_2D_index_stags(L::Int)
    i = 2*L + 1 # 列
    j = 2*L # 行
    sind_list = Array{Index}(undef, j, i) # j行 × i列
    for m in 1:i # 列優先でループを回す
        for n in 1:j
            sind_list[n, m] = Index(2, "s,j=$n,i=$m")
        end
    end
    return sind_list
end
function make_2D_index_σtags(L::Int)
    i = 2*L # 列
    j = 2*L + 1 # 行
    σind_list = Array{Index}(undef, j, i) # j行 × i列
    for m in 1:i 
        for n in 1:j
            σind_list[n, m] = Index(2, "σ,j=$n,i=$m")
        end
    end
    return σind_list
end 
# tagの張り替えの関数 replaceindがあったので没
#=
function change_tag(T::ITensor, old_ind::Index, new_ind::Index)
    new_T = ITensor(old_ind, new_ind)
    new_T[old_ind => 1, new_ind => 1] = 1
    new_T[old_ind => 2, new_ind => 2] = 1
    new_T[old_ind => 1, new_ind => 2] = 0
    new_T[old_ind => 2, new_ind => 1] = 0
    return new_T*T
end
=#