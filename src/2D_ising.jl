function ising_transfer_element_2D(i1::Int,i2::Int,j1::Int,j2::Int,
    β::Float64, J::Float64)
    σ1 = 2*i1 -1 #left spin
    σ2 = 2*i2 -1 #right spin
    s1 = 2*j1 -1 #down spin
    s2 = 2*j2 -1 #up spin
    return exp(β*J*(σ1*s1 + s1*σ2 + σ2*s2 + s2*σ1))
end
#2Dイジングの転送行列を定義する関数
function ising_transfer_tensor_2D(β::Float64, J::Float64,
    σ1::Index, σ2::Index, s1::Index, s2::Index)
    T = ITensor(σ1, σ2, s1, s2)
    for i1 in 0:1, i2 in 0:1,
        j1 in 0:1, j2 in 0:1
        T[σ1=>i1+1, σ2=>i2+1, s1=>j1+1, s2=>j2+1] =
            ising_transfer_element_2D(i1, i2, j1, j2, β, J)
    end
    return T
end
function make_boundary_tensors(L::Int, ind_list::Vector{<:Index})#Indexの全ての方を受け入れる
    l = length(ind_list)
    @assert l == 2*L "Length of ind_list must match L"
    boundary_list = Vector{ITensor}(undef, l)
    for (i, ind) in enumerate(ind_list)
        boundary_list[i] = fixend(ind)
    end
    return boundary_list
end
function make_inner_tensors(L::Int, β::Float64, J::Float64,
    sind_list::Array{<:Index},σind_list::Array{<:Index})
    length = 2*L
    T_list = Array{ITensor}(undef, length, length)
    for i in 1:length, j in 1:length
        T_list[i, j] = ising_transfer_tensor_2D(β, J,
            σind_list[j,i], σind_list[j+1,i],
            sind_list[j,i], sind_list[j,i+1])
    end
    return T_list
end 
#2Dの分配関数を作る関数
function setup_tensors_2D(L::Int, β::Float64, J::Float64)
    σ_list = make_2D_index_σtags(L)
    s_list = make_2D_index_stags(L)
    T_list = make_inner_tensors(L, β, J, s_list, σ_list)
    Dup_list = make_boundary_tensors(L, σ_list[2*L+1,:]) #上の境界
    Ddown_list = make_boundary_tensors(L, σ_list[1,:]) #下の境界
    L_list = make_boundary_tensors(L, s_list[:,1]) #左の境界
    R_list = make_boundary_tensors(L, s_list[:,2*L+1]) #右の境界
    return T_list, Dup_list, Ddown_list, L_list, R_list
end
function partition_function_2D(L::Int, β::Float64, J::Float64)
    T_list, Dup_list, Ddown_list, L_list, R_list = setup_tensors_2D(L, β, J)
    result = ITensor(true)
    # 内部テンソルをすべて掛ける
    for T in T_list
        result *= T
    end
    # 上端
    for T in Dup_list
        result *= T
    end
    # 左端
    for T in L_list
        result *= T
    end
    # 右端
    for T in R_list
        result *= T
    end
    # 下端
    for T in Ddown_list
        result *= T
    end
    return only(result)
end
function contruct_with_boundary(L::Int, β::Float64, J::Float64)
    T_list, Dup_list, Ddown_list, L_list, R_list = setup_tensors_2D(L, β, J)
    T1_list = copy(T_list) #テンソルのコピーを作成
    # 上端の境界条件を適用
    for (i, Ddown) in enumerate(Ddown_list)
        T1_list[i,1] = T_list[i,1] * Ddown
    end
    # 左端の境界条件を適用
    for (i, L) in enumerate(L_list)
        ord = order(T1_list[1,i])
        if ord == 0
            T1_list[1,i] = T_list[1,i] * L
        else
            T1_list[1,i] = T1_list[1,i] * L
        end
    end
    # 右端の境界条件を適用
    for (i, Dup) in enumerate(Dup_list)
        ord = order(T1_list[i,2*L])
        if ord == 0
            T1_list[i,2*L] = T_list[i,2*L] * Dup
        else
            T1_list[i,2*L] = T1_list[i,2*L] * Dup
        end
    end
    # 下端の境界条件を適用
    for (i, R) in enumerate(R_list)
        ord = order(T1_list[2*L,i])
        if ord == 0
            T1_list[2*L,i] = T_list[2*L,i] * R
        else
            T1_list[2*L,i] = T1_list[2*L,i] * R
        end
    end
    return T_list, T1_list
end
#期待値を取るサイトのインデックスを貼り替える関数
function replace_tag_ITensor(T::ITensor, tag::String)
    target = findfirst(ind -> hastags(ind, "$tag"), inds(T))
    myind = inds(T)[target]
    newind = prime(myind)
    T2 = replaceind(T, myind, newind)
    return T2, newind
end
#スピン期待値を計算する関数
function calc_spontaneous_magnetization(L::Int, β::Float64, J::Float64
    ,tag::String)
    Z = partition_function_2D(L, β, J)
    _, T1_list = contruct_with_boundary(L, β, J)
    # 上側と下側に分ける
    result_upper = ITensor(true)
    result_lower = ITensor(true)
    for T in T1_list[1:L,:]
        result_upper *= T
    end
    for T in T1_list[L+1:end,:]
        result_lower *= T
    end
    target = findfirst(ind -> hastags(ind, "$tag"), inds(result_upper))
    myind = inds(result_upper)[target]
    T2, newind = replace_tag_ITensor(result_upper, "$tag")
    σ_tensor = spin_tensor(myind, newind)
    σ = T2 * σ_tensor 
    return only(result_lower * σ) / Z
end
#スピン相関関数を計算するコード
function calc_spin_correlation(L::Int, β::Float64, J::Float64,
    tag1::String, tag2::String)
    Z = partition_function_2D(L, β, J)
    T_list, T1_list = contruct_with_boundary(L, β, J)
    # 上側と下側に分ける
    result_upper = ITensor(true)
    result_lower = ITensor(true)
    for T in T1_list[1:L,:]
        result_upper *= T
    end
    for T in T1_list[L+1:end,:]
        result_lower *= T
    end
    target1 = findfirst(ind -> hastags(ind, "$tag1"), inds(result_upper))
    myind1 = inds(result_upper)[target1]
    target2 = findfirst(ind -> hastags(ind, "$tag2"), inds(result_upper))
    myind2 = inds(result_upper)[target2]
    T2_1, newind1 = replace_tag_ITensor(result_upper, "$tag1")
    T2_2, newind2 = replace_tag_ITensor(T2_1, "$tag2")
    σ_tensor_1 = CTMRG.spin_tensor(myind1, newind1)
    σ_tensor_2 = CTMRG.spin_tensor(myind2, newind2)

    σσ = T2_2 * σ_tensor_1 * σ_tensor_2 
    return only(result_lower * σσ) / Z
end