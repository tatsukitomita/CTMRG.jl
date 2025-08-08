function ising_transfer_element_1D(i1::Int, i2::Int,β::Float64,J::Float64)
    s1 = 2*i1 - 1
    s2 = 2*i2 - 1
    expvalue = exp(β*J*s1*s2)
    return expvalue
end
function ising_transfer_matrix_1D(β::Float64, J::Float64,s1::Index, s2::Index)
    T = ITensor(s1, s2)
    for i in 0:1, j in 0:1
        T[s1=>i+1, s2=>j+1] = ising_transfer_element_1D(i, j, β, J)
    end
    return T
end
#1次元イジング模型の分配関数を計算する関数
function partition_function_1D(N::Int, β::Float64, J::Float64)
    s_list = Vector{Index}(undef, N) 
    T_list = Vector{ITensor}(undef, N-1)          # 転送行列を格納する配列
    boundary = Vector{ITensor}(undef, 2)
    for i in 1:N
        s_list[i] = Index(2, "s, n=$i")
    end
    for i in 1:N
        if i == 1
            boundary[1] = fixend(s_list[1])
            T_list[i] = ising_transfer_matrix_1D(β, J, s_list[i], s_list[i+1])
        elseif i == N
            boundary[2] = fixend(s_list[N])
        else 
            T_list[i] = ising_transfer_matrix_1D(β, J, s_list[i], s_list[i+1])
        end
    end
    return s_list, T_list, boundary
end
#分配関数の縮約をとる関数
function partition_function_contraction(N::Int)
    result = ITensor(true)
    s_list, T_list, boundary = partition_function_1D(N, β, J)
    L = boundary[1]
    R = boundary[2]
    for T in T_list
        result = result * T   # 順に縮約
    end
    return only(L * result * R) #順番は任意
end