function ising_element_1D(i1::Int, i2::Int, β::Float64, J::Float64)
    σ1 = 2*i1 -1 
    σ2 = 2*i2 -1
    return exp(β * J * σ1 * σ2)
end
function ising_partition_1D_loop(N::Int, β::Float64, J::Float64; left::Int=1, right::Int=1)
    if N < 2
        error("N must be ≥ 2")
    elseif N == 2
        # 端しかないとき
        return ising_element_1D(left, right, β, J)
    end

    n_middle = N - 2
    total = 0.0

    # 全ての中間スピン構成（2^(N-2) 通り）をループ
    for n in 0:(2^n_middle - 1)
         # 2進数に変換 
        middle = digits(n, base=2, pad=n_middle)
        # フルなスピン配列
        config = [left; middle; right]
        # 重み（ボルツマン因子）の積を計算
        weight = 1.0
        for i in 1:(N-1)
            weight *= ising_element_1D(config[i], config[i+1], β, J)
        end
        total += weight
    end

    return total
end