#３本足テンソルを作るコード
function P_edge(β, J, σind1, σind2, sind1)
    dumyind = Index(2, "dummy")
    T = ising_transfer_tensor_2D(β, J, σind1, σind2, sind1, dumyind)
    D = fixend(dumyind)
    P = T*D
    return P
end
function P_corner(β, J, σind1, sind1, sind2)
    dumyind = Index(2, "dummy")
    T = ising_transfer_tensor_2D(β, J, σind1, dumyind, sind1, sind2)
    D = fixend(dumyind)
    P = T*D
    return P
end
#２本足のテンソルを作るコード
function C_corner(β, J, σind1, sind1)
    dumyind_σ = Index(2, "dummy_σ")
    dumyind_s = Index(2, "dummy_s")
    T = CTMRG.ising_transfer_tensor_2D(β, J, σind1, dumyind_σ, sind1, dumyind_s)
    D_σ = fixend(dumyind_σ)
    D_s = fixend(dumyind_s)
    C = T*D_σ*D_s
    return C
end
#tagの作成
function make_CTMRG_ind_stags(L::Int)
    i = 2 * L # 行
    j = 2 * L -1 # 列
    sind_list = Array{Index}(undef, i, j)
    for m in 1:j, n in 1:i
        sind_list[n, m] = Index(2, "s, i=$n, j=$m")
    end
    return sind_list
end
function make_CTMRG_ind_σtags(L::Int)
    i = 2 * L - 1 #　行
    j = 2 * L # 列
    σind_list = Array{Index}(undef, i, j)
    for m in 1:j, n in 1:i
        σind_list[n, m] = Index(2, "σ, i=$n, j=$m")
    end
    return σind_list
end
#1ループ目の角転送行列
function setup_leftdown_corner_first(β::Float64, J::Float64,
    σindlist::Matrix{Index},sindlist::Matrix{Index})
    cindσ = combiner(σindlist[2, 1], σindlist[2, 2], tags="σ_all")
    cinds = combiner(sindlist[1, 2], sindlist[2, 2], tags="s_all")
    #2脚テンソルの定義
    C = CTMRG.C_corner(β, J, σindlist[1, 1], sindlist[1, 1])
    #3脚テンソルの定義
    P_edge = CTMRG.P_edge(β, J, σindlist[1, 1], σindlist[2, 1], sindlist[2, 1])
    P_corner = CTMRG.P_corner(β, J, σindlist[1, 2], sindlist[1, 1], sindlist[1, 2])
    #4脚テンソルの定義
    W =  CTMRG.ising_transfer_tensor_2D(β, J, σindlist[1, 2], σindlist[2, 2], sindlist[2, 1], sindlist[2, 2])
    #縮約
    _P_edge = P_edge * W
    C_expand = P_corner * C * _P_edge
    return C_expand,cindσ, cinds
end
function setup_leftup_corner_first(β::Float64, J::Float64,
    σindlist::Matrix{Index},sindlist::Matrix{Index})
    cindσ = combiner(σindlist[2, 1], σindlist[2, 2], tags="σ_all")
    cinds = combiner(sindlist[3, 2], sindlist[4, 2], tags="s_all")
    #2脚テンソルの定義
    C = CTMRG.C_corner(β, J, σindlist[3, 1], sindlist[4, 1])
    #3脚テンソルの定義
    P_edge = CTMRG.P_edge(β, J, σindlist[3, 1], σindlist[2, 1], sindlist[3, 1])
    P_corner = CTMRG.P_corner(β, J, σindlist[3, 2], sindlist[4, 1], sindlist[4, 2])
    #4脚テンソルの定義
    W =  CTMRG.ising_transfer_tensor_2D(β, J, σindlist[3, 2], σindlist[2, 2], sindlist[3, 1], sindlist[3, 2])
    #縮約
    _P_edge = P_edge * W
    C_expand = P_corner * C * _P_edge
    return C_expand, cindσ, cinds
end
#右上を作成
function setup_rightup_corner_first(β::Float64, J::Float64,
    σindlist::Matrix{Index},sindlist::Matrix{Index})
    cindσ = combiner(σindlist[2, 3], σindlist[2, 4], tags="σ_all")
    cinds = combiner(sindlist[3, 2], sindlist[4, 2], tags="s_all")
    #2脚テンソルの定義
    C = CTMRG.C_corner(β, J, σindlist[3, 4], sindlist[4, 3])
    #3脚テンソルの定義
    P_edge = CTMRG.P_edge(β, J, σindlist[3, 4], σindlist[2, 4], sindlist[3, 3])
    P_corner = CTMRG.P_corner(β, J, σindlist[3, 3], sindlist[4, 2], sindlist[4, 3])
    #4脚テンソルの定義
    W =  CTMRG.ising_transfer_tensor_2D(β, J, σindlist[2, 3], σindlist[3, 3], sindlist[3, 2], sindlist[3, 3])
    #縮約
    _P_edge = P_edge * W
    C_expand = P_corner * C * _P_edge
    return C_expand, cindσ, cinds
end
#右下を作成
function setup_rightdown_corner_first(β::Float64, J::Float64,
    σindlist::Matrix{Index},sindlist::Matrix{Index})
    cindσ = combiner(σindlist[2, 3], σindlist[2, 4], tags="σ_all")
    cinds = combiner(sindlist[1, 2], sindlist[2, 2], tags="s_all")
    #2脚テンソルの定義
    C = CTMRG.C_corner(β, J, σindlist[1, 4], sindlist[1, 3])
    #3脚テンソルの定義
    P_edge = CTMRG.P_edge(β, J, σindlist[2, 4], σindlist[1, 4], sindlist[2, 3])
    P_corner = CTMRG.P_corner(β, J, σindlist[1, 3], sindlist[1, 2], sindlist[1, 3])
    #4脚テンソルの定義
    W =  CTMRG.ising_transfer_tensor_2D(β, J, σindlist[1, 3], σindlist[2, 3], sindlist[2, 2], sindlist[2, 3])
    #縮約
    _P_edge = P_edge * W
    C_expand = P_corner * C * _P_edge
    return C_expand, cindσ, cinds
end
#対角化した関数を作成
function diag_coner_tensor(C_expand::ITensor,cindσ::ITensor,
    cinds::ITensor)
    _C_expand = C_expand * cindσ * cinds
    C_array = array(_C_expand)
    evals, evecs = eigen(C_array)
    Λ, A, iA = make_diag_Tesnor(evals, evecs, inds(cindσ)[1], inds(cinds)[1])
    return Λ, A, iA
end
#L+1にするテンソルの作成
function σ_edge_tensor(β::Float64, J::Float64,
    σindlist_up::Vector{Index{Int64}},σindlist_down::Vector{Index{Int64}},
    indname::String)
    Tlist = []
    ind = Index(2, "$indname")
    sdumyind_list = [Index(2, "dummy,$i") for i in 1:length(σindlist_up)-1]
    if σindlist_down[1] == σindlist_up[1]
        _σindlist_down = [prime(σindlist_down[i]) for i in 1:length(σindlist_down)]
    else
        _σindlist_down = σindlist_down
    end

    for i in 1:length(σindlist_up)
        if i == 1
            T1 = P_edge(β, J, σindlist_up[1], _σindlist_down[1], sdumyind_list[1])
            push!(Tlist, T1)
        elseif i == length(σindlist_up)
            T1 = ising_transfer_tensor_2D(β, J, σindlist_up[i], _σindlist_down[i], sdumyind_list[i-1],ind)
            push!(Tlist, T1)
        else
            T1 = ising_transfer_tensor_2D(β, J, σindlist_up[i], _σindlist_down[i], sdumyind_list[i-1], sdummyind_list[i])
            push!(Tlist, T1)
        end
    end
    result = ITensor(true)
    for T in Tlist
        result *= T
    end
    return result, ind, _σindlist_down
end