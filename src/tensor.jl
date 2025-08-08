function freeend(ind1::Index)
    D = ITensor(ind1)
    D[ind1 => 1] = 1.0
    D[ind1 => 2] = 1.0
    return D
end
function fixend(ind1::Index)
    D = ITensor(ind1)
    D[ind1 => 1] = 1.0
    D[ind1 => 2] = 0.0
    return D
end 
#=
function oneleg_tensor(idx::Index)
    T = ITensor(idx)
    return T,idx
end

function oneleg_tensor(s1::String,s1_value::Int, m::String,n::Int)
    sites_s = Index(s1_value, "$s1,$m=$n")
    T = ITensor(sites_s)
    return T,sites_s
end

function twoleg_tensor(idx1::Index, idx2::Index)
    T = ITensor(idx1, idx2)
    return T, idx1, idx2
end

function twoleg_tensor(s1::String, s2::String,
    s1_value::Int, s2_value::Int, 
    m::String, n1::Int,n2::Int)
    sites_s1 = Index(s1_value, "$s1,$m=$n1")
    sites_s2 = Index(s2_value, "$s2,$m=$n2")
    T = ITensor(sites_s1, sites_s2)
    return T, sites_s1, sites_s2
end

function threeleq_tensor(idx1::Index, idx2::Index, idx3::Index)
    T = ITensor(idx1, idx2, idx3)
    return T, idx1, idx2, idx3
end
function fourleg_tensor(idx1::Index, idx2::Index, idx3::Index, idx4::Index)
    T = ITensor(idx1, idx2, idx3, idx4)
    return T, idx1, idx2, idx3, idx4
end
=#


