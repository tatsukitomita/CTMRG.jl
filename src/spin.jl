function spin_tensor(ind1::Index, ind2::Index)
    S = ITensor(ind1, ind2)
    S[ind1[1], ind2[1]] = 1.0
    S[ind1[2], ind2[2]] = -1.0
    return S
end