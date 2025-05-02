A, m, ι = FunN(S5Irrep)

@planar AA_A[a123; a1 a2 a3] := m[a123; a12 a3] * m[a12; a1 a2]
@planar A_AA[a123; a1 a2 a3] := m[a123; a1 a23] * m[a23; a2 a3]

@assert isapprox(AA_A, A_AA, atol=1e-14) # check associativity

@planar ιA[a; b] := m[a; x b] * ι[x]

@assert isapprox(ιA, id(A), atol=1e-14) # check left unitality

@planar Aι[a; b] := m[a; b x] * ι[x]

@assert isapprox(Aι, id(A), atol=1e-14) # check right unitality

@planar mm[a; b] := m[a; x y] * m'[x y; b]

@assert isapprox(mm, id(A), atol=1e-14) # check unitarity