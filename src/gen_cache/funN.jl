function funN_rep_gen(N::Int)
    x1 = zeros(N, N)
    x2 = zeros(N, N)
    for i in 1:N
        x1[i,mod(i,N)+1] = 1
    end
    x2[1,2] = x2[2,1] = 1
    for i in 3:N
        x2[i,i] = 1
    end
    return [x1, x2]
end
function FunN(R::Type{SNIrrep{N}}) where {N}
    if R == S3Irrep
        irrep_gen = irreps_gen.S3
    elseif R == S4Irrep
        irrep_gen = irreps_gen.S4
    elseif R == S5Irrep
        irrep_gen = irreps_gen.S5
    else
        error("$R is not implemented.")
    end

    space = ()
    embedders = Dict()
    for (c, s) in enumerate(values(R))
        B = get_intertwiners(irrep_gen[c], funN_rep_gen(N))
        space = (space..., s=>length(B))
        embedders[s] = B
    end

    A = Vect[R](space...)


    mul_nosect = zeros(N,N,N)
    for i in 1:N
        mul_nosect[i, i, i] = 1
    end

    m = zeros(A←A⊗A)
    for i in sectors(A), j in sectors(A), k in sectors(A) # i←j⊗k
        for deg in 1:Nsymbol(j, k, i)
            for x in 1:dim(A, i), y in 1:dim(A, j), z in 1:dim(A, k) # x←y⊗z
                @tensor intertwiner_expansion[up; down] := mul_nosect[a; b c] * embedders[i][x]'[up;a] * embedders[j][y][b; left] * embedders[k][z][c; right] * CGC(j, k, i)[:,:,:,deg][left right; down]
                numb = intertwiner_expansion[1,1]
                m[FusionTree{R}((i,), i, (false,), ()), FusionTree{R}((j, k), i, (false, false), ())] .= reshape([numb], (dim(A, i), dim(A, j), dim(A, k))) # We used the fact that in this case, the fusion channel happens to be unique
            end
        end
    end


    ι = zeros(A)

    ι[FusionTree{R}((one(R),), one(R), (false,), ()), FusionTree{R}((), one(R), (), (), ())] .= reshape([sqrt(N)], (dim(A, one(R)), ))

    return A, m, ι
end


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
