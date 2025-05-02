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



