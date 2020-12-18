using ConeProgramDiff
using MathOptInterface
const MOI = MathOptInterface

function get_random_product_cone(n)
    zero_dim = rand(1:n)
    nonneg_dim = rand(1:n)
    soc_dim = rand(1:n)
    psd_dim = rand(1:n)
    cone_prod = [MOI.Zeros(zero_dim), MOI.Nonnegatives(nonneg_dim),
             MOI.SecondOrderCone(soc_dim), MOI.PositiveSemidefiniteConeTriangle(psd_dim)]
    return cone_prod
end


function check_col_major_to_row_major()
    a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    ix_map = ConeProgramDiff.index_map_to_row_major(4)
    b = a[ix_map]
    @assert all(b .== [1, 2, 4, 7, 3, 5, 8, 6, 9, 10])
end


function check_row_major_to_col_major()
    a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    ix_map = ConeProgramDiff.index_map_to_col_major(4)
    b = a[ix_map]
    @assert all(b .== [1, 2, 5, 3, 6, 8, 4, 7, 9, 10])
end

function check_rewrite_sdp_row_major()
    A_orig = rand(18, 5)
    b_orig = rand(18)
    A = copy(A_orig)
    b = copy(b_orig)
    cones = [
        MOI.PositiveSemidefiniteConeTriangle(3),
        MOI.Zeros(2),
        MOI.PositiveSemidefiniteConeTriangle(4)
    ]
    ConeProgramDiff.rewrite_sdps_to_row_major!(A, b, cones)
    @assert all(A[[1,2,3,4,5,6], :] .== A_orig[[1,2,4,3,5,6], :])
    @assert all(A[[7,8], :] .== A_orig[[7,8], :])
    @assert all(A[[9,10,11,12,13,14,15,16,17,18], :] .== A_orig[[9,10,12,15,11,13,16,14,17,18], :])
end


function check_rewrite_sdp_col_major()
    A_orig = rand(18, 5)
    b_orig = rand(18)
    A = copy(A_orig)
    b = copy(b_orig)
    cones = [
        MOI.PositiveSemidefiniteConeTriangle(3),
        MOI.Zeros(2),
        MOI.PositiveSemidefiniteConeTriangle(4)
    ]
    ConeProgramDiff.rewrite_sdps_to_col_major!(A, b, cones)
    @assert all(A[[1,2,3,4,5,6], :] .== A_orig[[1,2,4,3,5,6], :])
    @assert all(A[[7,8], :] .== A_orig[[7,8], :])
    @assert all(A[[9,10,11,12,13,14,15,16,17,18], :] .== A_orig[[9,10,13,11,14,16,12,15,17,18], :])
end

check_col_major_to_row_major()
check_row_major_to_col_major()
