function get_random_product_cone(n)
    zero_dim = rand(1:n)
    nonneg_dim = rand(1:n)
    soc_dim = rand(1:n)
    psd_dim = rand(1:n)
    cone_prod = [MOI.Zeros(zero_dim), MOI.Nonnegatives(nonneg_dim),
             MOI.SecondOrderCone(soc_dim), MOI.PositiveSemidefiniteConeTriangle(psd_dim)]
    return cone_prod
end
