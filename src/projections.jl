function _proj_custom(x, cone)
    error("Not Implemented")
    # TODO: Exponential Cone:
    if cone <: MOI.ExponentialCone
        # X = (r, s, t)
        # Projection: https://web.stanford.edu/~boyd/papers/pdf/prox_algs.pdf 6.3
        # if x is in exponential cone: return x
        # elif -x is in exponential cone: return 0
        # elif first two coordinates are less than zero, return (r, 0, max(t, 0))
        # otherwise, solve constrained optimization
        continue
    elseif typeof(cone) <: MOI.PowerCone
        # TODO: Power cone
        # https://link.springer.com/article/10.1007/s00186-015-0514-0
        # or http://www.optimization-online.org/DB_FILE/2014/08/4502.pdf (Proposition 2.3)
        # if x is in power cone: return x
        # elif -x is in power cone: return -x
        # elif z == 0: return (max(x, 0), max(y, 0), 0)
        # otherwise, solve constrained optimization
    end
end


function _d_proj_custom(x, cone)
    error("Not Implemented")
    if cone <: MOI.ExponentialCone
        # TODO: Exponential Cone:
        # Grad: http://proceedings.mlr.press/v70/ali17a/ali17a.pdf lemma 3.6
        # Supplement: https://arxiv.org/pdf/1705.00772.pdf
        continue
    elseif cone <: MOI.PowerCone
        # TODO: Power cone
        # x = (x, y, z) ∈ R × R × R^n
        # https://link.springer.com/article/10.1007/s00186-015-0514-0
        # or http://www.optimization-online.org/DB_FILE/2014/08/4502.pdf (Theorem 3.1)
        continue
    end
end
