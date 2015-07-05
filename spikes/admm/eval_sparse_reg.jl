mape(y, yhat) = mean(100.0 * abs((y .- yhat) ./ y))
medape(y, yhat) = median(100.0 * abs((y .- yhat) ./ y))

function eval_sparse_reg(A, b, xtrue, xhat)
    bhat = A * xhat
    strue = Set(find(xtrue))
    shat = Set(find(xhat))
    missing = setdiff(strue, shat)
    common = setdiff(strue, missing)
    additional = setdiff(shat, strue)
    mapeval = mape(b, bhat)
    medapeval = medape(b, bhat)
    println("MAPE: $(round(mapeval,1))%, MEDAPE: $(round(medapeval,1))%\n",
        "  Support: common = $(length(common)), missing = $(length(missing)), additional = $(length(additional))")
    return [mapeval, medapeval, length(common), length(missing), length(additional)]
end
