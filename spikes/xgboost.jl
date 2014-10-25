using XGBoost

D = 10
N = 200
NT = 100

x = float32(randn(N, D))
function rastrigin(x)
  D = length(x)
  10 * D + sum( x.^2 ) - 10 * sum( cos( 2 * π * x ) )
end
y = float32(vec(mapslices(rastrigin, x, 2)))
xx = float32(randn(NT, D))
yy = float32(vec(mapslices(rastrigin, xx, 2)))

bt = xgboost(x, 30, label = y, eta=1, max_depth=4, objective="reg:linear")
ypred = predict(bt, xx)
medape(ypred, y) = 100.0 * median(abs((ypred .- y) ./ y))
medape(ypred, yy)
