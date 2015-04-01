using DecisionTree

D = 10
N = 500
NT = 100

x = randn(N, D)
function rastrigin(x)
  D = length(x)
  10 * D + sum( x.^2 ) - 10 * sum( cos( 2 * π * x ) )
end
y = vec(mapslices(rastrigin, x, 2))
xx = randn(NT, D)
yy = vec(mapslices(rastrigin, xx, 2))

model = build_forest(y, x, 3, 1000, 5, 0.7)
ypred2 = apply_forest(model, xx)

medape(ypred, y) = 100.0 * median(abs((ypred .- y) ./ y))
medape(ypred2, yy)
