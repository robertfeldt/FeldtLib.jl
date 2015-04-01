
# R2, the coefficient of determination. The higher/closer to 1.0 it is the better
# the model is, i.e. the closer the predicted y values are to the actual y values.
function r2(ypred, y)
  SSR = sum( (y .- ypred).^2  )
  SST = sum( (y - mean(y)).^2 )
  1.0 - SSR / SST
end

mape(ypred, y) = 100.0 * mean(abs((ypred .- y) ./ y))
medape(ypred, y) = 100.0 * median(abs((ypred .- y) ./ y))

using XGBoost

function rank_rbfs(x::Array{Float64, 2}, y::Array{Float64, 2}, 
  xx::Array{Float64, 2}, yy::Array{Float64, 2}; rbfs = {
  "ThinPlate(1.0)"      => (rbf_thin_plate, 1.0),
  "ThinPlate(10.0)"      => (rbf_thin_plate, 10.0),
  "Cubic(1.0)"          => (rbf_cubic, 1.0),
  "Multiquadric(1.0)"   => (rbf_multiquadric, 1.0),
  "Gaussian(1/D^2)"     => (rbf_gaussian, 1/(size(x,1)^2))})

  rbfs = collect(rbfs)
  num_rbfs = length(rbfs)
  r2s = zeros(Float64, num_rbfs+1)
  medapes = zeros(Float64, num_rbfs+1)

  for i in 1:num_rbfs
    rbf, c = rbfs[i][2]
    ypred, rbfm = rbf_predict(x, y, xx; rbf = rbf, c = c)
    r2s[i] = round(r2(ypred, yy), 3)
    medapes[i] = round(medape(ypred, yy), 2)
  end

  # Compare to XGBoost
  xgm = xgboost(float32(x)', 2, label = vec(float32(y)), eta=1, max_depth=2, objective="reg:linear")
  ypred = XGBoost.predict(xgm, float32(xx)')
  r2s[num_rbfs+1] = round(r2(ypred, yy), 3)
  medapes[num_rbfs+1] = round(medape(ypred, yy), 2)

  # Compare to Random forest
  rfm = build_forest(y, x, 3, 1000, 5, 0.7)
  ypred2 = apply_forest(model, xx)
  r2s[num_rbfs+2] = round(r2(ypred2, yy), 3)
  medapes[num_rbfs+2] = round(medape(ypred2, yy), 2)

  perm = sortperm(-medapes; rev = true)

  for i in 1:num_rbfs
    j = perm[i]
    name = if j == num_rbs + 1
      "xgboost"
    elseif j == num_rbs + 2
      "randforest"
    else
      string(rbfs[j][1])
    end
    println("$(i). $(name), R2 = $(r2s[j]), MEDAPE = $(medapes[j])")
  end

end

# Test functions
function rastrigin(x)
  D = length(x)
  10 * D + sum( x.^2 ) - 10 * sum( cos( 2 * π * x ) )
end

test_functions = {
  "Sum of equally weighted quadratics" => (x) -> sum(x.^2),
  "Sum of equally weighted cubics" => (x) -> sum(x.^3),
  "Sum of non-equally weighted cubics" => (x) -> 10*x[1]^3 + sum(x[2:end].^3),
  "Sum of non-equally weighted quadratics" => (x) -> 10*x[1]^2 + sum(x[2:end].^2),
  "Rastrigin" => rastrigin
}
tfs = collect(test_functions)

D = 10
T = 50

Ns = [100, 500]

for i in 1:length(tfs)
  tfdesc, tf = tfs[i]

  for N in Ns
    xtraining = randn(D, N)
    ytraining = mapslices(tf, xtraining, 1)
    xtest = randn(D, T)
    ytest = mapslices(tf, xtest, 1)

    println("For test function $(tfdesc), N = $(N)")
    rank_rbfs(xtraining, ytraining, xtest, ytest)
    println("")
  end

end