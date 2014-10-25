# Some radial basis functions
rbf_thin_plate(r::Float64, c::Float64 = 1.0) = (r == 0.0) ? 0.0 : (r^2 * log(c*r))

rbf_gaussian(r::Float64, c::Float64 = 1.0) = exp(-c*r^2)

rbf_cubic(r::Float64, c::Float64 = 1.0) = r^3

rbf_multiquadric(r::Float64, c::Float64 = 1.0) = sqrt(r^2 + c^2)

identity_function(r::Float64) = r

euclidean_distance(x1::Array{Float64, 1}, x2::Array{Float64, 1}) = norm(x1 - x2)

# Create a kernel function from a RBF, a distance function and a c constant that
# is passed on to the RBF.
function create_kernel_func(rbffunc, distancefunc = euclidean_distance, c = 1.0)
  (x1::Array{Float64, 1}, x2::Array{Float64, 1}) -> rbffunc(distancefunc(x1, x2), c)
end

# Calculate the pairwise kernel matrix for the points in x taken to be the
# individual columns. The kernel matrix can be padded with zeros by specifying
# num_padding_columns and num_padding_rows
function calc_kernel_matrix(x::Array{Float64, 2}, centers::Array{Float64, 2},
  kernelfunc = create_kernel_func(rbf_gaussian, euclidean_distance, 1/size(x, 1)))

  num_points = size(x, 2)
  num_centers = size(centers, 2)
  k_matrix = zeros(Float64, num_points, num_centers)

  for i in 1:num_points
    xi = x[:,i]
    for j in 1:num_centers
      k_matrix[i, j] = kernelfunc(xi, centers[:,j])
    end
  end

  k_matrix

end

# Scale and center with the mean and standard deviation (or with specific 
# center and value vectors given as input (if given)).
function scale_and_center(x::Array{Float64, 2}, center = true, scale = true)

  if center == true
    center = mean(x, 2)
    xres = broadcast(-, x, center)
  elseif typeof(center) == Array{Float64, 2}
    if size(center, 1) != size(x, 1)
      throw("When giving a vector of center values they must have as many rows as the number of rows per point; here they differ!")
    end
    xres = broadcast(-, x, center)
  else
    xres = x
  end

  if scale == true
    scale = std(x, 2)
    xres = broadcast(/, xres, scale)
  elseif typeof(scale) == Array{Float64, 2}
    if size(scale, 1) != size(x, 1)
      throw("When giving a vector of scale values they must have as many rows as the number of rows per point; here they differ!")
    end
    xres = broadcast(/, xres, scale)
  end

  return xres, center, scale

end

type RadialBasisFunctionModel
  kernel::Function

  xcenters::Array{Float64, 2}     # X values for the centers. One center per column.
  x_row_means::Array{Float64, 2}  # Average value for the original, uncentered and unscaled xcenters.
  x_row_stds::Array{Float64, 2}   # Standard dev for the original, uncentered and unscaled xcenters.

  fcenters::Array{Float64, 2}     # Y values for the centers. One value per as many columns as in xcenters.
  f_mean::Float64                 # Mean value for the original, uncentered and unscaled function values.
  f_std::Float64                  # Standard dev for the original, uncentered and unscaled function values.

  lambdas::Array{Float64, 2}      # Coefficients in the RBF model.

  RadialBasisFunctionModel(x::Array{Float64, 2}, y::Array{Float64, 2}; 
    rbf = rbf_gaussian, c = nothing,
    distance_func = euclidean_distance) = begin

    # Set c to ok values unless it was specified
    if c == nothing && rbf == rbf_gaussian
      num_dims = size(x, 1)
      c = 1/(num_dims * num_dims)
    elseif c == nothing
      c = 1.0
    end

    kernel = create_kernel_func(rbf, distance_func, c)

    # Scale and center x and y (per row, i.e. covariate)
    xcenters, x_row_means, x_row_stds = scale_and_center(x)
    fcenters, f_mean, f_std = scale_and_center(y)

    # Create kernel matrix.
    A = calc_kernel_matrix(xcenters, xcenters, kernel)

    # Solve for the lambda coefficients.
    lambdas = A\(fcenters')

    new(kernel, xcenters, x_row_means, x_row_stds, fcenters, f_mean[1], f_std[1], lambdas)

  end
end

function predict(rbfm::RadialBasisFunctionModel, xx::Array{Float64, 2})

  # Scale xx
  xx_scaled, c, s = scale_and_center(xx, rbfm.x_row_means, rbfm.x_row_stds)

  # Predict for the test points in xx
  B = calc_kernel_matrix(xx_scaled, rbfm.xcenters, rbfm.kernel)
  ypred_scaled = (B * rbfm.lambdas)'

  # Transform back to original range
  ypred = rbfm.f_mean + rbfm.f_std * ypred_scaled

  return ypred

end

# R2, the coefficient of determination. The higher/closer to 1.0 it is the better
# the model is, i.e. the closer the predicted y values are to the actual y values.
function r2(ypred, y)
  SSR = sum( (y .- ypred).^2  )
  SST = sum( (y - mean(y)).^2 )
  1.0 - SSR / SST
end

mape(ypred, y) = 100.0 * mean(abs((ypred .- y) ./ y))
medape(ypred, y) = 100.0 * median(abs((ypred .- y) ./ y))

# Convenience function that builds the model from x and y and then uses the model
# to predict the y values for the points in xx.
function rbf_predict(x::Array{Float64, 2}, y::Array{Float64, 2}, 
  xx::Array{Float64, 2}; 
  c = nothing, rbf = rbf_gaussian, distance_func = euclidean_distance)

  rbfm = RadialBasisFunctionModel(x, y; 
    rbf = rbf, c = c, distance_func = distance_func)

  ypred = predict(rbfm, xx)

  return ypred, rbfm

end

function rank_rbfs(x::Array{Float64, 2}, y::Array{Float64, 2}, 
  xx::Array{Float64, 2}, yy::Array{Float64, 2}; rbfs = {
  "ThinPlate(1.0)"      => (rbf_thin_plate, 1.0),
  "ThinPlate(10.0)"      => (rbf_thin_plate, 10.0),
  "Cubic(1.0)"          => (rbf_cubic, 1.0),
  "Multiquadric(1.0)"   => (rbf_multiquadric, 1.0),
  "Gaussian(1/D^2)"     => (rbf_gaussian, 1/(size(x,1)^2))})

  rbfs = collect(rbfs)
  num_rbfs = length(rbfs)
  r2s = zeros(Float64, num_rbfs)
  medapes = zeros(Float64, num_rbfs)

  for i in 1:num_rbfs
    rbf, c = rbfs[i][2]
    ypred, rbfm = rbf_predict(x, y, xx; rbf = rbf, c = c)
    r2s[i] = round(r2(ypred, yy), 3)
    medapes[i] = round(medape(ypred, yy), 2)
  end

  perm = sortperm(-medapes; rev = true)

  for i in 1:num_rbfs
    j = perm[i]
    println("$(i). $(rbfs[j][1]), R2 = $(r2s[j]), MEDAPE = $(medapes[j])")
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