function prs_predict(x::Array{Float64, 2}, y::Array{Float64, 2}, 
  xx::Array{Float64, 2}; )

# Some radial basis functions
rbf_thin_plate(r::Float64, c::Float64 = 1.0) = (r == 0.0) ? 0.0 : (r^2 * log(c*r))

rbf_gaussian(r::Float64, c::Float64 = 1.0) = exp(-c*r^2)

rbf_cubic(r::Float64, c::Float64 = 1.0) = r^3

rbf_multiquadric(r::Float64, c::Float64 = 1.0) = sqrt(r^2 + c^2)

# Calculate the pairwise kernel matrix for the points in x taken to be the
# individual columns.
function calc_kernel_matrix(x::Array{Float64, 2}, centers::Array{Float64, 2}, 
  phifunc = rbf_gaussian, distfunc = norm, c = 1.0)
  num_points = size(x, 2)
  num_centers = size(centers, 2)
  k_matrix = zeros(Float64, num_points, num_centers)
  for i in 1:num_points
    xi = x[:,i]
    for j in 1:num_centers
      k_matrix[i, j] = phifunc(norm(xi - centers[:,j]), c)
    end
  end
  k_matrix
end

function rbf_predict(x::Array{Float64, 2}, y::Array{Float64, 2}, 
  xx::Array{Float64, 2}; c = nothing, rbf = rbf_gaussian)

  # Set c to ok values unless it was specified
  if c == nothing && rbf == rbf_gaussian
    D = size(x, 1)
    c = 1/(D*D)
  elseif c == nothing
    c = 1.0
  end

  # Scale x per row (i.e. covariate)
  x_row_means = mean(x, 2)
  x_row_std = std(x, 2)
  xt = broadcast(-, x, x_row_means)
  xsc = broadcast(/, xt, x_row_std)

  # Scale y per its mean and std
  ymean = mean(y)
  ystd = std(y)
  ysc = (y - ymean) / ystd

  # Dor LSR to get lambdas
  A = calc_kernel_matrix(xsc, xsc, rbf, norm, c)
  lambdas = A\(ysc')

  # Scale xx
  xxt = broadcast(-, xx, x_row_means)
  xxsc = broadcast(/, xxt, x_row_std)

  # Predict for the test points
  B = calc_kernel_matrix(xxsc, xsc, rbf, norm, c)
  ypredsc = (B * lambdas)'

  # Transform back to original range
  ypred = ymean + ystd * ypredsc

  return ypred

end

function r2(ypred, y)
  SSR = sum( (y .- ypred).^2  )
  SST = sum( (y - mean(y)).^2 )
  1.0 - SSR / SST
end

# Test function
tf(x) = sum(x.^4)

N = 100
D = 10
xtraining = randn(D, N)
ytraining = mapslices(tf, xtraining, 1)
xtest = randn(D, 5*D)
y = mapslices(tf, xtest, 1)

ypred_tp = rbf_predict(xtraining, ytraining, xtest; rbf = rbf_thin_plate)
r2_tp = r2(ypred_tp, y)

ypred_multiquadric = rbf_predict(xtraining, ytraining, xtest; rbf = rbf_multiquadric)
r2_multiquadric = r2(ypred_multiquadric, y)

ypred_cubic = rbf_predict(xtraining, ytraining, xtest; rbf = rbf_cubic)
r2_cubic = r2(ypred_cubic, y)

ypred_gaussiand = rbf_predict(xtraining, ytraining, xtest; rbf = rbf_gaussian, c = 1/D)
r2_gaussiand = r2(ypred_gaussiand, y)

ypred_gaussiand2 = rbf_predict(xtraining, ytraining, xtest; rbf = rbf_gaussian, c = 1/(D*D))
r2_gaussiand2 = r2(ypred_gaussiand2, y)