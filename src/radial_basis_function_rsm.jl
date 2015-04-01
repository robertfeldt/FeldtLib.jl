# Some radial basis functions
rbf_thin_plate{T <: FloatingPoint}(r::T, c::T = 1.0) = (r == 0.0) ? 0.0 : (r^2 * log(c*r))

rbf_gaussian{T <: FloatingPoint}(r::T, c::T = 1.0) = exp(-c*r^2)

rbf_cubic{T <: FloatingPoint}(r::T, c::T = 1.0) = r^3

rbf_multiquadric{T <: FloatingPoint}(r::T, c::T = 1.0) = sqrt(r^2 + c^2)

identity_function{T <: FloatingPoint}(r::T) = r

euclidean_distance{T <: FloatingPoint}(x1::Array{T, 1}, x2::Array{T, 1}) = norm(x1 - x2)

# Create a kernel function from a RBF, a distance function and a c constant that
# is passed on to the RBF.
function create_kernel_func(rbffunc, distancefunc = euclidean_distance, c = 1.0)
  (x1, x2) -> rbffunc(distancefunc(x1, x2), c)
end

# Calculate the pairwise kernel matrix for the points in x taken to be the
# individual columns. The kernel matrix can be padded with zeros by specifying
# num_padding_columns and num_padding_rows
function calc_kernel_matrix{T <: FloatingPoint}(x::Array{T, 2}, centers::Array{T, 2},
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
function scale_and_center{T <: FloatingPoint}(x::Array{T, 2}, center = true, scale = true)

  if center == true
    center = mean(x, 2)
    xres = broadcast(-, x, center)
  elseif typeof(center) == Array{T, 2}
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
  elseif typeof(scale) == Array{T, 2}
    if size(scale, 1) != size(x, 1)
      throw("When giving a vector of scale values they must have as many rows as the number of rows per point; here they differ!")
    end
    xres = broadcast(/, xres, scale)
  end

  return xres, center, scale

end

type RadialBasisFunctionModel{T <: FloatingPoint}
  kernel::Function

  xcenters::Array{T, 2}     # X values for the centers. One center per column.
  x_row_means::Array{T, 2}  # Average value for the original, uncentered and unscaled xcenters.
  x_row_stds::Array{T, 2}   # Standard dev for the original, uncentered and unscaled xcenters.

  fcenters::Array{T, 2}     # Y values for the centers. One value per as many columns as in xcenters.
  f_mean::T                 # Mean value for the original, uncentered and unscaled function values.
  f_std::T                  # Standard dev for the original, uncentered and unscaled function values.

  lambdas::Array{T, 2}      # Coefficients in the RBF model.
end

function RBFModel{T <: FloatingPoint}(x::Array{T,2}, y::Array{T,2}; 
  rbf = rbf_gaussian, 
  c = nothing,
  distance_func = euclidean_distance)

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

  RadialBasisFunctionModel(kernel, xcenters, x_row_means, x_row_stds, fcenters, f_mean[1], f_std[1], lambdas)

end


function predict{T <: FloatingPoint}(rbfm::RadialBasisFunctionModel, xx::Array{T, 2})
  # Scale xx
  xx_scaled, c, s = scale_and_center(xx, rbfm.x_row_means, rbfm.x_row_stds)

  # Predict for the test points in xx
  B = calc_kernel_matrix(xx_scaled, rbfm.xcenters, rbfm.kernel)
  ypred_scaled = (B * rbfm.lambdas)'

  # Transform back to original range
  ypred = rbfm.f_mean + rbfm.f_std * ypred_scaled

  return ypred
end

# Convenience function that builds the model from x and y and then uses the model
# to predict the y values for the points in xx.
function rbf_predict{T <: FloatingPoint}(x::Array{T, 2}, y::Array{T, 2}, 
  xx::Array{T, 2}; options...)

  rbfm = RBFModel(x, y; options...)
  ypred = predict(rbfm, xx)
  return ypred, rbfm

end
