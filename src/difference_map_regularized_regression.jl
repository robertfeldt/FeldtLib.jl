# This implements the Sparse Reconstruction using a Difference Map as described in the paper:
#
#  Will Landecker, Rick Chartrand, and Simon DeDeo. "Robust Compressed 
#  Sensing and Sparse Coding with the Difference Map." Submitted to the 
#  arXiv on 10/30/2013.
#
# (x, scores) = difference_map_sparse_reconstruction(d, y, s)
#
# Finds a sparse reconstruction 'x' using the Difference Map method. 
# 'x' is the solution as in y = d * x.
#
# Outputs:
#  'x' is an array with the sparse reconstruction of y given d using the 
#      Difference Map method.
#
#  'scores' is an array indicating the accuracy of 'x' at each iteration.
#
# Inputs:
#  'd' is an N by M matrix.
#
#  'y' is an N by 1 array.
#
#  'x0' is an N by 1 array with the starting point for the iteration. Set to
#       nothing to select a random starting point.
#
#  'max_degree' is the maximum number of non-zero values (the max L0-penalty) in the solution x.
#
#  'time_limit' is max number of time (in seconds) the algorithm is allowed to run.
#
#  'delta_treshold' is the treshold value below which the iteration is deemed to have converged.
function difference_map_sparse_reconstruction(d, y, max_degree;
	beta = -0.14,
	x0 = nothing,
	time_limit = 10.0,
	d_pseudo = nothing,
	delta_treshold = 1e-10,
	verbose = true
	)

	# Check and set up parameters
	M, N = size(d)

	if max_degree > M
		error("Max degree ($(max_degree)) is larger than number of available predictors/covariates ($(M))")
	end

	if x0 == nothing
		x0 = (randn(M) .- 0.5)'
	end

	if d_pseudo == nothing
		d_pseudo = pinv(d)
	end

	# Define the functions used for the update
	num_zero = M - max_degree

	# Init temp variables
	xn = x0
	delta = Inf
	iteration = 0
	start_time = time()

	# We use the stdev in the last max_var_steps iterations as an indication
	# of stability:
	# xs = zeros(Float64, M, max_var_steps)

	# Main loop
	while( delta > delta_treshold )

		iteration += 1
		xo = xn

		# Update x
		xn = xo + beta * ( A( fB(xo, beta, num_zero), d, d_pseudo, y )  -  
			B( fA(xo, beta, d, d_pseudo, y), num_zero ) )

		# Change in x in this step
		delta = sum((xo - xn) .^ 2)

		# How well do we reconstruct?
		reconstruction_distance = residual_sum_of_squares(y, B(xn, num_zero) * d)

		# Print some info
		if verbose
      println(@sprintf("DM. iter %s\n  reconstruction distance = %0.5e\n  L2 norm of update distance = %0.5e", iteration, reconstruction_distance, delta))
      show(xn)
      println()
    end

    # Break if over time limit
    delta_t = time() - start_time
    if delta_t > time_limit
    	break
    end

	end

	if verbose
		println("Finished after $(iteration) iterations.")
	end

	# Ensure the solution is sparse.
	x = B(xn, num_zero)

	return( x, iteration, aic(d, y, x) )
end

A(x, d, d_pseudo, y) = x .- ( (x*d .- y) * d_pseudo )

function fA(x, beta, d, d_pseudo, y)
	Ax = A(x, d, d_pseudo, y)
	Ax - (Ax - x)/beta
end

function B(x, num_zero)
	# Set the num_zero smallest values to 0
	res = copy(x)
	res[sortperm(x[:])[1:num_zero]] = 0.0
	res
end

function fB(x, beta, num_zero)
	Bx = B(x, num_zero)
	Bx + (Bx - x)/beta
end

residual_sum_of_squares(y, yhat) = sum( (yhat .- y) .^ 2 )

# AIC assuming independent errors and constant variance:
#   AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
function aic_from_rss(rss, d, coefs)
	k = sum(abs(coefs) .> 0)
	m, n = size(d)
	return( 2*k + n * (log(2*pi*rss/n) + 1) )
end

aic(d, y, coefs) = aic_from_rss( residual_sum_of_squares(y, coefs * d), d, coefs)

# Return the reconstruction that gives the lowest AIC assuming independent errors and
# constant variance.
function min_aic_difference_map_reconstruction(d, y;
	time_limit = 2.0,
	max_degree = nothing,
	min_degree = 1
	)

	M, N = size(d)

	if max_degree == nothing
		max_degree = int( 2*sqrt(M) )
	end

	# Pre-calc so we need not do it each time.
	dps = pinv(d)

	degrees = min_degree:max_degree
	num_degrees = Base.length(degrees)
	xs = zeros(Float64, M, num_degrees)
	its = zeros(Int64, num_degrees)
	aics = zeros(Float64, num_degrees)

	time_per_run = time_limit / num_degrees

	for i in 1:num_degrees
		deg = max_degree + 1 - i
		x, its[i], aics[i] = difference_map_sparse_reconstruction(d, y, deg; 
			time_limit = time_per_run, d_pseudo = dps, verbose = false)
		xs[:,i] = x'
		println("Min AIC DM, degree = $(deg), aic = $(aics[i]), its = $(its[i])")
	end

	p = sortperm(aics)
	bestindex = p[1]
	bestmodel = xs[:,bestindex]
	degs = reverse(degrees)[p]

	return( bestmodel, degs[1], xs[:,p], its[p], aics[p], degs )

end

# Zero center and unit variance of each row of a matrix.
function zero_center_and_normalize(d)
	d0 = d .- mean(d, 2)
	d0 ./ std(d, 2)
end

# Use stability selection (in the spirit of Meinshausen & Buhlman) with the 
# DM sparse reconstruction to find the covariates that are most consistently 
# in the selected set.
# 
#  inclusion_probability = stability_selection_dm(d, y)
#
#  
# Outputs:
#  'inclusion_probability' is the probability that a covariate/feature is
#                          selected when half the samples are used and the features
#                          reweighted by a random weight uniformly sampled in [alpha, 1].
#
function stability_selection_dm(d, y;
	num_bootstrap = 100,
	alpha = 0.20,
	min_degree = nothing,
	time_limit = 1.0,
	verbose = true
	)

	p, n = size(d)

	if min_degree == nothing
		min_degree = int(p/2)
	end

	#dn = zero_center_and_normalize(d)

	halfsize = int(n/2)
	counts = zeros(Float64, 1, p)	

	time_per_dm = time_limit / (2*num_bootstrap)

	for i in 1:num_bootstrap

		# Randomly reweight the covariate/predictor values.
		dr = dn .* (alpha .+ (1.0 - alpha) .* rand(p, n))

		# Randomly split the sample in two sets.
		permutation = shuffle(collect(1:n))
		set1 = permutation[1:halfsize]
		set2 = permutation[(halfsize+1):n]

		# Run the difference map sparse reconstructions on each set.
		degree = rand(min_degree:(p-1))
		dm1, its1 = difference_map_sparse_reconstruction(dr[:,set1], y[:,set1], degree; time_limit = time_per_dm, verbose = false)
		dm2, its2 = difference_map_sparse_reconstruction(dr[:,set2], y[:,set2], degree; time_limit = time_per_dm, verbose = false)

		# Update the counts for the selected features
		counts += abs(sign(dm1))
		counts += abs(sign(dm2))

		if verbose
			println("$(i):")
			#println("  dm1 = $(dm1)")
			#println("  dm2 = $(dm2)")
			println("  counts = $(counts/(2*i))")
		end
	end
			
	# Normalize frequency to a probability
	probs = counts / (2*num_bootstrap)

	# the final stability score is the maximum frequency over the steps
	# result <- apply(freq,2,max)
	return probs
end
