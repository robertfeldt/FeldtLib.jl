	beta = -0.14
	x0 = nothing
	time_limit = 10.0
	d_pseudo = nothing
	delta_treshold = 1e-10
	verbose = true

max_degree = 5
num_bootstrap = 100
alpha = 0.20


# Get a regression problem with 2 covars, expand it and use a model with 1 active covar.
rp, y, errors = random_regression_problem(2, 1)
# Use DM to reconstruct the coefficients:
m, k, xs, its, aics, degs = min_aic_difference_map_reconstruction(rp.d.data, y)
analyze_regression(rp, m)

# Test on random data
N = 30
M = 2
d = randn(M, N)


y = y1 = 4.5*d[1,:] + 17.0*d[2,:] + 8.2*d[3,:]
y2 = 1.2*d[1,:] + 2.0*d[2,:]
y3 = 1.5*d[1,:]

xbest1, k1, xs1, its1, aics1, degs1 = min_aic_difference_map_reconstruction(d, y1)
xbest2, k2, xs2, its2, aics2, degs2 = min_aic_difference_map_reconstruction(d, y2)
xbest3, k3, xs3, its3, aics3, degs3 = min_aic_difference_map_reconstruction(d, y3)


x4, i4, aic4 = difference_map_sparse_reconstruction(d, y, 4; time_limit = 1.0)
x3, i3, aic3 = difference_map_sparse_reconstruction(d, y, 3; time_limit = 1.0)

# The pattern seems to be that when we come down below the actual
# degree the reconstruction does not converge:
x2, i2, aic2 = difference_map_sparse_reconstruction(d, y, 2; time_limit = 1.0)
x1, i1, aic1 = difference_map_sparse_reconstruction(d, y, 1; time_limit = 1.0)

# There also seems to be a pattern that we get faster convergence if starting
# from an earlier fix point:
x4b, i4b = difference_map_sparse_reconstruction(d, y, 3; x0 = x5)
x3b, i3b = difference_map_sparse_reconstruction(d, y, 3; x0 = x4)
x2b, i2b = difference_map_sparse_reconstruction(d, y, 2; x0 = x3, time_limit = 1.0)
x1b, i1b = difference_map_sparse_reconstruction(d, y, 1; x0 = x2, time_limit = 1.0)

xr5 = repeated_difference_map_sparse_reconstruction(d, y, 5; time_limit = 1.0)
xr4 = repeated_difference_map_sparse_reconstruction(d, y, 4; time_limit = 1.0)
xr3 = repeated_difference_map_sparse_reconstruction(d, y, 3; time_limit = 1.0)

