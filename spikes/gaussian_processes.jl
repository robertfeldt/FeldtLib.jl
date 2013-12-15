# Gaussian Processes in Julia
# based on http://mrmartin.net/?p=223

#####################################################################
# 1. Univariate Normal
#####################################################################

# in order to generate a univariate random number, use any real numbers mu
# and sigma to define the distribution.
mu = 6.2          # mean
sigma_squared = 2 # variance

# randn returns a random number from the standard normal distribution, which
# is a univariate normal with mu = 0, sigma = 1, i.e. N(0, 1):
standard_normal_random_number = randn()

# which we can then convert to arbitrary univariate gaussian normal, N(mu, sigma):
normal_random_number = mu + sqrt(sigma_squared) * standard_normal_random_number

# We can generate many of them:
normal_random_numbers = mu + sqrt(sigma_squared) * randn(100000)

# and plot what is going on:
using Gadfly
using DataFrames
df = DataFrame(nums = normal_random_numbers);
draw(PNG("hist_normal_rand_nums.png", 50cm, 36cm),
  plot(df, x = "nums", Geom.histogram(bincount=50)));


#####################################################################
# 2. Multivariate Normal
#####################################################################

# The same trick works for the multivariate normal. Mu has to be a vector,
# and Sigma a symetric semidefinite matrix
Mu = [3, 0]

# a matrix is symetric iff M(a,b) = M(b,a), a positive semidefinite iff
# x*M*x'>=0 for any vector x.
Sigma =  [ 1.0000   -0.9195;  -0.9195    1.0000]

# to find A such that A'*A = Sigma, we find the eigendecompositon of Sigma:
# V'*D*V = Sigma
D, V = eig(Sigma)
A = V .* (D .^ 0.5);

# But a quicker way if Sigma has many dimensions is to use the Cholesky 
# decomposition instead:
# C = chol(Sigma)
# which then gives C' * C == Sigma

# The standard multivariate normal is just a series of independently drawn
# univariate normal random numbers
standard_random_vector = randn(2, 1)
 
# this is converted to an arbitrary multivariate normal by:
normal_random_vector = Mu .+ A * standard_random_vector
 
# so now we can generate many of them, and plot what's going on:
standard_random_vector = randn(2, 1000)
normal_random_vector = Mu .+ A * standard_random_vector

draw(PNG("scatter_std_gaussian.png", 50cm, 36cm),
  plot(x = standard_random_vector[1,:], y = standard_random_vector[2,:]));

draw(PNG("scatter_normal_gaussian.png", 50cm, 36cm),
  plot(x = normal_random_vector[1,:], y = normal_random_vector[2,:]));

#####################################################################
# 3. A Kernel Matrix demonstration Brownian motion
#####################################################################
n=40

# The brownian motion kernel is defined through its very simple inverse
inverse = 2 .* eye(n)
inverse = inverse - (triu(ones(n,n),1) - triu(ones(n,n),2))'
inverse = inverse - (triu(ones(n,n),1) - triu(ones(n,n),2))

# subplot(1,3,1);imagesc(inverse);title('Brownian motion K^{-1}')
# subplot(1,3,2);imagesc(inv(inverse));title('Brownian motion K')
#draw some samples
#r=randn(n,10)
#subplot(1,3,3);plot(inverse\r);title('10 samples')

