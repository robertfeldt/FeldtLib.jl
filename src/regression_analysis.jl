using Distributions

# A covariate expander is a function of one and more covariates and can describe 
# its function in a string.
type CovariateExpander
	name::Function # Map the name of a covar to a string describing the expanded covar
	func::Function # Map covar value(s) to expanded values
end

expand(ce::CovariateExpander, data) = ce.func(data)
new_name(ce::CovariateExpander, name) = ce.name(name)

# A DataSet is a float array with named covariates. Each covariate
# corresponds to a row of the data set.
type DataSet
	names::Vector{String}
	data::Array{Float64, 2}
end

num_covars(d::DataSet) = size(d.data, 1)
length(d::DataSet) = size(d.data, 2)

function random_data(covars = 2, n = 30, offset = 10.0)
	names = [string("x", i) for i in 1:covars]
	DataSet(names, offset .+ randn(covars, n))
end

function expand!(d::DataSet, ce::CovariateExpander, indexToExpand)
	d.data = cat(1, d.data, expand(ce, d.data[indexToExpand,:]))
	push!(d.names, new_name(ce, d.names[indexToExpand]))
end
expand!(d::DataSet, ces::Vector{CovariateExpander}, indexToExpand) = map((ce) -> expand!(d, ce, indexToExpand), ces)

LogCE = CovariateExpander((n) -> "log($(n))", (d) -> log(d))
SqrtCE = CovariateExpander((n) -> "sqrt($(n))", (d) -> sqrt(d))
NLogNCE = CovariateExpander((n) -> "$(n)*log($(n))", (d) -> d .* log(d))
Pow2CE = CovariateExpander((n) -> "$(n)^2", (d) -> d.^2)
Pow3CE = CovariateExpander((n) -> "$(n)^3", (d) -> d.^3)

InterceptName = "(Intercept)"
InterceptCE = CovariateExpander((n) -> InterceptName, (d) -> ones(d))

StandardCES = [LogCE, SqrtCE, NLogNCE, Pow2CE, Pow3CE]

function rand_in_range(rows, columns, minRange = 0.0, maxRange = 100.0)
	minRange .+ (maxRange - minRange) .* rand(rows, columns)
end

# Create a random data set sampled from a uniform range.
function random_data(minRange = 0.0, maxRange = 100.0, numCovars = 2, n = 30;
	expandCovars = true
	)

	names = [string("x", i) for i in 1:numCovars]
	d = DataSet([InterceptName, names], 
		cat(1, ones(Float64, 1, n), rand_in_range(numCovars, n, minRange, maxRange)))
	if expandCovars
		[expand!(d, StandardCES, i+1) for i in 1:numCovars]
	end
	d

end

# A regression problem is a dataset, a linear function of some subset of that data set and
# an error level (signal-to-noise ratio).
type RegressionProblem
	d::DataSet
	coefficients::Array{Float64,2}
	snr::Float64 # Mean of signal (y) over the standard deviation of the error
end

function random_regression_problem(numCovars = 2, numActiveCovars = 1; 
	n = 30, 
	minRange = 0.0, 
	maxRange = 100.0,
	expandCovars = true,
	snr = 20.0
	)

	d = random_data(minRange, maxRange, numCovars, n; expandCovars = expandCovars)

	coefs = zeros(Float64, 1, num_covars(d))
	[coefs[i] = 10.0*rand() for i in shuffle(collect(1:num_covars(d)))[1:numActiveCovars]]
	p = RegressionProblem(d, coefs, snr)
	y, errors = sample_y(p)
	return p, y, errors

end

# Get the model encoded in a regression problem as a string.
function model(rp::RegressionProblem, coefs = rp.coefficients)
	mstr = ["y"]
	sep() = last(mstr) == "y" ? " = " : " + "
	for i in 1:size(coefs, 2)
		n = rp.d.names[i]
		c = coefs[i]
		if c != 0.0
			if n == InterceptName
				mstr = [mstr, sep(), @sprintf("%.4f", c)]
			else
				mstr = [mstr, sep(), @sprintf("%.4f", c), " * ", n]
			end
		end
	end
	join(mstr, "")
end

import Base.length
length(r::RegressionProblem) = size(r.d.data, 2)

y_actual(r::RegressionProblem) = r.coefficients * r.d.data

function sample_y(r::RegressionProblem)
	y = y_actual(r)
	mu = mean(y)
	sd = mu / r.snr
	errors = rand(Normal(0.0, sd), length(r))'
	return( y .+ errors, errors )
end

# Analyze a regression model proposed for a certain problem. Will calculate
# RSS, AIC, False positives (FP) and False negatives (FN) for the model.
function analyze_regression(rp::RegressionProblem, model::Array{Float64,2})
	report = [
		"M:             $(num_covars(rp.d))",
		"N:             $(length(rp))",
		"Active covars: $(num_nonzero(rp.coefficients))",
		"SNR (target):  $(rp.snr)"
	]
	println(join(report))
end