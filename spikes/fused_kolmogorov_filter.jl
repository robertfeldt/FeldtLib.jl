# Fused Kolmogorov Filter for Variable Screening as described in the paper Mai2014:
#  http://arxiv.org/pdf/1403.7701v2.pdf
#
# Julia implementation by Robert Feldt, robert.feldt@gmail.com
#

# Kolmogorov-Smirnov Statistic, i.e. supremum of differences between empirical cdfs.
function ksstatistic{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})
  nx, ny = length(x), length(y)
  sortidx = sortperm([x; y])
  pdf_diffs = [ones(nx)/nx; -ones(ny)/ny][sortidx]
  cdf_diffs = cumsum(pdf_diffs)
  deltap = maximum(cdf_diffs)
  deltan = -minimum(cdf_diffs)
  delta = max(deltap, deltan)
  return delta
end

# Rank variables by their fused kolmogorov filter statistic.
#   varrank, fkfstats, kfstats = fused_kolmogorov_filter(X, y)
function fused_kolmogorov_filter{T <: Real, S <: Real}(X::AbstractMatrix{T}, y::AbstractVector{S};
  partitions = 3:6)

  # Convert vector of bools to an int vector
  if typeof(y[1]) <: Bool
    y = map(b -> b ? 2 : 1, y)
  end

  # Calc actual number of partitions to use based on number of unique y values
  num_ys = length(unique(y))
  minpartitions = min(minimum(partitions), num_ys)
  maxpartitions = min(maximum(partitions), num_ys)
  @assert minpartitions > 1
  numpartitions = maxpartitions - minpartitions + 1

  n,p = size(X)

  # Now loop and calc the fkfstatistics per variable (column) in X
  kfstats = zeros(numpartitions, p)
  for G in minpartitions:maxpartitions
    splitvalues = quantile(y, cumsum(ones(G)/G)[1:end-1])
    slice_idxs = slice_indices(y, splitvalues)
    idxsperslice = map(s -> find(i -> i == s, slice_idxs), 1:G)
    gindex = G - minpartitions + 1
    for pi in 1:p
      xs = X[:,pi]
      # For all pairwise combinations of slices we calc the ksstatistic, the max of them
      # is the Kolmogorov filter statistic
      kfs = map(combinations(1:G, 2)) do slices
        lidxs, midxs = idxsperslice[slices[1]], idxsperslice[slices[2]]
        ksstatistic(xs[lidxs], xs[midxs])
      end
      kfstats[gindex, pi] =  maximum(kfs)
    end
  end

  # Fuse the kfstats by summing
  fkfstats = sum(kfstats, 1)[:]

  return sortperm(fkfstats, rev=true), fkfstats, kfstats
end

# For each value in y find the slice it belongs to given the max slice value
# in splitvalues.
function slice_indices{T <: Real, S <: Real}(y::AbstractVector{T}, splitvalues::AbstractVector{S})
  n = length(y)
  nsplits = length(splitvalues)
  idxs = zeros(Int, n)
  for i in 1:n
    j = 1
    while j <= nsplits && y[i] >= splitvalues[j]
      j += 1 
    end
    idxs[i] = j
  end
  idxs
end

# Test problem is a simple linear model:
N = 100
P = 1000
X = randn(N, P)
y = 10*X[:,1] + 10*X[:,2] + randn(N)
@time xrank, fkfs, kfs = fused_kolmogorov_filter(X, y; partitions = 3:7);
@show xrank[1:5]
