# We want to explore fast ways to subset a large, square N*N matrix in julia dependent on
#   1. N = number of dimensions of large matrix / data set
#   2. S = subset size = floor(log(N)) - floor(10*log(N))
#   3. Method used to extract subset of dimensions/features

function extract_with_getindex(matrix, subset)
  matrix[subset,:][:,subset]
end

function extract_with_getindex2(matrix, subset)
  matrix[subset,subset]
end

function extract_in_loop(matrix, subset)
  m = length(subset)
  res = zeros(m, m)
  for i in 1:m
    si = subset[i]
    for j in 1:m
      res[i, j] = matrix[si, subset[j]]
    end
  end
  res
end

# Do some runs to ensure everything is compiled
m = randn(10,10)
ss = sort(shuffle(collect(1:10))[1:4])
r1 = extract_in_loop(m, ss)
r2 = extract_with_getindex(m, ss)
r3 = extract_with_getindex2(m, ss)

num_features(n) = rand(int(floor(log(n))):int(floor(10*log(n))))

NumRepeats = 10
for n in [100, 1000, 10000]
  times1 = zeros(NumRepeats)
  times2 = zeros(NumRepeats)
  times3 = zeros(NumRepeats)
  for i in 1:NumRepeats
    matrix = randn(n, n)
    ss = sort(shuffle(collect(1:n))[1:num_features(n)])
    tic()
      res1 = extract_in_loop(matrix, ss)
    times1[i] = toq()
    tic()
      res2 = extract_with_getindex(matrix, ss)
    times2[i] = toq()
    tic()
      res3 = extract_with_getindex2(matrix, ss)
    times3[i] = toq()
    if res1 != res2 || res1 != res3 || res2 != res3
      println("Results differ!")
    end
  end
  println("n = $(n), extract_in_loop = $(mean(times1)), extract_with_getindex = $(mean(times2)), extract_with_getindex2 = $(mean(times3)), getindex2 factor faster = $(mean(times1)/mean(times3))")
end