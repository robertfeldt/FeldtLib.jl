const Rmath = "libRmath-julia"

function rbeta(alpha::Float64, beta::Float64)
  ccall( (:rbeta, Rmath), Float64, (Float64, Float64), alpha, beta )
end

function rbeta(n::Int, alpha::Float64, beta::Float64)
  res = zeros(Float64, n)
  for i in 1:n
    res[i] = rbeta(alpha, beta)
  end
  res
end

#function pbeta(x::Float64, alpha::Float64, beta::Float64)
#  ccall( (:pbeta, Rmath), Float64, (Float64, Float64, Float64, Int64, Int64), x, alpha, beta, 0, 0)
#end

# Incomplete beta function.
#function betainc(x, a, b)
#  pbeta(float(x), float(a), float(b))
#end
