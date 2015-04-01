function ranks(v::AbstractArray)
  n     = length(v)
  place = sortperm(v)
  ord   = Array(Int64, n)

  for(i in 1:n)
    ord[place[i]] = i
  end

  ord
end

function factorial(n::Integer)
  prod(1:n)
end

function combinations(k::Integer, n::Integer)
  factorial(n) / (factorial(k) * factorial(n-k))
end