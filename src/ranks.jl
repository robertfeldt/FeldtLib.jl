function ranks(v::AbstractArray)
  n     = length(v)
  place = sortperm(v)
  ord   = Array(Float64, n)

  for(i in 1:n)
    ord[place[i]] = i
  end

  ord
end