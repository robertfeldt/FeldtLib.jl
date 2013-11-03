function ranks(v::AbstractArray)
  n     = length(v)
  place = sortperm(v)
  ord   = Array(Float64, n)

  for(i in 1:n)
    ord[place[i]] = i
  end

  ord
end

function baumgartner_weis_schindler_statistic(x, y)
  # Create a single array for sorting the values.
  values = vcat(x, y)

  # Get the ranks. The first 1:n are ranks for x, (n+1):(n+m) for y.
  rs = ranks(values)

  FeldtLib.bws_statistic_from_ranks(rs, length(x), length(y))
end

function bws_statistic_from_ranks(ranks, n, m)
  # Pre-calc some values needed in both b calcs.
  nm = n+m
  np1 = n + 1
  mp1 = m + 1

  # Now calc the two b values
  b_x = (1/n) * (n / (m*nm)) * calc_bws_sum(ranks,   1,  n, n, np1, nm, 0)
  b_y = (1/m) * (m / (n*nm)) * calc_bws_sum(ranks, np1, nm, m, mp1, nm, n)

  (b_x + b_y) / 2
end

function calc_bws_sum(ranks, mini, maxi, n, np1, nm, subi = 0)
  sum = 0.0
  for(i in mini:maxi)
    k = i - subi
    nom = ( ranks[i] - (nm/n*k) )^2
    denom = (k/np1) * (1 - (k/np1))
    sum += (nom / denom)
  end
  sum
end