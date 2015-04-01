# In a recent IEEE Trans. on Evolutionary Computation paper by Gary Yen et al
# it is proposed that a double elimination tournament is held between pareto
# approximation fronts in order to rank the EA's producing these fronts.
# Even though the proposed method is interesting and seems to give robust
# comparison results the paper seems to lack a comparison to some base method
# for ranking. Here I want to study this problem in the abstract using
# simulations and statistics. My base hypothesis is that there are simpler
# and more straightforward ways to accomplish this type ranking.

# The basic setting is that of a pool of performance metrics on which a pool
# of solutions should be compared. Each metric calculates a single, real value
# for each front. Thus we can summarize this information in a matrix where
# each row corresponds to a metric and each column to a front that originates
# from one EA. Lets fix some numbers for our simulations:
Nmetrics = 4
NEAs = 5
Nfronts = 10 # 50 fronts per EA

# Without loss of generality, assume that we normalize all the performance metrics
# to the [0.0, 1.0] range and that higher values indicate better performance.
# Thus we can summarize all the info we have in a randomly sampled matrix where
# each column comes from one of the EA's:
rand_performance() = (rand(Nmetrics, NEAs * Nfronts), [rand(1:5) for i in 1:(NEAs*Nfronts)])

# To hold one tournament among N candidates we randomly select one metric and 
# then have N/2 binary tournaments, a pool is an array of column indices into an
# original performance matrix
function hold_tournament(P, pool)
  N = length(pool)
  pool_indices = shuffle(collect(1:N))
  winners = zeros(Int64, int(ceil(N/2)))
  losers = zeros(Int64, int(ceil(N/2)))
  for i in 1:int(floor(N/2))
    i1 = pool[pool_indices[2*i-1]]
    i2 = pool[pool_indices[2*i]]
    metric = rand(1:Nmetrics)
    if P[metric, i1] > P[metric, i2]
      winners[i] = i1
      losers[i] = i2
    else
      losers[i] = i1
      winners[i] = i2
    end
  end

  # Add the odd man out to both winners and losers bracket as stated in their paper.
  if mod(N, 2) == 1
    losers[end] = winners[end] = pool[pool_indices[end]]
  end

  return winners, losers
end

# Lets generate one example:
P, EAs = rand_performance()

# And lets run a double elimination tournament to rank the EAs:
function double_elimination_tournament(P, EAs)
  # Initial round
  winners, losers = hold_tournament(P, collect(1:length(EAs)))
  while length(winners) > 1
    winners, wl = hold_tournament(P, winners, metric)
    lw, eliminated = hold_tournament(P, losers, metric)
    losers, eliminated = hold_tournament(P, [wl, lw], metric)
  end
  winners[1]
end

# Lets run many tournaments and collect stats on who wins the most:
function run_many_tournaments(P, EAs, numtournaments = 1000)
  counts = Dict{Any,Int64}()
  for i in 1:numtournaments
    winner = double_elimination_tournament(P, EAs)
    key = (winner, EAs[winner])
    counts[key] = get(counts, key, 0) + 1
  end
  counts
end

function rank_vector(v, rev=false)
  p = sortperm(v, rev = rev)
  ranks = zeros(Int64, size(v))
  for i in 1:length(v)
    ranks[p[i]] = i
  end
  ranks
end

# We can compare these results to methods based on ranks per metric.
function rank_per_row(matrix, rev = false)
  ranks = zeros(Int64, size(matrix))
  for row in 1:size(matrix, 1)
    ranks[row, :] = rank_vector(matrix[row, :][:], rev)
  end
  ranks
end

ranks = rank_per_row(P, true)
meanranks = mean(ranks, 1)[:]
