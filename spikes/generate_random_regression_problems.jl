using Distributions
using DataFrames

# Randomly generate covariates, response and errors for a regression problem.
#
#  (A, b, errors) = gen_gaussian_regression_problem(N, P, coefs, probdist, snrDB)
#
# Outputs:
#  A - N*P matrix where each row has P covariates
#  b - N*1 response vector (actual response + errors)
#  errors - N*1 error vector
#
# Inputs:
#  coefs - coefficients for the regression problem
#  probdist - probability ditribution for generating the covariates
#  snrDB - Signal-noise ratio in decibels
#
function gen_regression_problem(N, P, coefs, probdist = Normal(0.0, 1.0), snrDB = 20.0)
  A = zeros(Float64, N, P)
  for i in 1:N
    A[i,:] = rand(probdist, P)
  end
  A

  y = A * coefs

  noisestd = sqrt( noisevar(snrDB, var(y)) )
  noisedist = Normal(0.0, noisestd)
  errors = rand(noisedist, N)

  b = y .+ errors

  return (A, b, errors)
end

noisevar(snrDB, signalvar) = signalvar / (10 ^ (snrDB / 10))

# coefs = rand(5)
# A, b, e = gen_regression_problem(20, 5, coefs)

### Evaluating quality of regression coefficients ###

function confusionmatrix(expected, actual)
  @assert length(expected) == length(actual)
  cm = zeros(2, 2)
  for i in 1:length(expected)
    e = abs(expected[i]) > 0.0 ? 1 : 2
    a = abs(actual[i])   > 0.0 ? 1 : 2
    cm[a, e] += 1.0
  end
  return cm
end

tp(cm) = cm[1,1]
fp(cm) = cm[1,2]
fn(cm) = cm[2,1]
tn(cm) = cm[2,2]
sensitivity(cm) = cm[1,1] / (cm[1,1] + cm[2,1])
accuracy(cm) = (tp(cm) + tn(cm)) / (sum(cm))
mcc(cm) = (tp(cm) * tn(cm) - fp(cm) * fn(cm)) / sqrt( (tp(cm) + fp(cm)) * (tp(cm) + fn(cm)) * (tn(cm) + fp(cm)) * (tn(cm) + fn(cm)))

numactive(coefs) = sum(abs(coefs) .> 0.0)
positivecoefs(p) = 1.0 + rand(p)
normalcoefs(p) = randn(p)

gencoefs(p, k, coefgen = positivecoefs) = shuffle([coefgen(k), zeros(p-k)])

# regressors is a dict from a regressor name to functions from (matrix, vector) to vector of coefficients.
function compare_regressors(regressors, coefgenerator = positivecoefs;
  Ns = [10, 50, 100], 
  Pratios = [0.50, 1.0, 2.0], 
  Ks = [1, 5, (p) -> sqrt(p)],
  SNRs = [10, 50], 
  NumReps = 1,
  CovarProbDists = [Normal(0.0, 1.0)]
  )

  startime = time()

  problemnum = 0

  # For each experiment with a regressor we save its Accuracy, MCC, MSE, CoeffNorm, NumSelectedCoefs
  df = DataFrame()

  for rep in 1:NumReps
    for N in Ns
      for Pratio in Pratios
        P = ifloor(Pratio*N)
        ks = map(Ks) do K
          if typeof(K) <: Function
            K(P)
          else
            K
          end
        end

        for k in ks
          k = ifloor(k)
          # Skip if k is not smaller than P
          if k >= P
            next
          end

          for SNR in SNRs
            for covarprobdist in CovarProbDists

              activecoefs = coefgenerator(k)
              coefs = shuffle([activecoefs, zeros(P-k)])
              A, b, errors = gen_regression_problem(N, P, coefs, covarprobdist, SNR)
              problemnum += 1

              for (name, reg) in regressors

                println("$name: N = $N, P = $P, k = $k, SNR = $SNR")
                tic()
                  betas = reg(A, b)
                t = toq()
                msev = mse(b, A * betas)
                println("  time taken: $t, mse = $msev")

                cm = confusionmatrix(coefs, betas)

                res = DataFrame(
                  # Problem description
                  N = N, P = P, k = k, SNRdb = SNR, CovariateGen = string(covarprobdist), PId = problemnum,

                  Time = t,
                  Accuracy = accuracy(cm),
                  MCC = mcc(cm),
                  nSelected = numactive(betas),
                  AvgAbsNorm = norm(coefs .- betas, 1)/P,
                  Regressor = name,
                  MSE = msev

                )

                df = rbind(df, res)

              end
            end
          end
        end
      end
    end
  end

  filename = strftime("compare_regressors_%Y%m%d_%H%M%S.csv", startime)
  writetable(filename, df)

  return df
end
