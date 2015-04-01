include("adaptive_iterative_tresholding.jl")
include("generate_random_regression_problems.jl")

N, P = (50, 20)
coefs = gencoefs(20, 6)
A, b, errors = gen_regression_problem(50, 20, coefs)

ait_fixed_scad_reg(A, b) = adaptive_iterative_tresholding(A, b; f = scad)[1]
ait_fixed_hard_reg(A, b) = adaptive_iterative_tresholding(A, b; f = hard)[1]
ait_fixed_scad_reg(rand(10, 2), rand(10)) # Ensure compiled
ait_fixed_hard_reg(rand(10, 2), rand(10)) # Ensure compiled

ait_scad_bestk_reg(A, b) = begin
  res = ait_best_k(A, b; f = scad)
  for k in sort(collect(keys(res[3])))
    println(k, " => ", res[3][k])
  end
  res[1]
end
ait_scad_bestk_reg(rand(100, 10), rand(100)) # Ensure compiled

df = compare_regressors({
#  "AIT/scad fixed k=sqrt(P)" => ait_fixed_scad_reg,
  "AIT/scad best k" => ait_scad_bestk_reg,
#  "AIT/hard fixed k=sqrt(P)" => ait_fixed_hard_reg,
})
