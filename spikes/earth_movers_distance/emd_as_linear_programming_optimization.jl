using JuMP
using SCS

w1 = [0.4, 0.2, 0.2, 0.1, 0.1]
w2 = [0.6, 0.2, 0.2]

C = [ 3 5 2;
      0 2 5;
      1 1 3;
      8 4 3;
      7 6 5]

# The EMD is the sum(T .* C) where T[i, j] from each word/elem i in w1 to 
# each word/elem j in w2 where T is the matrix that minimizes sum(T .* C)
# and sum(T, 1) = w2 and sum(T, 2) = w1.

# T = flow = emd(w1, w2, C)

m = Model()
N = length(w1)
M = length(w2)
@variable(m, 0 <= T[1:N,1:M] <= 1)
for j in 1:M
    @constraint(m, sum(T[i,j] for i=1:N) == w2[j])
end
for i in 1:N
    @constraint(m, sum(T[i,j] for j=1:M) == w1[i])
end
@objective(m, Min, sum(T .* C))
print(m)
@time status = solve(m)
println("Objective value: ", getobjectivevalue(m))
getvalue(T)