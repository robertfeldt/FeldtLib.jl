# This only works if v is a n*1 matrix and the op is *.
function Base.broadcast{Tv,Ti}(op, v::Array{Tv,2}, A::SparseMatrixCSC{Tv,Ti})
  I, J = findn(A)
  V = zeros(nnz(A))
  for(l in 1:nnz(A))
    row = I[l]
    V[l] = op(v[row], A.nzval[l])
  end
  sparse(I, J, V)
end

n = 100
m = 42
v = randn(n,1)
s = sprandn(n,m,0.75)
res = broadcast(*, v, s)

for(row in 1:n)
  for(col in 1:m)
    @assert res[row, col] == (v[row] * s[row, col])
  end
end