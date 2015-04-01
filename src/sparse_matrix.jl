# This only works if v is a n*1 matrix and the op is *.
function Base.broadcast{Tv,Ti}(op, v::Array{Tv,2}, A::SparseMatrixCSC{Tv,Ti})
  I, J = findn(A)
  V = zeros(nnz(A))
  vn, vm = size(v)
  if vn >= 1 && vm == 1
    for(l in 1:nnz(A))
      row = I[l]
      V[l] = op(v[row], A.nzval[l])
    end
  elseif vn == 1 && vm >= 1
    for(l in 1:nnz(A))
      col = J[l]
      V[l] = op(v[col], A.nzval[l])
    end
  else
    throw(ArgumentError("invalid dimensions"))
  end
  sparse(I, J, V)
end

n = 100
m = 42

v = randn(n,1)
s = sprandn(n,m,0.75)
@assert full(broadcast(*, v, s)) == broadcast(*, v, full(s))

v = randn(1,n)
s = sprandn(m,n,0.42)
res = broadcast(*, v, s)
@assert full(broadcast(*, v, s)) == broadcast(*, v, full(s))