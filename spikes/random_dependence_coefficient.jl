using CrossDecomposition

# This implements the RDC as described in the paper:
#  http://arxiv.org/pdf/1304.7717v2.pdf
function rdc(x::AbstractMatrix, y::AbstractMatrix, k = 10, s = 1/6, f1 = cos, f2 = sin)
  tx = hcat(mapslices((u) -> sortperm(u)/length(u), x, 1), ones(size(x,1)))
  ty = hcat(mapslices((u) -> sortperm(u)/length(u), y, 1), ones(size(y,1)))
  wx = s/size(tx,2) * randn(size(tx,2), k)
  wy = s/size(ty,2) * randn(size(ty,2), k)
  cor(canoncor(hcat(f1(tx*wx), f2(tx*wx)), hcat(f1(tx*wx), f2(tx*wx))))[1]
end

x = randn(100,1)
y = 10x
rdc(x,y,10,1/6)

# R code:
#rdc <- function(x,y,k,s) {
#  tx <- cbind( apply(as.matrix(x),2,function(u) rank(u)/length(u)), 1 )
#  ty <- cbind( apply(as.matrix(y),2,function(u) rank(u)/length(u)), 1 )
#  wx <- matrix(rnorm(ncol(tx)*k,0,s),ncol(tx),k)
#  wy <- matrix(rnorm(ncol(ty)*k,0,s),ncol(ty),k)
#  cancor(cbind(cos(tx%*%wx),sin(tx%*%wx)), cbind(cos(ty%*%wy),sin(ty%*%wy)))$cor[1]
#}