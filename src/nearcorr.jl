# This is based on Higham's own Matlab implementation based on his paper:
#   By N. J. Higham, 13/6/01, updated 30/1/13.
#   Reference:  N. J. Higham, Computing the nearest correlation
#   matrix---A problem from finance. IMA J. Numer. Anal.,
#   22(3):329-343, 2002.
#
# as presented on: http://nickhigham.wordpress.com/2013/02/13/the-nearest-correlation-matrix/
# It is typically not the fastest algorithm but should be ok for small (<500) dimensions...

# NEARCORR    Nearest correlation matrix.
#    (X, iter) = NEARCORR(A,TOL,FLAG,MAXITS,N_POS_EIG,W,PRNT)
#    finds the nearest correlation matrix to the symmetric matrix A.
#    TOL is a convergence tolerance, which defaults to 16*EPS.
#    If using FLAG == 1, TOL must be a 2-vector, with first component
#    the convergence tolerance and second component a tolerance
#    for defining "sufficiently positive" eigenvalues.
#    FLAG = 0: solve using full eigendecomposition (EIG).
#    FLAG = 1: treat as "highly non-positive definite A" and solve
#              using partial eigendecomposition (EIGS).
#    MAXITS is the maximum number of iterations (default 100, but may
#    need to be increased).
#    N_POS_EIG (optional) is the known number of positive eigenvalues of A.
#    W is a vector defining a diagonal weight matrix diag(W).
function nearcorr(A, tol = [16*Base.eps()], flag = 0, maxits = 100, 
   n_pos_eig = false, w = ones(size(A, 1),1))

   if A != A'
      throw( ArgumentError("The matrix must be symmetric.") )
   end
   
   if length(tol) == 1 && flag == 1
      tol = size(A, 1) * Base.eps() * [1.0, 1.0]
   end
 
   n = size(A, 1)
   if flag == 1 && n_pos_eig == false
      d, V = eig(A)
      n_pos_eig = sum(d >= tol[2] * d[n])
   end

   X = copy(A)
   Y = copy(A)
   iter = 1
   rel_diffXY = rel_diffY = rel_diffX = Inf
   dS = zeros(size(A)) 
   w = w[:]
   Whalf = sqrt(w*w')

   while max(rel_diffX, rel_diffY, rel_diffXY) > tol[1]
 
      Xold = X
      R = X - dS
      R_wtd = Whalf .* R

      if flag == 0
         X = proj_spd(R_wtd)
      elseif flag == 1
         X, np = proj_spd_eigs(R_wtd, n_pos_eig, tol[2])
      end

      X =  X ./ Whalf
      dS = X - R
      Yold = Y
      Y = proj_unitdiag(X)
      normy = normfro(Y)
      rel_diffX = normfro(X-Xold)/normfro(X)
      rel_diffY = normfro(X-Xold)/normy
      rel_diffXY = normfro(Y-X)/normy

      iter = iter + 1;
      if iter > maxits
         throw( OverflowError("Stopped after $(maxits) iterations. Try increasing maxits.") )
      end

      X = Y
 
   end

   return X, iter
end

function proj_spd(A)
   d, V = eig(A)
   d[d .< 0.0] = 0.0 # Set negative eigen values to 0
   B = V * diagm(d) * V'
   (B+B')/2 # Ensure symmetric
end

function proj_spd_eigs(A, n_pos_eig, tol)
   k = n_pos_eig + 10 # 10 is safety factor.
   if k > size(A, 1)
      k = n_pos_eig
   end
   d, V = eigs(A, k; which = "LM") 
   j = d .> tol*maximum(d)
   n_pos_eig_found = sum(j)
   B = V * diagm(d[j]) * V'
   B = (B+B')/2 # Ensure symmetric
   return B, n_pos_eig_found
end

function proj_unitdiag(A)
   B = copy(A)
   B[diagind(B)] = 1.0
   B
end
