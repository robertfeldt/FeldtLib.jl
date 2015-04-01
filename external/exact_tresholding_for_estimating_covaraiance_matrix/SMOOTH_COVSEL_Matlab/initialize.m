function [alpha,beta,L,N,sigma1] = initialize(U0,Sigma,rho,in_alpha,in_beta,eps)

n = max(size(Sigma));
s = eig(Sigma);
trace = sum(s);

% estimate alpha
if (in_beta == inf) && (in_alpha ==0)
  alpha = 1.0/(s(n)+n*rho);
else
  alpha = in_alpha;
end

% estimate beta
if in_beta == inf
  beta = min([(n-alpha*trace)/rho, (n-rho*sqrt(n)*alpha)/s(1)]);
else
  beta = in_beta;
end

% setup Lipschitz constant
L = (rho*beta)^2;

% d_1(U) = |U|^2/2, D1 = max |U-U0|^2/2 %
sigma1 = 1.0;
D1 = 0;
for i=1:n
  for j=1:n
     D1 = D1 + max([(1-U0(i,j))^2, (1+U0(i,j))^2]);
  end
end
D1 = D1/2.0;

% worst-case iterate complexity based on rarely rough beta
N = 2*sqrt(L*D1/(sigma1*eps))-1;
