function [Time, X,S,C] = compute_GL_exact_smacs(Sigma_cov,rho,MAX_COMP_SIZE,maxiter,freqprint,epstol)



p=size(Sigma_cov,1);

cov_thr_graph= sparse(abs(Sigma_cov) > rho);
[S, C] = graphconncomp(cov_thr_graph,'Directed',false);

max_comp_size=0;
for ii=1:S
   max_comp_size=max(sum(C==ii),max_comp_size);
end
   

if (max_comp_size > MAX_COMP_SIZE)
   disp(['Size of largest connected component exceeds',MAX_COMP_SIZE])
   disp('Exiting code, try larger value of rho')
   Time=[];X=[]; S=[]; C=[]; 
   return    
end


Time=zeros(S,1); Iter=Time;
X=sparse(p,p);


for ii=1:S
   ids=(C==ii);
 
   % if size of conn-comp is larger than 1, then call smacs
   if (sum(ids)>1)
   t=tic;
   [tmp, Obj, Gap, Iter(ii), T] =smoothcovsel_main(Sigma_cov(ids,ids), rho,0,inf, maxiter, freqprint, epstol);
   X(ids,ids)= sparse(tmp); clear tmp;
   Time(ii)=toc(t);
   
   % if size of conn-comp is 1, directly write what value of diagonal is
   else
   X(ids,ids)= Sigma_cov(ids,ids) + rho; 
   end
   
 
   clear ids
      
end
