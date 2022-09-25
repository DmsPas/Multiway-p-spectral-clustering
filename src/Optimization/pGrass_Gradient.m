function [out, V] = pGrass_Gradient(V,W,p,constraint_first,normalized)
% Gradient of the clustering algorithm presented in
% 
% "Multiway p-spectral graph cuts on Grassmann manifolds",
% by Pasadakis, D., Alappat, C.L., Schenk, O., Wellein, G,
%
% published in Machine Learning 111, 791–829 (2022).


global Hess_store;

N = size(V.main,1);
K = size(V.main,2);
gradient = zeros(N,K);
lambda   = 1E-10;
% max_scaling = 0;

deg=full(sum(W));% row vector
% loop over k
for k = 1:K
    x = V.main(:,k);
    % call Functional as per B&H
    objective     = V.funcarray(k);
    if(k~=1 || ~constraint_first)        
        gradient(:,k) = Gradient(x,W,p,normalized,deg,objective);
    end
    Hess            = Hessian(x,W,p,normalized,deg);    
    scaling         = 1;
    Hess            = Hess/scaling;    
    Hess_scale      = Hess + lambda*speye(N);
    Hess_store{k}   = Hess_scale; 
end


if(constraint_first)
    gradient(:,1) = V.main * V.main' * gradient(:,1);
end

out.main = gradient;

end
