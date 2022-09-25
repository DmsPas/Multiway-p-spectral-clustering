function [out,V] = pGrass_Hessian(V,eta,W,p)
% Hessian of the clustering algorithm presented in
% 
% "Multiway p-spectral graph cuts on Grassmann manifolds",
% by Pasadakis, D., Alappat, C.L., Schenk, O., Wellein, G,
%
% published in Machine Learning 111, 791–829 (2022).



global Hess_store;

N = size(V.main,1);
K = size(V.main,2);
hessian = zeros(N,K);

% loop over k
for k = 1:K
% For a symmetric matrix Ax = A'x. The right hand side is 
% more efficient in MATLAB, due to its CSC format.
    hessian(:,k) = Hess_store{k}'*eta.main(:,k);    
end

out.main = hessian;

end