function [objective, V] = pGrass_Functional(V,W,p,normalized)
% Functional of the clustering algorithm presented in
% 
% "Multiway p-spectral graph cuts on Grassmann manifolds",
% by Pasadakis, D., Alappat, C.L., Schenk, O., Wellein, G,
%
% published in Machine Learning 111, 791–829 (2022).


N = size(V.main,1);
K = size(V.main,2);
objective = 0;
functional_array = zeros(K,1);

deg=full(sum(W));% row vector
% loop over k
for k = 1:K
    x = V.main(:,k);
    % call Functional as per B&H
    temp = Functional(x,W,p,normalized,deg);
    objective = objective + temp;
    functional_array(k) = temp;
end

V.funcarray = functional_array;
V.objective = objective;

end
