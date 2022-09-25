function [L,Diag,vw] = CreateLapl(W, normalized)
% Create the Laplacian matrix of a graph (normalized or combinatorial)
% Input
% W    : Adjacency matrix
% Output
% L    : Laplacian



% Degree matrix
Diag = zeros(size(W,1));
for i = 1:size(W,1)
   Diag(i,i) = sum(W(:,i)); 
end
% number of points
NumNodes  = size(W,1);
% Vertex weights
vw         = sum(W,2); 

if normalized == 0
    % Combinatorial Laplacian L = D-W
    L = Diag - W;
    L = sparse(L);
    elseif normalized == 1
    % Symmetrically normalized Laplacian L_sym = D^{-1/2} * L * D^{-1/2}   
    VWsqrt = spdiags(1./(vw.^0.5), 0, NumNodes, NumNodes);
    L      = VWsqrt*(spdiags(vw,0,NumNodes,NumNodes)-W)*VWsqrt; 
end
