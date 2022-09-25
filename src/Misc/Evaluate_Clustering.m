function [RCut,clusters] = Evaluate_Clustering(V,W,label)
% Calculate the accuracy of the clustering and the RatioCut 
% for a decreasing p value
% Input:
% Num_Eigs: # of eigenvectors to be considered
% V: the p-spectral eigenvectors
% W: adjacency matrix of the graph
% label: true labels of the nodes
% Output:
% idx: indices of these solutions
% RCut: Ratio cut 
% clusters: vector of node assignments
% 

%% Initialize
K        = size(V,2);

%% perform orthogonal k-means
num_ortho  = 20;
num_random = 10; 

[RCut,clusters]  = kmeans_orthogonal(V(:,1:K)', K, num_ortho, num_random,label,W);


end