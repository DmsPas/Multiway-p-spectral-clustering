function [Dunn_curr] = dunn_index_calculation(V_transpose,x_kmeans,cluster_centers)

Num_Clust = length(unique(x_kmeans));
k         = Num_Clust;

% Dunn index computation
inertia        = 0;
max_intra_dist = 0;
%Calculate intra-cluster distance
for clust=1:Num_Clust
    idx  = find(x_kmeans==clust);
    pts  = V_transpose(:,idx);
    pts  = pts';
    dist = pdist2(pts, cluster_centers(clust,:));
    intra_dist = sum(dist);
    inertia    = inertia + intra_dist;
    max_intra_dist = max(max_intra_dist, intra_dist);
end
%Calculate inter-cluster dist
inter_dist     = pdist2(cluster_centers, cluster_centers);
max_inter_dist = max(inter_dist(:));
%add this to diag so it doesn't make self cluster small
inter_dist     = inter_dist + eye(k)*max_inter_dist;
min_inter_dist = min(inter_dist(:));

Dunn_curr   = min_inter_dist / max_intra_dist;

end
