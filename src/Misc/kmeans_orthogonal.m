function [Cut,cluster,cluster_centers] = kmeans_orthogonal(xx, k, num_ortho, num_random, label, W, normalized)
% V(:,1:K)', K, num_ortho, num_random, label, W


p = size(xx,1);
n = size(xx,2);
total_num = num_random+num_ortho;

norm_xx  = normalize_2nd(xx);
Cut_all = [];

cluster        = zeros(n,1);
Cut_min        = intmax;
cluster_centers = []; 


for i=1:total_num
    
    %assign the center index ortho or randomly
    if (i <= num_ortho)
        center_index = gen_orthogonal_centers(norm_xx(1:k,:));
    else
        is_ran_center      = zeros(1,n+1);
        is_ran_center(n+1) = 1;
        center_index       = zeros(1,k);
        for j=1:k
            ran_index=n+1;
            while is_ran_center(ran_index)
                ran_index=floor(1+n*rand(1));
            end
            center_index(j)=ran_index;
        end
    end
    % run the kmeans on them
    kcenter    = xx(:,center_index);
    %  unnormalized    
    [x_kmeans,cluster_centers_curr] = kmeans(xx',k,'Start',kcenter');
    
    % Ratio Cut, Ratio Cheeger Cut
    % Normalized Cut, Normalized Cheeger Cut
    [Cut_kmeans,CCut_kmeans] = computeRCutValue(x_kmeans,W,normalized);
            
    if Cut_kmeans < Cut_min
        Cut_min =  Cut_kmeans;
        cluster  = x_kmeans;
        cluster_centers = cluster_centers_curr; 
    end
            
    Cut_all = [Cut_all, Cut_kmeans];
end

[Cut,idx.p_minCut]  = min(Cut_all);



end
