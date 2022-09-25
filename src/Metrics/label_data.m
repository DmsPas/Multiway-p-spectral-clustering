function [x_inferred, inferred_label] = label_data(x,labels,method)

% Number of clusters
K = max(x);

L = size(unique(labels),1);

inferred_label = -1*ones(K,1);

x_inferred = x;


if method == 1
    
    for k=1:K
        index             = find(x==k);
        cur_labels        = labels(index);
        actual_label      = mode(cur_labels);
        inferred_label(k) = actual_label;
        x_inferred(index) = actual_label;
    end
    
%         if size(unique(inferred_label),1) ~= size(inferred_label,1)
%             fprintf('Clusters are merging into 1 label\n');
%             fprintf('Probably unbalanced dataset\n');
%         end
    
    
else
    
    x_labeled     = x;       
    confusion_mat = confusionmat(labels, x_labeled);
    Cost_mat      = zeros(L,K);
    
    for i = 1:L
        for j = 1:K
            %             delta = confusion_mat(i,:)-confusion_mat(i,j);
            Cost_mat(i,j) = sum(confusion_mat(i,:))-confusion_mat(i,j);
        end
    end
    costofnonassignment = 2*max(max(Cost_mat));
%     costofnonassignment = .2;
    [assignments, unassignedTracks, unassignedDetections] = ...
        assignmunkres(Cost_mat,costofnonassignment);
    
    
    %% Munkres vectorized implementation
    %     [assignments_Matrix,Cost_assgn]    = munkres(Cost_mat);
    %     assignments = zeros(K,2);
    %     for i = 1:K
    %         assignments(i,1) = i;
    %         assignments(i,2) = find(assignments_Matrix(i,:));
    %     end
    %%
    
    inferred_label(assignments(:,2)) = assignments(:,1);
       
    for j = 1:L
        index             = find(x_labeled==j);
        x_inferred(index) = inferred_label(j);
    end
    
    if size(unique(inferred_label),1) ~= size(inferred_label,1)
        fprintf('Clusters are merging into 1 label\n');
        fprintf('Probably unbalanced dataset\n');
    end
    
    
end



end