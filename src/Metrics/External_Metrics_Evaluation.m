function [ACC,Confusion,VI,CE,NMI,F_score] = External_Metrics_Evaluation(clusters,labels,method)
% Input: 
% clusters
% labels
% method:1 most frequent
% method = 2; % munkres




[inferred_labels,~] = label_data(clusters,labels,method);


%% Confusion Matrix
[Confusion] = confusionmat(labels,inferred_labels);

%% ACC
n = length(clusters);
diff = (labels - inferred_labels);
hits = length(find(diff==0));
ACC  = hits/n;

%% VI 
[VI] = compare_clusterings(labels,clusters);

%% CE
[CE,~] = clustering_error(labels, clusters);
 
%% NMI
[NMI] = nmi(labels, clusters);
%% f-score
[Scores] = evaluate_scores(labels,inferred_labels);
F_score = Scores(6);


end
