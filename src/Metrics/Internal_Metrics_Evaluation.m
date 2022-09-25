function [RCut,RCCut,Modul,Dunn,NCut,NCCut] = Internal_Metrics_Evaluation(clusters,W,cluster_centers,V)

%% Dunn calculation
Dunn = -1;
if nargin > 2
    V = V';
    [Dunn] = dunn_index_calculation(V,clusters,cluster_centers);
end

%% Modularity calculation
Modul  = QFModul(clusters,W);

%% Ratio Cut, Ratio Cheeger Cut
normalized   = false;
[RCut,RCCut] = computeRCutValue(clusters,W,normalized);

%% Normalized Cut
normalized   = true;
[NCut,NCCut] = computeRCutValue(clusters,W,normalized);

end
