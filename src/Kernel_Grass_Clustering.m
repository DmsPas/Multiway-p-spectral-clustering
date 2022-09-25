function [V,Levels,clusters] =  Kernel_Grass_Clustering(W,V,p_final,factor,label,normalized,use_p_contin)
% Kernel of the paper
% 
% "Multiway p-spectral graph cuts on Grassmann manifolds",
% by Pasadakis, D., Alappat, C.L., Schenk, O., Wellein, G,
%
% published in Machine Learning 111, 791–829 (2022).
% 
% Input:
%       W: nxn adjacency matrix.
%       V: starting nxk eigenvectors (usually at p=2) with trad. spectral
%          clustering.
%       p_final: final level of the p-norm where the graph will be
%       clustered.
%       factor: if use use_p_contin == 0 the factor with which the value of
%       p is reduced.
%       label: the label assignments (true solution).
%       normalized: 1/0 for normalized/unnormalized multiway p-spectral
%       clustering.
%       use_p_contin: 1/0 decrease the level of p in a pseudo continuous
%       way.
% Output:
%       V: the nxk p-eigenvectors at the last p-level
%       Levels: Information regarding the progress of the continuous
%               (obj, grad) and discrete (Cut, ACC, NMI) metrics for the
%               various p-levels.
%       clusters: nx1 vector of cluster assignments.
% 



%%
tic;
% Initializations
maxIter     = 20;
V_spec.main = V;
p_old       = 2;
% Metrics for intermediate levels
Levels.ACC       = [];
Levels.NMI       = [];
Levels.RCut      = [];
Levels.Modul      = [];
Levels.p_monotone = [];
Levels.RCut_monotone = [];
Levels.p         = 2;
Levels.obj       = [];
Levels.grad_norm = 1;
Levels.obj_all     = zeros(20,20);
Levels.grads_all   = zeros(20,20);
Clusters_all       = {};
K                  = size(V,2);
%

[obj_curr] = pGrass_Functional(V_spec,W,p_old,normalized);
Levels.obj = [Levels.obj;obj_curr];



%     EVALUATIONS at  p = 2  %
% -------------------------------------- %
p_levels = 1;
num_ortho  = 50;
num_random = 100; 

[RCut_init,clusters_init]  = kmeans_orthogonal(V(:,1:K)',K, num_ortho, num_random,label,W,normalized);
Levels.RCut                = [Levels.RCut; RCut_init];
Clusters_all{p_levels}     = clusters_init;
% -------------------------------------- %

if p_final<2

    if use_p_contin == 0
        p_new = p_old * factor;
    else
         [p_new] = p_continuation(p_old,1e-2);
    end
    %% main loop
    while p_new>p_final
        
        p_levels = p_levels+1;
                
        if (p_new*factor<p_final)
            p_new = p_old*sqrt(p_final/p_old);
        end
        
        B = V;
        V = [];
        V.main = B;
        
        fprintf('====================\n');
        fprintf('Current p-value: %f \n',p_new);
        fprintf('====================\n');
        
        [V, fv, gfv, gfgf0, funs, grads] = ROPTLib_Grassmann(V,W,p_new,maxIter,normalized);
        Levels.obj         = [Levels.obj;fv];
        Levels.grad_norm   = [Levels.grad_norm;gfv];
        
        Levels.obj_all(1:length(funs),p_levels)   = funs;
        Levels.grads_all(1:length(funs),p_levels) = grads;
                

        V = V.main;
        %       INTERMEDIATE EVALUATIONS         %
        % -------------------------------------- %        
        [RCut_inter,clusters]  = kmeans_orthogonal(V(:,1:K)', K, num_ortho, num_random,label,W,normalized);               
        Levels.RCut            = [Levels.RCut; RCut_inter];
        Levels.p               = [Levels.p;p_new];
        Clusters_all{p_levels} = clusters;
        % -------------------------------------- %
        
        
        % Update the value of p
        p_old    = p_new;        
        [p_new] = p_continuation(p_old,1e-2);

    end
    
    p_levels = p_levels+1;
    fprintf('====================\n')
    fprintf('Final p-value: %f \n',p_final)
    fprintf('====================\n')
    B = V;
    V = [];
    V.main = B;
    
    %% Evaluation directly at p = p_final
    [V, fv, gfv, gfgf0,funs, grads] = ROPTLib_Grassmann(V,W,p_final,maxIter,normalized);
    Levels.obj         = [Levels.obj;fv];
    Levels.grad_norm   = [Levels.grad_norm;gfv];
    Levels.obj_all(1:length(funs),p_levels)   = funs;
    Levels.grads_all(1:length(funs),p_levels) = grads;
        
        
    V      = V.main;
    
    %       FINAL  EVALUATIONS, p = p_final  %
    % -------------------------------------- %    
    [RCut_inter,clusters]  = kmeans_orthogonal(V(:,1:K)', K, num_ortho, num_random,label,W,normalized);               
    Levels.RCut            = [Levels.RCut; RCut_inter];
    Levels.p               = [Levels.p;p_final];
    Clusters_all{p_levels} = clusters;
    % -------------------------------------- %
end

% Find minimum RCut values
[min_val,min_idx]    = min(Levels.RCut);
Levels.p_best        = Levels.p(min_idx);
clusters             = Clusters_all{min_idx};

temp           = [Levels.RCut(1); Levels.RCut(1:end-1)];
monotone_ratio = Levels.RCut./temp;
[val_p_monotone,idx_p_monotone] = find(monotone_ratio>1.05);
if isempty(idx_p_monotone)
    idx_p_monotone = min_idx;
end

Levels.p_monotone    = Levels.p(idx_p_monotone);
Levels.RCut_monotone = Levels.RCut(idx_p_monotone);

end
