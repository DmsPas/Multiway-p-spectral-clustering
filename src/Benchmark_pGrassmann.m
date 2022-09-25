% Main execution script for the paper
%
% "Multiway p-spectral graph cuts on Grassmann manifolds",
% by Pasadakis, D., Alappat, C.L., Schenk, O., Wellein, G,
%
% published in Machine Learning 111, 791–829 (2022).
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%% Add paths and rng
clear all;close all;
rng(1991);
addpath_pGrass;
%% Declare global variables
global func_counter;
global grad_counter;
global objold;
global NoFuncChange;
global ctr_in_hess;
global Hess_store;
global vec_check;
global x_old;
global NoSolChange;
global Adjacency;
global true_sol;
global GLB_NumClusters;
global GLB_ACC;
global GLB_ACC_orth;
global GLB_RCut;
global GLB_RCut_orth;
global GLB_NoRCutChange;
%% Output to file
write_output_to_file = false;
if write_output_to_file
    diary on;
    diary('pGrass_output.txt');
end

%% Define minimization (or max) target
Targets = {'Graph Cut'};

%% Strings Initialization

StringRCut     = strings([3,1]);
StringACC      = strings([3,1]);
StringACC_max  = strings([3,1]);
StringVI       = strings([3,1]);
StringCE       = strings([3,1]);
StringNMI      = strings([3,1]);
StringModul    = strings([3,1]);
StringDunn     = strings([3,1]);
StringTime     = strings([3,1]);
StringFscore   = strings([3,1]);
Stringpbest    = strings([3,1]);
StringRCCut    = strings([3,1]);
StringNCut     = strings([3,1]);
StringNCCut    = strings([3,1]);
%%
cases = {

    % --------------------------- %
    %         LFR Data           %
    % --------------------------- %
    %     'LFR_10';
    %     'LFR_12';
    %     'LFR_14';
    %     'LFR_16';
    %     'LFR_18';
    %     'LFR_20';
    %     'LFR_22';
    %     'LFR_24';
    %     'LFR_26';
    %     'LFR_28';
    %     'LFR_30';
    %     'LFR_32';
    %     'LFR_34';
    %     'LFR_36';
    %     'LFR_38';
    %     'LFR_40';

    % -------------------------------------------- %
    %    Gauss synthetic cases at Input/Graphs     %
    % -------------------------------------------- %
    % ------------------------------------------- %
    %         Omniglot digits at Input/Graphs     %
    % ------------------------------------------- %
    'Omniglot_AngloSaxon_10NN.mat';

    };
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n');
fprintf('|                                                    |\n');
fprintf('|     p-spectral Clustering on Riemannian manifolds  |\n');
fprintf('|                                                    |\n');
fprintf('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n\n');
% Parameters
p_final      = 1.100; % final value of p
factor       = .9;    % factor of p-reduction (only if use_p_contin == 0)
normalized   = 1;     % 1/0 normalized/ unnormalized p-clustering
use_p_contin = 1;     % 1/0 use pseudo-continuation for the reduction of p


fprintf('Normalized: %d \n',normalized);
fprintf('Use p continuation: %d \n',use_p_contin);
nc = length(cases);
maxlen = 0;
for c = 1:nc
    if length(cases{c}) > maxlen
        maxlen = length(cases{c});
    end
end


num_cases = 0;
for c = 1:nc
    fprintf('.');
    temp = load(cases{c});
    if ~isfield(temp,'coords')
        temp.coords = [];
    end
    Graph_n_labels(c) = temp;
    num_cases = num_cases + 1;
end


node_cluster_ratio = zeros(nc,1);

fprintf('\n\n Report Cases - %6d %12s %9s %10s %20s\n',num_cases,'Nodes','Edges','Clusters','Clusters/Nodes');
fprintf('------------------------------------------------------\n');
for c = 1:nc
    % get number of vertices
    params.numberOfVertices = size(Graph_n_labels(c).W,1);
    % get number of edges
    params.numberOfEdges    = nnz(Graph_n_labels(c).W)/2;
    % get number of clusters
    params.numberOfClusters = size(unique(Graph_n_labels(c).label),1);
    % calculate the cl/node ratio
    node_cluster_ratio(c) = params.numberOfClusters/params.numberOfVertices;
end

[~,ratio_id] = sort(node_cluster_ratio);

for orig_ID = 1:nc
    % reorder according to cl/node ratio
    c = ratio_id(orig_ID);
    % cancel the reordering
    c = orig_ID;
    % get number of vertices
    params.numberOfVertices = size(Graph_n_labels(c).W,1);
    % get number of edges
    params.numberOfEdges    = nnz(Graph_n_labels(c).W)/2;
    % get number of clusters
    params.numberOfClusters = size(unique(Graph_n_labels(c).label),1);
    % spacing
    spacers = repmat('.', 1, maxlen+3-length(cases{c}));
    % print headers
    fprintf('%s %s %10d %10d %10d %20.3f\n',cases{c},spacers,params.numberOfVertices, ...
        params.numberOfEdges,params.numberOfClusters,node_cluster_ratio(c));
end

fprintf('------------------------------------------------------\n\n');

%% Loop over cases
for orig_ID = 1:nc
    c = orig_ID;
    spacers = repmat('.', 1, maxlen+3-length(cases{c}));
    W                = Graph_n_labels(c).W;
    labels           = Graph_n_labels(c).label;
    if min(labels) == 0
        labels = labels + 1;
    elseif min(labels) == 1
        labels = labels;
    else
        fprintf('Labelling error, min(labels) = %f\n',min(labels));
    end

    GLB_NumClusters  = size(unique(Graph_n_labels(c).label),1);
    Num_Eigs         = GLB_NumClusters; % # of eigenvecs to be considered
    GLB_NoRCutChange = 0;
    GLB_ACC          = [];
    GLB_ACC_orth     = [];
    GLB_RCut         = [];
    GLB_RCut_orth    = [];
    Adjacency        = W;
    objold           = -1;
    ctr_in_hess      = 0;
    NoFuncChange     = 0;
    func_counter     = 0;
    grad_counter     = 0;
    NoSolChange      = 0;
    Hess_store       = cell(Num_Eigs,1);
    true_sol         = labels;
    %% spectral clustering in the 2-norm
    [L,Diag,vw] = CreateLapl(W,normalized);
    [V,lambda] = eigs(L,Num_Eigs,'SM');
    if normalized == 1

        % Transform eigs of L_sym to those of L_rw, by multiplying the
        % entries of the eigenvectors by sqrt(d). 
        % Afterwards renorm them s.t. norm = 1.
        for i=1:Num_Eigs
            V(:,i) = V(:,i)./sqrt(vw);
            V(:,i) = V(:,i)/norm(V(:,i));
        end
    end

    % Initialize x_old with the p=2 eigenvecs
    x_old      = V(:,1:GLB_NumClusters);
    V_spec     = V;

    %% p-spectral clustering on Grassmann manifold
    tStart_pGrass = tic;    
    [V_Gr,Levels_Gr,clusters_Grass] =  Kernel_Grass_Clustering(W,V_spec,p_final,factor,labels,normalized, use_p_contin); labels = labels';
    
    clusters_Grass = clusters_Grass';
    Metrics_pGrass.time   = toc(tStart_pGrass);
    %  Internal pGrassmann Metrics Evaluation
    [Metrics_pGrass.RCut,Metrics_pGrass.RCCut,Metrics_pGrass.Modul,Metrics_pGrass.Dunn,Metrics_pGrass.NCut,Metrics_pGrass.NCCut]...
        = Internal_Metrics_Evaluation(clusters_Grass,W);

    %  External pGrassmann Metrics
    method = 1; % most frequent
    [Metrics_pGrass.ACC,Metrics_pGrass.Confusion,Metrics_pGrass.VI,Metrics_pGrass.CE,Metrics_pGrass.NMI,Metrics_pGrass.Fscore] = External_Metrics_Evaluation(clusters_Grass,labels,method);
         
    for target_id = 1:length(Targets)
        Target = Targets{target_id};

        fprintf('%s %s', cases{c}, spacers);

        [Stringpbest(target_id)]     = String_Metrics_pGrass(Stringpbest(target_id),cases{c},Levels_Gr.p_best);

        [StringRCut(target_id)]      = String_Metrics_pGrass(StringRCut(target_id),cases{c},Metrics_pGrass.RCut);

        [StringRCCut(target_id)]     = String_Metrics_pGrass(StringRCCut(target_id),cases{c},Metrics_pGrass.RCCut);

        [StringNCut(target_id)]      = String_Metrics_pGrass(StringNCut(target_id),cases{c},Metrics_pGrass.NCut);

        [StringNCCut(target_id)]     = String_Metrics_pGrass(StringNCCut(target_id),cases{c},Metrics_pGrass.NCCut);

        [StringModul(target_id)]     = String_Metrics_pGrass(StringModul(target_id),cases{c},Metrics_pGrass.Modul);

        [StringDunn(target_id)]      = String_Metrics_pGrass(StringDunn(target_id),cases{c},Metrics_pGrass.Dunn);

        [StringACC(target_id)]       = String_Metrics_pGrass(StringACC(target_id),cases{c},Metrics_pGrass.ACC);

        [StringVI(target_id)]        = String_Metrics_pGrass(StringVI(target_id),cases{c},Metrics_pGrass.VI);

        [StringCE(target_id)]        = String_Metrics_pGrass(StringCE(target_id),cases{c},Metrics_pGrass.CE);

        [StringNMI(target_id)]       = String_Metrics_pGrass(StringNMI(target_id),cases{c},Metrics_pGrass.NMI);

        [StringFscore(target_id)]    = String_Metrics_pGrass(StringFscore(target_id),cases{c},Metrics_pGrass.Fscore);

        [StringTime(target_id)]      = String_Metrics_pGrass(StringTime(target_id),cases{c},Metrics_pGrass.time);


    end
end

Header1 = sprintf('%30s %10s %10s %10s', 'Case', 'pGrass');

for target_id = 1:length(Targets)
    Target = Targets{target_id};

    fprintf('\n=========================\n');
    fprintf('\n@@@ Target: %s @@@\n', Target);
    fprintf('\n=========================\n');
    fprintf('\n@@@ INTERNAL METRICS @@@\n');

    fprintf('\n-------\n');
    fprintf('pBest \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(Stringpbest(target_id));

    fprintf('\n-------\n');
    fprintf('RCut (smaller the better) \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringRCut(target_id));

    fprintf('\n-------\n');
    fprintf('RCCut (smaller the better) \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringRCCut(target_id));

    fprintf('\n-------\n');
    fprintf('NCut (smaller the better) \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringNCut(target_id));

    fprintf('\n-------\n');
    fprintf('NCCut (smaller the better) \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringNCCut(target_id));

    fprintf('\n-------\n');
    fprintf('Modul (bigger the better) \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringModul(target_id));

    fprintf('\n@@@ EXTERNAL METRICS @@@\n');

    fprintf('-------\n');
    fprintf('ACC (bigger the better) \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringACC(target_id));

    fprintf('-------\n');
    fprintf('VI (smaller the better) \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringVI(target_id));

    fprintf('-------\n');
    fprintf('CE (smaller the better) \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringCE(target_id));

    fprintf('-------\n');
    fprintf('NMI (bigger the better) \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringNMI(target_id));

    fprintf('-------\n');
    fprintf('Time \n');
    fprintf('-------\n\n');
    fprintf('%s \n',Header1);
    fprintf(StringTime(target_id));

end
