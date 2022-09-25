function [FinalX, fv, gfv, gfgf0, funs, grads] = ROPTLib_Grassmann(X,W,p,maxIter,normalized)
% Interface between ROPTLIB library and the clustering algorithm pGrass
% 
% "Multiway p-spectral graph cuts on Grassmann manifolds",
% by Pasadakis, D., Alappat, C.L., Schenk, O., Wellein, G,
%
% published in Machine Learning 111, 791–829 (2022).

k = size(X.main,2);
n = size(X.main,1);

constraint_first = false;

% Cost function evaluation
fhandle  =  @(X) pGrass_Functional(X,W,p,normalized);


% Euclidean gradient evaluation
ghandle  = @(X) pGrass_Gradient(X,W,p,constraint_first,normalized);
hhandle  = @(X,eta) pGrass_Hessian(X,eta,W,p);

SolverParams.IsCheckParams = 1;
SolverParams.IsCheckGradHess = 0;
SolverParams.method = 'RNewton';

SolverParams.DEBUG = 3;
if p == 1.1
    SolverParams.Max_Iteration = maxIter;
else
    SolverParams.Max_Iteration = maxIter;
end
SolverParams.OutputGap = 10;
SolverParams.LineSearch_LS = 0; % Armijo
SolverParams.Minstepsize = 1e-10;
SolverParams.RCGmethod = 2; %H-S


ManiParams.name = 'Grassmann';
ManiParams.n = n;
ManiParams.p = k;
ManiParams.ParamSet = 1;
HasHHR = 0;

initialX.main = X.main;


[FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times] = ...
    DriverOPT(fhandle, ghandle, hhandle, SolverParams, ManiParams, HasHHR, initialX);

end


function FLAG = IsStopped(x, gf, f, ngf, ngf0)

global objold;
global NoFuncChange;
global vec_check;
global x_old;
global NoSolChange;
global Adjacency;
global true_sol;
global GLB_NoRCutChange;
global GLB_ACC;
global GLB_ACC_orth;
global GLB_RCut;
global GLB_RCut_orth;
global GLB_NumClusters;

FLAG = 0;
eigennorm_diff = norm(x.main - x_old);
fprintf('---------------------------------------\n');
fprintf('eigennorm_diff: %.15f\n',eigennorm_diff);

if eigennorm_diff < 1e-6
    NoSolChange = NoSolChange + 1;
else
    NoSolChange = 0;
end

if(NoSolChange == 2)
    fprintf('---------------------------------------\n');
    fprintf('.      Eigenvectors converged         .\n');
    fprintf('---------------------------------------\n');
    FLAG = 1;
    NoSolChange = 0;
end

if (vec_check)
    functChange = 1-x.funcarray./objold;
else
    functChange=1-x.objective/objold;
end

if (max(abs(functChange))< 1e-8)
    NoFuncChange = NoFuncChange + 1;
else
    NoFuncChange = 0;
end

if (NoFuncChange==10)
    fprintf('---------------------------------------\n');
    fprintf('.         Objective converged         .\n');
    fprintf('---------------------------------------\n');
    FLAG = 1;
    NoFuncChange = 0;
end

% Check if converged
if (ngf / sqrt(size(gf,1)) < 1e-6)
    fprintf('---------------------------------------\n');
    fprintf('.          Gradient converged         .\n');
    fprintf('---------------------------------------\n');
    FLAG=1;
end

if (vec_check)
    objold = x.funcarray;
else
    objold = x.objective;
end


W = Adjacency;
purpose = 2;
[Confusion_p,ACC_p,RCut_p,~,x_p,ACC_p_orth,RCut_p_orth] = Quality_Clustering(x.main,GLB_NumClusters,W,true_sol,purpose);
% [Confusion_p,ACC_p,RCut_p,~,x_p,ACC_p_orth,RCut_p_orth] = Quality_Clustering(x.main,GLB_NumClusters,W,true_sol,purpose);


% ACC  = [ACC;ACC_p];
% ACC_orth  = [ACC_orth;ACC_p_orth];
GLB_RCut       = [GLB_RCut;RCut_p];
GLB_RCut_orth  = [GLB_RCut_orth;RCut_p_orth];

index = size(GLB_RCut_orth,1);
if(index > 1)
    RCutChange = GLB_RCut_orth(index) - GLB_RCut_orth(index-1);
    fprintf('==================\n');
    fprintf('RCutChange = %f\n',RCutChange);
else
    RCutChange = 1000;
end

if (max(abs(RCutChange))< 1e-3)
    GLB_NoRCutChange = GLB_NoRCutChange + 1;
else
    GLB_NoRCutChange = 0;
end

if (GLB_NoRCutChange==5)
    fprintf('---------------------------------------\n');
    fprintf('.             RCut converged          .\n');
    fprintf('---------------------------------------\n');
    FLAG = 1;
    GLB_NoRCutChange = 0;
end

x_old = x.main;

end
