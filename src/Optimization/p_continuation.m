function [p] = p_continuation(p_val,tol)
% IPOPT inspired pseudo-continuation function for the decrease of p
% Input
% p_val, tol: Value of p and tol are needed
% Output
% p: the level of p
% 
% pasadd@usi.ch


kappa    = 0.9;
theta    = 1.25;
p        = p_val;     

p        = 1 + max(tol,min(kappa*(p-1),(p-1)^theta));


end
