function [U_new, T_new, center_new, obj_fcn] = step_simpfcml(data, U, T,...
    cluster_n, m, eta, a, b, N_k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mf = U .^ m;
tf = T.^ eta;
tfo = (1-T) .^ eta;

end

