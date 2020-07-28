function [U_new, T_new, obj_fcn] = step_simpfcml(data, U, T,...
    cluster_n, m, eta, a, b, N_k, H, img_size)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

beta = 0.2; % depends on the noise
mf = U .^ m;
tf = T.^ eta;
tfo = (1-T) .^ eta;

if img_size(1) == 0 && img_size(2) == 0
    tau = simMeasure(data, cluster_n, U, T, m, eta, a, b, H);
    f = labelInfo(data, cluster_n, U, N_k, beta);
else
    tau = simMeasure(data, cluster_n, U, T, m, eta, a, b, H, img_size);
    f = labelInfo(data, cluster_n, U, N_k, beta, img_size);
end

tmp = (f .* tau) .^ (2 / (m - 1));
den = repmat(sum(tmp, 1), cluster_n, 1);
U_new = tmp ./ den;

gamma = 0.08; % a user defined constant
f2 = f .^ 2;
tau2 = tau .^ 2;
den2 = ((b/gamma) ./ (f2 .* tau2)).^(1 / (eta - 1));
T_new = 1 ./ ( 1 + den2);

obj_fcn = sum(sum((1./f2 + 1./tau2).*(a*mf + b*tf))) + sum(gamma*sum(tfo));


end

