function [U_new, T_new, obj_fcn] = step_simpfcml(data, U, T,...
    cluster_n, m, eta, a, b, N_k, H, img_size)
%STEP_SIMPFCML One step in clustering algorithm
%   Perform one iteration of the algorithm proposed in "Similarity 
%   Measure-Based Possibilistic FCM With Label Information for Brain MRI
%   Segmentation" where
%   DATA: matrix of data to be clustered. (Each row is a data point.)
%   U: partition matrix. 
%   T: typicality matrix.
%   CLUSTER_N: number of clusters.
%   M: exponent (> 1) for the partition matrix.
%   ETA: exponent (> 1) for the typicality matrix.
%   A: parameter a
%   B: parameter b
%   N_K: neighborhood window size in measuring local label information
%   H: bandwidth
%   IMG_SIZE: input image size (2D) 
%   U_NEW: new partition matrix.
%   T_NEW: new typicality matrix.
%   OBJ_FCN: objective function for partition U.

beta = 0.2; % depends on the noise
mf = U .^ m;
tf = T.^ eta;
tfo = (1-T) .^ eta;

if img_size(1) == 0 && img_size(2) == 0
    tau = simMeasure(data, cluster_n, U, T, m, eta, a, b, H);
    % f = labelInfo(data, cluster_n, U, N_k, beta);
    f = ones(size(U));  % if input data is a vector of 2D points we don't use local label info.
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

obj_fcn = sum(sum(((a*mf + b*tf)./(f2 .* tau2)))) + sum(gamma*sum(tfo));


end

