function [center, U, T, obj_fcn] = sim_pfcm_l(data, cluster_n, options)
%sim_pfcm_l Dataset clustering using similarity measure-based possibilistic
% FCM with label information
%   [CENTER, U, T, OBJ_FCN] = SIM_PFCM_L(DATA, CLUSTER_N) finds CLUSTER_N 
%   number of clusters in the dataset DATA. DATA is size M-by-N, where M is
%   the number of data points and N is the number of coordinates for each
%   data point. The coordinates for each cluster center are returned in the
%   rows of the matrix CENTER. The membership function matrix U contains
%   the grade of membership of each data point in each cluster. The values
%   0 and 1 indicate no membership and full membership respectively. Grades
%   between 0 and 1 indicate that the data point has partial membership in
%   a cluster. The typicality matrix T is contains the typicality of each
%   data point in each cluster and indicates each data point's intraclass
%   similarity and interclass dissimilarity. At each iteration, an 
%   objective function is minimized to find the best location for the
%   clusters and its values are returned in OBJ_FCN.
%
%   [CENTER, U, T, OBJ_FCN] = SIM_PFCM_L(DATA, CLUSTER_N, OPTIONS)
%   specifies a vector of options for the clustering process:
%       OPTIONS(1): exponent for the matrix U -> m      (default: 2.0)
%       OPTIONS(2): exponent for the matrix T -> eta    (default: 2.0)
%       OPTIONS(3): window size N_k                     (default: 3)
%       OPTIONS(4): parameter a                         (default: 0.65)
%       OPTIONS(5): parameter b                         (default: 0.35)
%       OPTIONS(6): epsilon in |U(t+1)-U(t)| < epsilon  (default: 1e-5)
%       OPTIONS(7): maximum number of iterations        (default: 100)
%       OPTIONS(8): info display during iteration       (default: 1)
if nargin ~= 2 && nargin ~= 3
	error(message("Incorrect Number of Input Arguments"))
end

default_options = [2;	% exponent for the partition matrix U
        2;	% exponent for the typicality matrix T
        3;  % window size N_k
        0.65;   % parameter a
        0.35;   % parameter b
        1e-5;	% min. amount of improvement
		100;	% max. number of iteration
		1];	% info display during iteration 

if nargin == 2
	options = default_options;
else
	% If "options" is not fully specified, pad it with default values.
	if length(options) < 8
		tmp = default_options;
		tmp(1:length(options)) = options;
		options = tmp;
	end
	% If some entries of "options" are nan's, replace them with defaults.
	nan_index = find(isnan(options));
	options(nan_index) = default_options(nan_index);
	if options(1) <= 1 || options(2) <= 1
		error(message("Fuzziness Parameter Should Be Greater Than 1"))
	end
end

m = options(1);         % exponent for U
eta = options(2);       % exponent for T
N_k = options(3);       % window size
a = options(4);         
b = options(5);
e = options(6);         % stop criteria epsilon
max_iter = options(7);  % maximum number of iterations
display = options(8);   % display info or not

obj_fcn = zeros(max_iter, 1);   % array for objective function

U = initU(data, cluster_n);
T = initT(data, cluster_n);

for i = 1 : max_iter
    pre_U = U;
    [U, T, center, obj_fcn(i)] = step_simpfcml(data, U, T, cluster_n, ...
        m, eta, a, b, N_k);
    if display
        fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
    end
    % check termination condition
    if i > 1
		if norm((U-pre_U),'fro') < e, break; end
    end
end

iter_n = i;	% Actual number of iterations 
obj_fcn(iter_n+1:max_iter) = [];
end

