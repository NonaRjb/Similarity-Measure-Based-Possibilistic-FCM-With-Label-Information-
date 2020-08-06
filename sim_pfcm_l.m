function [U, T, obj_fcn] = sim_pfcm_l(data, cluster_n, options)
%sim_pfcm_l Dataset clustering using similarity measure-based possibilistic
% FCM with label information
%   [U, T, OBJ_FCN] = SIM_PFCM_L(DATA, CLUSTER_N) finds CLUSTER_N 
%   number of clusters in the dataset DATA. DATA is size M-by-N, where M is
%   the number of data points and N is the number of coordinates for each
%   data point. The membership function matrix U contains
%   the grade of membership of each data point in each cluster. The values
%   0 and 1 indicate no membership and full membership respectively. Grades
%   between 0 and 1 indicate that the data point has partial membership in
%   a cluster. The typicality matrix T is contains the typicality of each
%   data point in each cluster and indicates each data point's intraclass
%   similarity and interclass dissimilarity. At each iteration, an 
%   objective function is minimized to find the best location for the
%   clusters and its values are returned in OBJ_FCN.
%
%   [U, T, OBJ_FCN] = SIM_PFCM_L(DATA, CLUSTER_N, OPTIONS)
%   specifies a vector of options for the clustering process:
%       OPTIONS(1): exponent for the matrix U -> m      (default: 2.0)
%       OPTIONS(2): exponent for the matrix T -> eta    (default: 2.0)
%       OPTIONS(3): window size N_k                     (default: 3)
%       OPTIONS(4): parameter a                         (default: 0.65)
%       OPTIONS(5): parameter b                         (default: 0.35)
%       OPTIONS(6): epsilon in |U(t+1)-U(t)| < epsilon  (default: 1e-5)
%       OPTIONS(7): maximum number of iterations        (default: 100)
%       OPTIONS(8): info display during iteration       (default: 1)
%       OPTIONS(9): image size [row col]                (default: [0 0])
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
		1;      % info display during iteration
        0;      % image size row
        0];     % image size col

if nargin == 2
	options = default_options;
else
	% If "options" is not fully specified, pad it with default values.
	if length(options) < 10
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
img_size = [options(9) options(10)];
H = h_bandwidth(data);

obj_fcn = zeros(max_iter, 1);   % array for objective function

U = initU(data, cluster_n);
T = initT(size(data,1), cluster_n);

[~, i] = max(U, [], 1);
figure
if img_size(1) > 0 && img_size(2) > 0 
    imshow(reshape(i, img_size(1), img_size(2)),[])
    title('FCM output')
else
    gscatter(data(:,1), data(:,2), i')
    title('FCM output')
end

for i = 1 : max_iter
    pre_U = U;
    [U, T, obj_fcn(i)] = step_simpfcml(data, U, T, cluster_n, ...
        m, eta, a, b, N_k, H, img_size);
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

