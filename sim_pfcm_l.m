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


end

