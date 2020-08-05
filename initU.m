function U = initU(data, cluster_n)
%initU Initializes membership matrix, U, using the results of FCM
%   U = initU(DATA, CLUSTER_N) finds CLUSTER_N number of clusters in the
%   data set DATA. DATA is size M-by-N, where M is the number of data
%   points and N is the number of coordinates for each data point. The
%   The membership function matrix U contains the grade of membership of
%   each DATA point in each cluster. The values 0 and 1 indicate no
%   membership and full membership respectively. Grades between 0 and 1
%   indicate that the data point has partial membership in a cluster.
options = [NaN NaN NaN false];
[~, U] = fcm(data, cluster_n, options);
[~, i] = max(U, [], 1);
% U = rand(cluster_n, length(data));
% col_sum = sum(U);
% U = U./col_sum(ones(cluster_n, 1), :);
% [~, i] = max(U, [], 1);
figure
% gscatter(data(:,1), data(:,2), i')
imshow(reshape(i, 128, 128), [])
end

