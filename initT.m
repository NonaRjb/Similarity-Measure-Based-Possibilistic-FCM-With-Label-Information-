function T = initT(data_n, cluster_n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
T = rand(cluster_n, data_n);
col_sum = sum(T);
T = T./col_sum(ones(cluster_n, 1), :);
end

