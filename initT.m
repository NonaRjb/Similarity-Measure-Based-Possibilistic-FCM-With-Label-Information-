function T = initT(data_n, cluster_n)
%INITT Initializes typicality matrix with random numbers
T = rand(cluster_n, data_n);
col_sum = sum(T);
T = T./col_sum(ones(cluster_n, 1), :);
end

