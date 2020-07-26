function f = labelInfo(data, cluster_n, U, N_k, beta, img_size)
%LABELINFO Computes a measure for local label information
%   f(ik) = 1 + beta * P(ik)  i:cluster_id and k:data point index
if nargin < 6
    P = neighbour_label(data, cluster_n, U, N_k);
else
    P = neighbour_label(data, cluster_n, U, N_k, img_size);
end
f = 1 + beta * P;
end

function P = neighbour_label(data, cluster_n, U, N_k, img_size)

[~, cluster_id] = max(U, [], 1);

if nargin < 5
    xy = max(data, [], 1);  % we assume that if data is not an image it is a set of 2D points
else
    xy = img_size;
end

if mod(N_k, 2) == 0
    error('The Window Size Should Be Odd')
end

cluster_2d = reshape(cluster_id, xy(1), xy(2));
Y = zeros(cluster_n, size(data, 1));
for i = 1 : cluster_num
    for k = 1 : size(data, 1)
       [row, col] = ind2sub(xy, k);
       v = row - (N_k-1)/2 : row + (N_k-1)/2;
       h = col - (N_k-1)/2 : col + (N_k-1)/2;
       v = v(v > 0);
       h = h(h > 0); 
       for r = v
           for c = h
               if cluster_2d(r, v) == i
                  Y(i, k) = Y(i, k) + 1 / sqrt((v-row)^2+(h-col)^2);
               end
           end
       end
    end
end
den0 = sum(Y, 1);
den = repmat(den0, cluster_n, 1);
P = Y ./ den;
end