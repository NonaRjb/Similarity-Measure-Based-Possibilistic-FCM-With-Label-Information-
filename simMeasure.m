function S = simMeasure(data, cluster_n, a, b, U, T, img_size)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

end

function S = similarity(data, U, T, m, eta, a, b, img_size)
end

function rho = rel_sim(data, U, T, m , eta, a, b, img_size)
if nargin < 2
    W = point_sim(data);
else
    W = point_sim(data, img_size);
end
den = sum(a * U .^ m + b * T .^ eta, 2);
rho = (a * U .^ m + b * T .^ eta) * W ./ den; % todo: d_lk<w is not considered
end

function W = point_sim(data, img_size)
if nargin < 2
    dis = coordinate_dist(data);
else
    dis = coordinate_dist(data, img_size);
end
d = intensity_dist(data);
W = exp(-1*d.^2 - dis);
end

function d = intensity_dist(data)
%INTENSITY_DIST measures the Euclidean distance between 2 data points
d = zeros(size(data, 1), size(data, 1));
for l = 1 : size(data, 1)
    for k = l : size(data, 1)
        d(l, k) = norm(data(l) - data(k));
    end
end
d = d + d'; % todo: still not sure about what H is in the paper (should be devided by H)
end

function dis = coordinate_dist(data, img_size)
if nargin == 2
    data_ind = 1 : length(data);
    data_new(:, 1:2) = ind2sub(img_size, data_ind);
else
    data_new = data;
end
dis = zeros(size(data, 1), size(data, 1));
for l = 1 : size(data, 1)
    for k = l : size(data, 1)
        dis(l, k) = norm(data_new(l) - data_new(k));
    end
end
diag_data = max(dis(:));
dis = dis + dis';
dis = dis ./ diag_data;
end
