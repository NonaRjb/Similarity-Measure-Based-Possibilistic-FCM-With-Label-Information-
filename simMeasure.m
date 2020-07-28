function tau = simMeasure(data, cluster_n, U, T, m, eta, a, b, img_size)
%SIMMEASURE measures the interclass and intraclass similarity
%   tau = simMeasure(data, cluster_n, U, T, m, eta, a, b, img_size)
%   measures interclass and intraclass similarity for each data point in
%   a cluster. 
%   tau(ik) = sum( rho(jk) * S(ij) ); j = 1 : C

if nargin < 9
    rho = rel_sim(data, U, T, m, eta, a, b);
    S = similarity(data, cluster_n, U, T, m, eta, a, b);
else
    rho = rel_sim(data, U, T, m, eta, a, b, img_size);
    S = similarity(data, cluster_n, U, T, m, eta, a, b, img_size);
end

tau = S * rho;

end

function S = similarity(data, cluster_n, U, T, m, eta, a, b, img_size)

if nargin < 9
    rho = rel_sim(data, U, T, m, eta, a, b, img_size);
else
    rho = rel_sim(data, U, T, m, eta, a, b);
end

max_vals = max(rho, [], 2);
indices = [];
% flag = zeros(1, size(rho, 1));
% pre_flag = 0;

for i = 1 : size(rho, 1)
    indices = [indices, find(rho(i, :) == max_vals(i))];
    % flag(1, i) = length(indices) + pre_flag;
    % pre_flag = flag(1, i);
end
points = unique(indices);

S = zeros(cluster_n, cluster_n);
for i = 1 : cluster_n
    for j = i : cluster_n
        num = 0;
        den1 = 0;
        den2 = 0;
        for k = points
            num = num + rho(i, k) * rho(j, k);
            den1 = den1 + rho(i, k);
            den2 = den2 + rho(j, k);
        end
        den = den1 * den2;
        den = sqrt(den);
        S(i, j) = num / den;
    end
end

end

function rho = rel_sim(data, U, T, m , eta, a, b, img_size)

if nargin < 8
    W = point_sim(data);
else
    W = point_sim(data, img_size);
end
den = sum(a * U .^ m + b * T .^ eta, 2);
rho = (a * U .^ m + b * T .^ eta) * W ./ den; % todo: d_lk < w is not considered
end

function w = w_thresh(data)
mean_data = mean(data, 1);
data_white = data - repmat(mean_data, size(data, 1), 1);
data_white = sqrt(sum(data_white .^ 2, 2));
d_bar = mean(data_white, 1);
delta = 3;
w = delta * d_bar;
end

function W = point_sim(data, img_size)
if nargin < 2
    dis = coordinate_dist(data);
else
    dis = coordinate_dist(data, img_size);
end
d = intensity_dist(data);
wth = w_thresh(data);
ind = d > wth;

if max(d(:)) > 0
    d = d / max(d(:)); % todo: still not sure about what H is in the paper (should be devided by H)
end

W = zeros(size(d));
W(ind) = exp(-1*d(ind).^2 - dis(ind));
end

function d = intensity_dist(data)
%INTENSITY_DIST measures the Euclidean distance between 2 data points
d = zeros(size(data, 1), size(data, 1));
for l = 1 : size(data, 1)
    for k = l : size(data, 1)
        d(l, k) = norm(data(l) - data(k));
    end
end
% moved to point_sim function
% if max(d(:)) > 0
%     d = d / max(d(:)); % todo: still not sure about what H is in the paper (should be devided by H)
% end
d = d + d'; 
end

function dis = coordinate_dist(data, img_size)
if nargin == 2
    data_ind = 1 : length(data);
    data_new = zeros(size(data,1), 2);
    [data_new(:, 1), data_new(:,2)] = ind2sub(img_size, data_ind);
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
