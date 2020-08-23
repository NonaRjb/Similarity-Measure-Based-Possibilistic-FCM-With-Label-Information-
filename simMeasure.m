function tau = simMeasure(data, cluster_n, U, T, m, eta, a, b, H, img_size)
%SIMMEASURE measures the interclass and intraclass similarity
%   tau = simMeasure(data, cluster_n, U, T, m, eta, a, b, img_size)
%   measures interclass and intraclass similarity for each data point in
%   a cluster. 
%   tau(ik) = sum( rho(jk) * S(ij) ); j = 1 : C

if nargin < 10
    rho = rel_sim(data, U, T, m, eta, a, b, H);
    S = similarity(data, cluster_n, U, T, m, eta, a, b, H);
else
    rho = rel_sim(data, U, T, m, eta, a, b, H, img_size);
    S = similarity(data, cluster_n, U, T, m, eta, a, b, H, img_size);
end

tau = S * rho;

end

function S = similarity(data, cluster_n, U, T, m, eta, a, b, H, img_size)
%SIMILARITY Computes the similarity measure S introduced in the paper
%   S(ij) = sum(rho(ik)*rho(jk))/sqrt(sum(rho(ik))*sum(rho(jk))) 
%   k = argmax{rho(ik) or rho(jk)}

if nargin < 9
    rho = rel_sim(data, U, T, m, eta, a, b, H, img_size);
else
    rho = rel_sim(data, U, T, m, eta, a, b, H);
end

% finding all of the points which maximize rho
max_vals = max(rho, [], 2);
indices = zeros(size(rho));

for i = 1 : size(rho, 1)
    indices(i, rho(i, :) == max_vals(i)) = 1;
end


S = zeros(cluster_n, cluster_n);
for i = 1 : cluster_n
    for j = i : cluster_n
        num = sum(rho(i, indices(i,:) == 1 | indices(j,:) == 1) .* ...
            rho(j, indices(i,:) == 1 | indices(j,:) == 1));
        den1 = sum(rho(i, indices(i,:) == 1 | indices(j,:) == 1));
        den2 = sum(rho(j, indices(i,:) == 1 | indices(j,:) == 1));
        den = den1 * den2;
        den = sqrt(den);
        S(i, j) = num / den;
    end
end

S = S + S' - diag(diag(S));

end

function rho = rel_sim(data, U, T, m , eta, a, b, H, img_size)
%REL_SIM Computes relative similarity measure, rho, introduced in the paper
%   rho(ik) = sum((a*u(il)^m+b*t(il)^eta)W(lk))/sum((a*u(il)^m+b*t(il))
%   W(lk) = exp(-(d(lk)/H)^2 - dis(lk)/Diag)
%   d(lk) = Euclidean distance between pixel intensities (in case of image
%   input) and Euclidean distance between point coordinates otherwise
%   dis(lk) = Euclidean distance between coordinates
%   Diag is the farthest distance of two points among
%   all points

if nargin == 9
    data_ind = 1 : length(data);
    data_new = zeros(size(data,1), 2);
    % extraction of 2D coordinate of each pixel
    [data_new(:, 1), data_new(:,2)] = ind2sub(img_size, data_ind);
    Diag = sqrt((img_size(1)-1)^2 + (img_size(2)-1)^2);
else
    data_new = data;
    Diag = H;
end
rho = zeros(size(U, 1), size(data, 1));
den = sum(a * U .^ m + b * T .^ eta, 2);
wth = w_thresh(data);

% according to the size of data choose which method to use for calculation
% of distances
if log2(size(data, 1)) <= 12
    d = sqrt(reshape(sum((repmat(data, size(data, 1), 1) - ...
        repelem(data, size(data, 1), 1)).^2, 2), size(data, 1), size(data, 1)));
    dis = sqrt(reshape(sum((repmat(data_new, size(data_new, 1), 1) - ...
        repelem(data_new, size(data_new, 1), 1)).^2, 2), size(data, 1), size(data, 1)));
    ind = d < wth;
    d = d / H;
    dis = dis / Diag;
    W = zeros(size(data, 1), size(data, 1));
    W(ind) = exp(-1*d(ind).^2 - dis(ind));
    clear d
    clear dis
    rho = (a * U .^ m + b * T .^ eta) * W ./ repmat(den, 1, size(data, 1));
    clear W
elseif 12 < log2(size(data, 1)) && log2(size(data, 1)) < 16
    if mod(size(data, 1), 4) == 0
        n0 = size(data, 1)/4;
        n = [n0, 2*n0, 3*n0, 4*n0];
    else
        n0 = floor(size(data, 1) / 4);
        n = [n0, 2*n0, 3*n0, size(data, 1)];
    end
    for i = 1 : length(n)
	if i ==  1
		d = sqrt(reshape(sum((repmat(data, n(i), 1) - ...
			repelem(data(1:n(i), :), size(data,1), 1)).^2, 2), size(data,1), n(i)));
		dis = sqrt(reshape(sum((repmat(data_new, n(i), 1) - ...
			repelem(data_new(1:n(i), :), size(data,1), 1)).^2, 2), size(data,1), n(i)));
		ind = d < wth;
		d = d / H;
		dis = dis / Diag;
		W = zeros(size(data, 1), n(i));
		W(ind) = exp(-1*d(ind).^2 - dis(ind));
		clear d
		clear dis
		rho(:, 1:n(i)) = (a * U .^ m + b * T .^ eta) * W ./ repmat(den, 1, n(i));
		clear W
	else
		d = sqrt(reshape(sum((repmat(data, n(i)-n(i-1), 1) - ...
        repelem(data(n(i-1)+1:n(i), :), size(data,1), 1)).^2, 2), size(data,1), n(i)-n(i-1)));
		dis = sqrt(reshape(sum((repmat(data_new, n(i)-n(i-1), 1) - ...
			repelem(data_new(n(i-1)+1:n(i), :), size(data,1), 1)).^2, 2), size(data,1), n(i)-n(i-1)));
		ind = d < wth;
		d = d / H;
		dis = dis / Diag;
		W = zeros(size(data, 1), n(i)-n(i-1));
		W(ind) = exp(-1*d(ind).^2 - dis(ind));
		clear d
		clear dis
		rho(:, n(i-1)+1:n(i)) = (a * U .^ m + b * T .^ eta) * W ./ repmat(den, 1, n(i)-n(i-1));
		clear W
	end
    end
else
	for k = 1 : size(data, 1)
	    d = sqrt(sum((repmat(data(k,:), size(data, 1), 1)- data).^2, 2));
	    dis = sqrt(sum((repmat(data_new(k,:), size(data_new, 1), 1)- data_new).^2, 2));
	    ind = d < wth;
	    d = d / H;
	    dis = dis / Diag;
	    W = zeros(size(d));
	    W(ind) = exp(-1*d(ind).^2 - dis(ind));
	    rho(:, k) = (a * U .^ m + b * T .^ eta) * W ./ den ;
	end
    
end

end


function w = w_thresh(data)
%W_THRESH Computes the threshold for d(lk) in d(lk)<w in calculation of rho

mean_data = mean(data, 1);
data_white = data - repmat(mean_data, size(data, 1), 1);
data_white = sqrt(sum(data_white .^ 2, 2));
d_bar = mean(data_white, 1);
delta = 0.3;
w = delta * d_bar;
end
