function H = h_bandwidth(data)
H = 0.01;
% D = repmat(data, size(data,1), 1);
% data_rpt = repelem(data, size(data,1), 1);
% d = sqrt(sum((D - data_rpt).^2, 2));
% H = max(d);
for k = 1 : size(data, 1)
    d = sqrt(sum((repmat(data(k), size(data, 1), 1)- data).^2, 2));
    max_d = max(d);
    if max_d > H
        H = max_d;
    end
end
end
