function H = h_bandwidth(data)
H = 0.01;
for k = 1 : size(data, 1)
    d = sqrt(sum((repmat(data(k), size(data, 1), 1)- data).^2, 2));
    max_d = max(d);
    if max_d > H
        H = max_d;
    end
%     for k = l : size(data, 1)
%         d = norm(data(l) - data(k));
%         if max_d < d
%             max_d = d;
%         end
%     end
end
end
