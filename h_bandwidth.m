function H = h_bandwidth(data)
H = 0.01;
if log2(size(data, 1)) <= 12
    d = sqrt(reshape(sum((repmat(data, size(data, 1), 1) - ...
        repelem(data, size(data, 1), 1)).^2, 2), size(data, 1), size(data, 1)));
    H = max(d(:));
    clear d
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
			max_d = max(d(:));
			clear d
		else
			d = sqrt(reshape(sum((repmat(data, n(i)-n(i-1), 1) - ...
			repelem(data(n(i-1)+1:n(i), :), size(data,1), 1)).^2, 2), size(data,1), n(i)-n(i-1)));
			if max(d(:)) > max_d
				max_d = max(d(:));
			end
			clear d
		end
    end
    H = max_d;
else
	for k = 1 : size(data, 1)
		d = sqrt(sum((repmat(data(k), size(data, 1), 1)- data).^2, 2));
		max_d = max(d);
		if max_d > H
			H = max_d;
		end
	end
end
end
