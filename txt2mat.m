function mat = txt2mat(filename)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fun = @(s)regexp(s,'\t','split');
[fid,msg] = fopen(filename,'rt');
assert(fid>=3,msg)
% hdr = fun(fgetl(fid));
% val = fun(fgetl(fid));
str = fgetl(fid);
% while isempty(str)
% 	str = fgetl(fid);
% end
col = numel(fun(str));
out = {};
while ~feof(fid)
	%str = fgetl(fid);
	tmp = str2double(fun(str));
	tmp(end+1:col) = NaN;
	out{end+1} = tmp;
    str = fgetl(fid);
end
fclose(fid);
mat = vertcat(out{:});
end

