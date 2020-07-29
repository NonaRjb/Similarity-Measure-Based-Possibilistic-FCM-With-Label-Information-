clc 
close all 
clear variables

img = imread('Mri2.bmp');
img = imresize(img, 0.5);
[m, n] = size(img);
img_vec = reshape(im2double(img), [], 1);

mat = txt2mat('Aggregation.txt');

options = [NaN, NaN, 3, NaN, NaN, NaN, 100, NaN, 0, 0]; 
cluster_n = 7;
[U, T, obj_fcn] = sim_pfcm_l(mat(:, 1:2), cluster_n, options);
[~, i] = max(U, [], 1);
figure
gscatter(mat(:,1), mat(:,2), i')
% for k = 1 : cluster_n
%     figure
%     imshow(reshape(U(k,:), m, n))
% end
