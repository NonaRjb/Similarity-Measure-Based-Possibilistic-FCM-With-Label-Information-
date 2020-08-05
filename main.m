clc 
close all 
clear variables

img = imread('Mri2.bmp');
img = imresize(img, 0.5);
[m, n] = size(img);
img_vec = reshape(im2double(img), [], 1);

mat = txt2mat('pathbased.txt');

options = [NaN, NaN, 3, 1, 0, NaN, 5, NaN, m, n]; 
cluster_n = 5;
[U, T, obj_fcn] = sim_pfcm_l(img_vec, cluster_n, options);
[~, i] = max(U, [], 1);
figure
imshow(reshape(i, m, n), [])
% figure
% gscatter(mat(:,1), mat(:,2), i')
% for k = 1 : cluster_n
%     figure
%     imshow(reshape(U(k,:), m, n))
% end
